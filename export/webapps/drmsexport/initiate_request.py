#!/usr/bin/env python3

from argparse import Action as ArgsAction
import asyncio
from check_dbserver import StatusCode as CdbStatusCode
from parse_specification import StatusCode as PsStatusCode
from drms_export import Error as ExportError, ErrorCode as ExportErrorCode, ErrorResponse, Response, securedrms
from drms_parameters import DRMSParams, DPMissingParameterError
from drms_utils import Arguments as Args, ArgumentsError as ArgsError, Choices, CmdlParser, Formatter as DrmsLogFormatter, Log as DrmsLog, LogLevel as DrmsLogLevel, LogLevelAction as DrmsLogLevelAction, MakeObject, StatusCode as SC
from json import dumps as json_dumps, loads as json_loads
from os.path import join as path_join
import re
from sys import exc_info as sys_exc_info, exit as sys_exit, stdout as sys_stdout
from time import sleep
from threading import Event, Timer

DEFAULT_LOG_FILE = 'ir_log.txt'

FILE_NAME_SIZE = 256
DUMMY_FITS_FILE_NAME = 'data.FITS'

class StatusCode(SC):
    # fetch success status codes
    REQUEST_COMPLETE = (0, 'request has been completely processed')
    REQUEST_PROCESSING = (1, 'processing request')
    REQUEST_QUEUED = (2, 'request has been queued for processing') # status == 2 exists only in jsoc.export_new; jsoc_fetch op=exp_request creates a record in this table and sets the status value to 2 (or 12)
    # this should go away - if in jsoc.export_new, but not in jsoc.export, then this is the same as REQUEST_QUEUED
    REQUEST_NOT_QUEUED = (6, 'request has not been queued yet') # can no longer happen with jsoc_fetch op=exp_status call
    REQUEST_QUEUED_DEBUG = (12, 'request has been queued for processing')

class ErrorCode(ExportErrorCode):
    # fetch error codes
    REQUEST_TOO_LARGE = (3, 'requested payload is too large') # can only happen with jsoc_fetch op=exp_request call
    REQUEST_FATAL_ERROR = (4, 'a fatal error occurred during processing')
    REQUEST_NOT_ONLINE = (5, 'exported files are no longer online') # as far as I can tell, the export code never sets status to 5
    REQUEST_TOO_MANY = (7, 'too many simulatneous exports') # can only happen with jsoc_fetch op=exp_request call

    # other IR errors
    PARAMETERS = (101, 'failure locating DRMS parameters')
    ARGUMENTS = (102, 'bad arguments')
    LOGGING = (103, 'failure logging messages')
    DRMS_CLIENT = (104, 'drms client error')
    STREAM_FORMAT = (105, 'format error in downloaded payload')
    EXPORT_ACTION = (106, 'failure calling export action')
    RESPONSE = (107, 'unable to generate valid response')

class IrBaseError(ExportError):
    def __init__(self, *, exc_info=None, error_message=None):
        if exc_info is not None:
            self.exc_info = exc_info
            e_type, e_obj, e_tb = exc_info

            if error_message is None:
                error_message = f'{e_type.__name__}: {str(e_obj)}'
            else:
                error_message = f'{error_message} [ {e_type.__name__}: {str(e_obj)} ]'

        super().__init__(error_message=error_message)

class ParametersError(IrBaseError):
    _error_code = ErrorCode.PARAMETERS

class ArgumentsError(IrBaseError):
    _error_code = ErrorCode.ARGUMENTS

class LoggingError(IrBaseError):
    _error_code = ErrorCode.LOGGING

class DRMSClientError(IrBaseError):
    _error_code = ErrorCode.DRMS_CLIENT

class StreamFormatError(IrBaseError):
    _error_code = ErrorCode.STREAM_FORMAT

class ExportActionError(IrBaseError):
    _error_code = ErrorCode.EXPORT_ACTION

class ExportArgumentsAction(ArgsAction):
    def __call__(self, parser, namespace, value, option_string=None):
        # convert json arguments to dict
        export_arguments = json_loads(value)
        if 'specification' not in export_arguments:
            raise ArgumentsError(error_message=f'missing required export argument `specification`')
        setattr(namespace, self.dest, export_arguments)

class ResponseError(IrBaseError):
    _error_code = ErrorCode(ErrorCode.RESPONSE)

def webserver_action_constructor(self, drms_params):
    self._drms_params = drms_params

def webserver_action(self, parser, namespace, value, option_string=None):
    webserver_dict = {}
    webserver_dict['host'] = value
    webserver_dict['public'] = True if value.lower() != self._drms_params.get_required('WEB_DOMAIN_PRIVATE') else False
    webserver_obj = MakeObject(name='webserver', data=webserver_dict)()
    setattr(namespace, self.dest, webserver_obj)

def create_webserver_action(drms_params):
    cls = type('WebserverAction', (ArgsAction,),
    {
        '_drms_params' : drms_params,
        '__call__' : webserver_action,
    })

    return cls

class Arguments(Args):
    _arguments = None

    @classmethod
    def get_arguments(cls, *, is_program, program_name=None, program_args=None, module_args=None, drms_params):
        if cls._arguments is None:
            try:
                log_file = path_join(drms_params.get_required('EXPORT_LOG_DIR'), DEFAULT_LOG_FILE)
                private_db_host = drms_params.get_required('SERVER')
                db_port = int(drms_params.get_required('DRMSPGPORT'))
                db_name = drms_params.get_required('DBNAME')
                db_user = drms_params.get_required('WEB_DBUSER')
            except DPMissingParameterError as exc:
                raise ParametersError(exc_info=sys_exc_info(), error_message=str(exc))

            if is_program:
                args = None

                if program_args is not None and len(program_args) > 0:
                    args = program_args

                if program_name is not None and len(program_name) > 0:
                    parser = CmdlParser(prog=program_name, usage='%(prog)s [ --dbhost=<db host> ] [ --dbport=<db port> ] [ --dbname=<db name> ] [ --dbuser=<db user>] [ --skiptar=<boolean value> ] [ --compression=<boolean value ] [ --file_name_format=<> ] webserver=<> spec=<>')
                else:
                    parser = CmdlParser(usage='%(prog)s [ --dbhost=<db host> ] [ --dbport=<db port> ] [ --dbname=<db name> ] [ --dbuser=<db user>] [ --skiptar=<boolean value> ] [ --compression=<boolean value ] [ --file_name_format=<> ] webserver=<> spec=<>')

                # optional
                parser.add_argument('-c', '--drms-client-type', help='securedrms client type (ssh, http)', choices=[ 'ssh', 'http' ], dest='drms_client_type', default='ssh')
                parser.add_argument('-l', '--log_file', help='the path to the log file', metavar='<log file>', dest='log_file', default=log_file)
                parser.add_argument('-L', '--logging_level', help='the amount of logging to perform; in order of increasing verbosity: critical, error, warning, info, debug', metavar='<logging level>', dest='logging_level', action=DrmsLogLevelAction, default=DrmsLogLevel.ERROR)
                parser.add_argument('-N', '--dbname', help='the name of the database that contains export requests', metavar='<db name>', dest='db_name', default=db_name)
                parser.add_argument('-P', '--dbport', help='the port on the host machine that is accepting connections for the database', metavar='<db host port>', dest='db_port', type=int, default=db_port)
                parser.add_argument('-r', '--requestor', help='the name of the export user', metavar='<requestor>', dest='requestor', default=None)
                parser.add_argument('-U', '--dbuser', help='the name of the database user account', metavar='<db user>', dest='db_user', default=db_user)

                # required
                parser.add_argument('address', help='the email addressed registered for export', metavar='<email address>', dest='address', required=True)
                parser.add_argument('export_type', help='the export type: premium, mini, or streamed', metavar='<webserver>', choices=Choices(['premium', 'mini', 'streamed']), dest='export_type', required=True)
                parser.add_argument('db_host', help='the machine hosting the database that contains export requests from this site', metavar='<db host>', dest='db_host', required=True)
                parser.add_argument('webserver', help='the webserver invoking this script', metavar='<webserver>', action=create_webserver_action(drms_params), dest='webserver', required=True)
                parser.add_argument('arguments', help='export arguments', action=ExportArgumentsAction, dest='export_arguments', required=True)

                arguments = Arguments(parser=parser, args=args)
                arguments.drms_client = None
            else:
                # `program_args` has all `arguments` values, in final form; validate them
                def extract_module_args(*, address, export_type, db_host, webserver, export_arguments, drms_client_type='ssh', drms_client=None, requestor=None, log_file=log_file, logging_level=DrmsLogLevel.ERROR, db_port=db_port, db_name=db_name, db_user=db_user):
                    arguments = {}

                    arguments['address'] = address
                    arguments['export_type'] = export_type
                    arguments['db_host'] = db_host
                    arguments['webserver'] = webserver
                    arguments['export_arguments'] = export_arguments
                    arguments['drms_client_type'] = drms_client_type
                    arguments['drms_client'] = drms_client
                    arguments['requestor'] = requestor
                    arguments['log_file'] = log_file
                    arguments['logging_level'] = logging_level
                    arguments['db_port'] = db_port
                    arguments['db_name'] = db_name
                    arguments['db_user'] = db_user

                    return arguments

                # dict
                module_args_dict = extract_module_args(**module_args)
                arguments = Arguments(parser=None, args=module_args_dict)

            # if the caller has specified a public webserver, make sure that the db specified is not the private one
            if arguments.webserver.public and arguments.db_host == private_db_host:
                raise ArgumentsError(error_message=f'cannot specify private db server to handle public webserver requests')

            arguments.private_db_host = private_db_host

            cls.validate_export_arguments(arguments=arguments)
            cls._arguments = cls.make_objects(arguments=arguments)

        return cls._arguments

    @classmethod
    def validate_export_arguments(cls, *, arguments):
        if arguments.export_type == 'premium':
            required_args = ('access', 'package', 'specification',)
            optional_args = ('file_format', 'file_format_args', 'file_name_format', 'number_records', 'processing',)
        elif arguments.export_type == 'mini':
            required_args = ('specification',)
            optional_args = ('file_name_format', 'number_records',)
        elif arguments.export_type == 'streamed':
            required_args = ('specification',)
            optional_args = ('file_name_format',)

        all_args = []
        all_args.extend(required_args)
        all_args.extend(optional_args)

        for required in required_args:
            if not arguments.export_arguments.get(required, False):
                raise ArgumentsError(error_message=f'missing required export argument `{required}`')

        for optional in optional_args:
            # set all missing optional argument values to None
            if not arguments.export_arguments.get(optional, False):
                arguments.export_arguments[optional] = None

        for arg in arguments.export_arguments.keys():
            if arg not in all_args:
                raise ArgumentsError(error_message=f'unexpected export argument `{arg}`')

    @classmethod
    def make_objects(cls, *, arguments):
        if 'package' in arguments.export_arguments:
            data = arguments.export_arguments['package']
            if not isinstance(data, dict):
                raise ArgumentsError(error_message=f'the `package` export arugment must be a `dict`')
            if 'type' not in data:
                raise ArgumentsError(error_message=f'`package` dict argument must have a `type` property')
            arguments.export_arguments['package'] = MakeObject(name='package', data=data)()

        return arguments

def jsonify(data_frame):
    # input:
    #         record                    filename
    #   0   hmi.M_720s[2017.1.23]       hmi__m_720s_20170123.fits
    #   1   hmi.M_720s[2017.1.24]       hmi__m_720s_20170124.fits
    #   2   hmi.M_720s[2017.1.25]       hmi__m_720s_20170125.fits
    #
    # output:
    # [
    #    {
    #      "record" : "hmi.M_720s[2017.1.23]",
    #      "filename" : "hmi__m_720s_20170123.fits"
    #    },
    #    {
    #      "record" : "hmi.M_720s[2017.1.24]",
    #      "filename" : "hmi__m_720s_20170124.fits"
    #    },
    #    {
    #      "record" : "hmi.M_720s[2017.1.25]",
    #      "filename" : "hmi__m_720s_20170125.fits"
    #    }
    # ]
    return data_frame.to_json(orient='records')

def get_response_dict(export_request, log):
    log.write_debug([ f'[ get_response_dict ]' ])
    error_code, status_code = get_request_status(export_request) # 2/12 ==> asynchronous processing, 0 ==> synchronous processing
    data = None
    export_directory = None
    keywords_file= None
    tar_file = None
    request_url = None

    fetch_error = None

    # None, unless asynchronous processing (but None if error)
    request_id = export_request.request_id if export_request.request_id is not None else export_request.REQUEST_ID

    if error_code is not None:
        log.write_debug([ f'[ get_response_dict ] error calling fetch `{error_code.description()}` for request `{request_id}`'])
        # a fetch error occurred (`fetch_error` has the error message returned by fetch, `status_description` has the IR error message)
        fetch_error = export_request.error # None, unless an error occurred
        status_code = error_code
    else:
        # no fetch error occurred
        log.write_debug([ f'[ get_response_dict ] NO error calling fetch, status `{status_code.description()}` for request `{request_id}`'])

        if status_code == StatusCode.REQUEST_COMPLETE:
            log.write_debug([ f'[ get_response_dict ] exported data were generated for request `{request_id}`'])
            data = json_dumps(list(zip(export_request.urls.record.to_list(), export_request.urls.url.to_list()))) # None, unless synchronous processing
            export_directory = export_request.dir
            keywords_file = export_request.keywords
            tar_file = export_request.tarfile
            request_url = export_request.request_url


    count = export_request.file_count # None, unless synchronous processing
    rcount = export_request.record_count
    size = export_request.size
    method = export_request.method
    protocol = export_request.file_format

    if isinstance(status_code, ErrorCode):
        response_dict = { 'error_code' : error_code, 'error_message' : fetch_error }
    else:
        response_dict = { 'status_code' : status_code, 'fetch_status' : int(export_request.status), 'fetch_error' : fetch_error, 'request_id' : request_id, 'count' : count, 'rcount' : rcount, 'size' : size, 'method' : method, 'protocol' : protocol, 'data' : data, 'export_directory' : export_directory, 'keywords_file' : keywords_file, 'tar_file' : tar_file, 'request_url' : request_url }

    log.write_debug([ f'[ get_response_dict] response dictionary `{str(response_dict)}`'])

    return (error_code is not None, response_dict)

def get_response(export_request, log):
    # was there an error?
    error_code = None

    error_occurred, response_dict = get_response_dict(export_request, log)

    if error_occurred:
        # an error occurred in jsoc_fetch; this is not the same thing as an error happening in MR so we have to manually
        # create an error response (the status_code will not be a MR error code, but the fetch_status will be a fetch error code)
        response = ErrorResponse.generate_response(**response_dict)
    else:
        response = Response.generate_response(**response_dict)

    return response

def get_request_status(export_request):
    error_code = None
    status_code = None

    try:
        error_code = ErrorCode(int(export_request.status))
    except KeyError:
        pass

    if error_code is None:
        try:
            status_code = StatusCode(int(export_request.status))
        except KeyError:
            raise DRMSClientError(exc_info=sys_exc_info(), error_message=f'unexpected fetch status returned {str(export_request.status)}')

    return (error_code, status_code)

def get_request_status_code(export_request):
    error_code, status_code = get_request_status(export_request)
    code = error_code if error_code is not None else status_code

    return code

def request_is_complete(export_request):
    error_code, status_code = get_request_status(export_request)
    return status_code == StatusCode.REQUEST_COMPLETE

def request_is_pending(export_request):
    error_code, status_code = get_request_status(export_request)
    return status_code == StatusCode.REQUEST_PROCESSING or status_code == StatusCode.REQUEST_QUEUED or status_code == StatusCode.REQUEST_NOT_QUEUED or status_code == REQUEST_QUEUED_DEBUG

def stop_loop(do_loop):
    do_loop = False

def export_premium(*, drms_client, address, requestor=None, log, **export_arguments):
    '''
    export_arguments:
    {
        access : <http, ftp>, # required
        file_format <exported file type>, # optional
        file_format_args <exported file-type arguments>, # optional
        file_name_format : <format string for name of exported file>, # optional
        number_records : <maximum number of records exported>, # optional
        processing : <dict of processing argument>, # optional
        package : <type, package_file_name>, # required
        specification : <DRMS record-set specification> # required
    }
    '''
    response = None

    try:
        if export_arguments['package'].type == 'tar':
            if export_arguments['access'] == 'ftp':
                method = 'ftp-tar'
            else:
                method = 'url-tar'
        else:
            if export_arguments['access'] == 'ftp':
                method = 'ftp'
            else:
                method = 'url'

        export_request = drms_client.export(email=address, requestor=requestor, ds=export_arguments['specification'], process=export_arguments['processing'], method=method, protocol=export_arguments['file_format'], protocol_args=export_arguments['file_format_args'], filename_fmt=export_arguments['file_name_format'], n=export_arguments['number_records'], synchronous_export=False)
        status_code = get_request_status_code(export_request)

        if isinstance(status_code, ErrorCode):
            # GRRRRR! if an error occurred, then export_request contains almost no info, not even the request ID; add
            # property that does contain the request id value
            export_request.REQUEST_ID = request_id

        log.write_debug([ f'[ export_premium ] request {str(export_request.request_id)} status is `{status_code.description()}`' ])

        # if jsoc_fetch processed the request synchronously, then it will have returned status == 0 to `drms_client` - there is no request id created, so `request` will not have an `id` property; if jsoc_fetch processed the request asynchronously (because data were offline), then `request` will have an `id` property that contains the request ID needed by exportdata.html so it can poll for export-processing completion
        if request_is_complete(export_request):
            log.write_info([ f'[ export_premium ] request {str(export_request.request_id)} was processed synchronously and is complete' ])
        elif request_is_pending(export_request):
            # jsoc_fetch executed the `url` branch of code, so the request is now asynchronous
            log.write_info([ f'[ export_premium ] request {str(export_request.request_id)} is being processed asynchronously' ])

            # we want to wait for the request to appear in the request db table (jsoc.export) - this should no longer happen since jsoc_fetch now will return status == 2 if the request is in jsoc.export_new, but not yet in jsoc.export
            do_loop = True
            timer = Timer(8, stop_loop, args=(do_loop,))
            timer.start()
            while do_loop:
                export_request = drms_client.export_from_id(export_request.request_id) # updates status with jsoc_fetch exp_status call
                status_code = get_request_status_code(export_request)

                if isinstance(status_code, ErrorCode):
                    # GRRRRR! if an error occurred, then export_request contains almost no info, not even the request ID; add
                    # property that does contain the request id value
                    export_request.REQUEST_ID = request_id

                log.write_debug([ f'[ export_premium ] request {str(export_request.request_id)} status is `{status_code.description()}`' ])
                if status_code != StatusCode.REQUEST_NOT_QUEUED:
                    timer.cancel()
                    break

                sleep(1)
        else:
            # must be an error
            raise DRMSClientError(error_message=f'failure processing premium export ({status_code.description()})')
    except securedrms.SecureDRMSError as exc:
        raise DRMSClientError(exc_info=sys_exc_info())
    except Exception as exc:
        raise DRMSClientError(exc_info=sys_exc_info())

    try:
        response = get_response(export_request, log)
    except Exception as exc:
        raise ResponseError(exc_info=sys_exc_info())

    return response

def export_mini(*, drms_client, address, requestor=None, log, **export_arguments):
    '''
    export_arguments:
    {
      specification : <DRMS record-set specification>, # required
      file_name_format : <format string for name of exported file>, # optional
      number_records : <maximum number of records exported> # optional
    }
    '''
    response = None

    try:
        # export() runs the export job asynchronously, returning a SecureExportRequest instance; the following attributes are returned and used by the export web page to display a table of links to the exported data files:
        # ExportRequest.urls - a pandas.DataFrame with three columns; the third column has the urls needed by exportdata.html:
        #   >>> req.urls
        #                                                 record          filename                                                url
        #   0  hmi.M_720s[2017.10.10_00:00:00_TAI][3]{magneto...  magnetogram.fits  http://jsoc.stanford.edu/SUM91/D980961841/S000...
        #   1  hmi.M_720s[2017.10.10_00:12:00_TAI][3]{magneto...  magnetogram.fits  http://jsoc.stanford.edu/SUM91/D980961841/S000...
        #
        #   ExportRequest.id - is None for `url-quick` IF the requested data are online; if any are offline, then jsoc_fetch executes the `url` branch of code, but keeps the method as `url_quick`

        # jsonify(request.urls) returns:
        # >>> req.urls.to_json()
        # { "record" :
        #    { "0" : "hmi.M_720s[2017.10.10_00:00:00_TAI][3]{magnetogram}",
        #      "1" : "hmi.M_720s[2017.10.10_00:12:00_TAI][3]{magnetogram}"
        #    },
        #   "filename" :
        #    { "0" : "magnetogram.fits",
        #      "1" : "magnetogram.fits"},
        #    "url" :
        #    { "0" : "http:\\/\\/jsoc.stanford.edu\\/SUM91\\/D980961841\\/S00003\\/magnetogram.fits",
        #      "1" : "http:\\/\\/jsoc.stanford.edu\\/SUM91\\/D980961841\\/S00004\\/magnetogram.fits"
        #    }
        # }
        method = 'url-quick'

        export_request = drms_client.export(email=address, requestor=requestor, ds=export_arguments['specification'], method='url_quick', protocol='as-is', protocol_args=None, filename_fmt=export_arguments['file_name_format'], n=export_arguments['number_records'], synchronous_export=False)
        status_code = get_request_status_code(export_request)

        if isinstance(status_code, ErrorCode):
            export_request.REQUEST_ID = None

        log.write_debug([ f'[ export_mini ] request export status for user {address} is `{status_code.description()}`' ])

        # if jsoc_fetch processed the request synchronously, then it will have returned status == 0 to `drms_client` - there is no request id created, so `request` will not have an `id` property; if jsoc_fetch processed the request asynchronously (because data were offline), then `request` will have an `id` property that contains the request ID needed by exportdata.html so it can poll for export-processing completion
        if request_is_complete(export_request):
            log.write_info([ f'[ export_mini ] request for `{export_arguments["specification"]}` was processed synchronously and is now complete' ])
        elif request_is_pending(export_request):
            # jsoc_fetch executed the `url` branch of code without error, so the request is now asynchronous, and request_id property exists
            log.write_info([ f'[ export_mini ] data offline for request `{str(export_request.request_id)}`; must perform premium export - export request is being processed asynchronously' ])

            # we want to wait for the request to appear in the request db table (jsoc.export); this should happen relatively quickly, so have a small timeout and if the timeout event occurs, return an error to the caller; otherwise, examine the retured status code - if it is QUEUED or PROCESSING, then make a 'pending' response; if it is COMPLETE, then make a 'complete' response; otherwise make an 'error' response

            # set timeout timer for 8 seconds
            do_loop = True
            timer = Timer(8, stop_loop, args=(do_loop,))
            timer.start()
            while do_loop:
                export_request = drms_client.export_from_id(export_request.request_id) # updates status with jsoc_fetch exp_status call
                status_code = get_request_status_code(export_request)

                log.write_debug([ f'[ export_mini ] request `{str(export_request.request_id)}` status is `{status_code.description()}`' ])
                if status_code != StatusCode.REQUEST_NOT_QUEUED:
                    timer.cancel()
                    break

                sleep(1)
        else:
            # must be an error
            raise DRMSClientError(error_message=f'failure processing mini export ({status_code.description()})')
    except securedrms.SecureDRMSError as exc:
        raise DRMSClientError(exc_info=sys_exc_info())
    except Exception as exc:
        raise DRMSClientError(exc_info=sys_exc_info())

    try:
        response = get_response(export_request, log)
    except Exception as exc:
        raise ResponseError(exc_info=sys_exc_info())

    return response

async def read_from_proc(destination, proc):
    if not destination.file_name_read:
        # a FITS file prepended with filename; read the filename first
        bytes_read = await proc.stdout.read(FILE_NAME_SIZE)

        # truncate padding (0 bytes)
        file_name = bytes_read.rstrip(b'\x00').decode()

        # force alphanumeric (plus '_'), preserving the file extension
        matches = re.search(r'(^.+)[.]([^.]+$)', file_name)
        if matches is not None:
            base = matches.group(1)
            extension = matches.group(2)
            file_name = f'{re.sub(r"[^a-zA-Z0-9_]", "_", base)}.{extension}'

        disposition = f'Content-Disposition: attachment; filename={file_name}\n'

        # return data back to securedrms to write it to stdout
        data = b'Content-type: application/octet-stream\n' + disposition.encode() + b'Content-transfer-encoding: binary\n\n'

        destination.file_name_read = True
    else:
        # dump the remainder of the FITS file
        bytes_read = await proc.stdout.read(4096)
        data = bytes_read

    return data

def export_streamed(*, drms_client, address, requestor=None, log, **export_arguments):
    '''
    export_arguments:
    {
      specification : <DRMS record-set specification>, # required
      file_name_format : <format string for name of exported file> # optional
    }
    '''
    response = None

    try:
        method = 'url_direct'

        # spawn a thread to download and write to `download_stream`; will return while streaming is occurring;
        # drms-export-to-stdout will verify email address
        export_request = drms_client.export_and_stream_file(spec=export_arguments['specification'], filename_fmt=export_arguments['file_name_format'])

        # asyncio.run(stream_file_data(request, download_stream))

        destination = securedrms.OnTheFlyDownloader.StreamDestination(sys_stdout.buffer)
        destination.file_name_read = False
        destination.stream_reader = read_from_proc
        export_request.download(destination=destination)
    except securedrms.SecureDRMSError as exc:
        raise DRMSClientError(exc_info=sys_exc_info())

    response = Response.generate_response(status_code=StatusCode.REQUEST_COMPLETE)
    return response

def perform_action(is_program, program_name=None, **kwargs):
    # catch all expections so we can always generate a response
    response = None
    log = None

    try:
        args = None
        module_args = None

        if is_program:
            args = []
            for key, val in kwargs.items():
                if val is not None:
                    if key == 'options':
                        for option, option_val in val.items():
                            args.append(f'--{option}={option_val}')
                    else:
                        args.append(f'{key}={val}')
        else:
            # a loaded module
            module_args = {}
            for key, val in kwargs.items():
                if val is not None:
                    if key == 'options':
                        for option, option_val in val.items():
                            module_args[option] = option_val
                    else:
                        module_args[key] = val

        try:
            drms_params = DRMSParams()

            if drms_params is None:
                raise ParametersError(error_message='unable to locate DRMS parameters package')

            arguments = Arguments.get_arguments(is_program=is_program, program_name=program_name, program_args=args, module_args=module_args, drms_params=drms_params)
        except ArgsError as exc:
            raise ArgumentsError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')
        except Exception as exc:
            raise ArgumentsError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

        try:
            formatter = DrmsLogFormatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
            log = DrmsLog(arguments.log_file, arguments.logging_level, formatter)
        except Exception as exc:
            raise LoggingError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

        if is_program:
            log.write_debug([ f'[ perform_action ] program invocation' ])
        else:
            log.write_debug([ f'[ perform_action ] module invocation' ])

        log.write_debug([ f'[ perform_action ] action arguments: {str(arguments)}' ])

        # parse specification to obtain series so we can check for pass-through series (relevant only if the user is on a public webserver);
        # parsing checks syntax, it does not check for the existence of series in the specification, so either a public or private drms client can be used
        factory = None
        public_drms_client = None
        private_drms_client = None
        if arguments.webserver.public:
            # arguments.db_host must be a public db server; use external drms client
            log.write_debug([ f'[ perform_action ] public webserver {arguments.webserver.host} initiating request' ])
            try:
                if arguments.drms_client is None:
                    log.write_debug([ f'[ perform_action ] no securedrms client provided; creating public one' ])
                    debug = True if arguments.logging_level == DrmsLogLevel.DEBUG else False
                    factory = securedrms.SecureClientFactory(debug=debug, email=arguments.address)
                    use_ssh = True if arguments.drms_client_type == 'ssh' else False
                    connection_info = { 'dbhost' : arguments.db_host, 'dbport' : arguments.db_port, 'dbname' : arguments.db_name, 'dbuser' : arguments.db_user }

                    public_drms_client = factory.create_client(server='jsoc_external', use_ssh=use_ssh, use_internal=False, connection_info=connection_info)
                else:
                    # public client (since the webserver is public)
                    log.write_debug([ f'[ perform_action ] public securedrms client provided' ])
                    public_drms_client = arguments.drms_client

                log.write_debug([ f'[ perform_action ] parsing record-set specification' ])
                response_dict = public_drms_client.parse_spec(arguments.export_arguments['specification'])

                if response_dict['errMsg'] is not None:
                    raise DRMSClientError(error_message=f'failure parsing record-set specification {arguments.export_arguments["specification"]}: {response_dict["errMsg"]}')

                # `subsets` exists if status == PsStatusCode.SUCCESS
                subsets = response_dict['subsets']
            except Exception as exc:
                raise ExportActionError(exc_info=sys_exc_info(), error_message=f'unable to parse record-set specification')

            try:
                # need to determine if pass-through series have been specified; if so, use securedrms client that uses private db
                series_dict = { 'series' : [] }

                for subset in subsets:
                    series_dict['series'].append(subset['seriesname'])

                log.write_debug([ f'[ perform_action ] determining DB server suitable for requested data from series {",".join(series_dict["series"])}' ])

                action_type = 'determine_db_server'
                action_args = { 'public_dbhost' : arguments.db_host, 'series' : series_dict, 'drms_client' : public_drms_client }
                action = Action.action(action_type=action_type, args=action_args)
                response = action()
                if response['status_code'] != CdbStatusCode.SUCCESS:
                    raise ExportActionError(error_message=f'failure calling `{action_type}` action; status: `{response.status.description()}`')
                if response['server'] is None:
                    raise ArgumentsError(error_message=f'cannot service any series in `{", ".join(arguments.series)}`')
                db_host = response['server']
            except Exception as exc:
                raise ExportActionError(exc_info=sys_exc_info(), error_message=f'unable to determine supporting db server')
        else:
            log.write_debug([ f'[ perform_action ] private webserver {arguments.webserver.host} initiating request' ])
            private_drms_client = arguments.drms_client # could be None

        # now we know whether we should be using a public or private drms client
        use_public_db_host = True if arguments.webserver.public and db_host != arguments.private_db_host else False
        if use_public_db_host:
            log.write_debug([ f'[ perform_action ] public db host will be used to service request' ])
            if public_drms_client is not None:
                drms_client = public_drms_client
            else:
                # we must have had a public client - we used it to parse the record-set specification
                raise DRMSClientError(error_message=f'missing public drms client')
        else:
            log.write_debug([ f'[ perform_action ] private db host must be used to service request' ])
            # private drms client needed; if a public one was passed in, then private_drms_client is None
            if private_drms_client is not None:
                log.write_debug([ f'[ perform_action ] no securedrms client provided; creating private one' ])
                drms_client = private_drms_client
            else:
                try:
                    debug = True if arguments.logging_level == DrmsLogLevel.DEBUG else False
                    if factory is None:
                        factory = securedrms.SecureClientFactory(debug=debug, email=arguments.address)
                    use_ssh = True if arguments.drms_client_type == 'ssh' else False
                    connection_info = { 'dbhost' : arguments.private_db_host, 'dbport' : arguments.db_port, 'dbname' : arguments.db_name, 'dbuser' : arguments.db_user }
                    log.write_debug([ f'[ perform_action ] creating private securedrms client' ])
                    private_drms_client = factory.create_client(server='jsoc_internal', use_ssh=use_ssh, use_internal=False, connection_info=connection_info)
                    drms_client = private_drms_client
                except securedrms.SecureDRMSError as exc:
                    raise DRMSClientError(exc_info=sys_exc_info())
                except Exception as exc:
                    raise DRMSClientError(exc_info=sys_exc_info())

        # there are four potential attributes an export request can have:
        #   processing (processed ==> original image is modified); NA for keyword-only exports
        #   package (tar ==> exported files are packaged into a tar file)
        #   access (ftp ==> exported produce is accessible from ftp server; http ==> exported produce is accessible from http server; stream ==> data are streamed back to user)
        #   file format (NATIVE, FITS, JPEG, MPEG, MP4); NA for keyword-only exports
        # depending on the operation selected, some of these attributes may be determined and immutable; also, some file-format conversions are not possible (for example, if an export contains metadata, but no image data, then none of the image file formats are relevant)

        # possible export arguments (from webapp):
        #   specification [RecordSet] (record-set specification)
        #   address [Notify] (registered email address)
        #   requestor [Requestor] export user's name, either entered in text edit, or from cookie
        #   export_type [Method - url, url_direct, url_quick, ftp, url-tar, ftp-tar] these map to...
        #     url ==> export_premium (no tar)
        #     url_direct ==> export_streamed
        #     url_quick ==> export_mini
        #     ftp ==> export_premium (ftp-access)
        #     url-tar ==> export_premium (tar)
        #     ftp-tar ==> export_premium (tar, ftp access)
        #   processing_info [Processing] (processing-program-specific arguments)
        #   tar [from Method (name contains '-tar' string)]
        #   access [from Method]
        #   file_format [Protocol - as-is, FITS, JPEG, MPEG, MP4] (exported file type)
        #   file_name_format [Filename Format]
        #   size_ratio [not shown] determined by code in processing.js plus keyword metadata

        # the response to each of these calls contains these attributes:
        #   count : the number of files exported
        #   rcount : the number of DRMS records in the export request
        #   size : the size, in MB, of the file payload
        #   method (request.method) : the export method (`url`, `url_quick`, ...)
        #   protocol : the exported file format (`as-is`, `fits`, `jpeg`, ...)
        #   status : export status (1, 2, 4, 6, 7)
        #   error : a string describing a jsoc_fetch error or null
        #   requestid : the export Request ID or null (for synchronous requests)

        try:
            export_arguments = arguments.export_arguments # dict

            if arguments.export_type == 'premium':
                # supports all export four export attributes; use securedrms.SecureClient.export(method='url')
                log.write_info([ f'[ perform_action ] servicing premium request for user `{arguments.address}`: {str(export_arguments)}' ])
                response = export_premium(drms_client=drms_client, address=arguments.address, requestor=arguments.requestor, log=log, **export_arguments)
            elif arguments.export_type == 'mini':
                # supports no processing, no tar, http access, native file format attributes; use securedrms.SecureClient.export(method='url_quick')
                log.write_info([ f'[ perform_action ] servicing mini request for user `{arguments.address}`: {str(export_arguments)}' ])
                response = export_mini(drms_client=drms_client, address=arguments.address, requestor=arguments.requestor, log=log, **export_arguments)
            elif arguments.export_type == 'streamed':
                # supports no-processing, no-tar, stream-access, native file format attributes; exports a single file only, all SUs must be online; use securedrms.SecureClient.export_package()
                log.write_info([ f'[ perform_action ] servicing streamed request for user `{arguments.address}`: {str(export_arguments)}' ])
                response = export_streamed(drms_client=drms_client, address=arguments.address, requestor=arguments.requestor, log=log, **export_arguments)
        except Exception as exc:
            raise ExportActionError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')
    except ExportError as exc:
        response = exc.response
        if log:
            e_type, e_obj, e_tb = exc.exc_info
            log.write_error([ f'[ perform_action ] ERROR LINE {str(e_tb.tb_lineno)}: {exc.message}' ])

    log.write_info([f'[ perform_action ] request complete; status {str(response)}'])
    return response

# for use in export web app
from action import Action
class InitiateRequestAction(Action):
    actions = [ 'start_premium_export', 'start_mini_export', 'start_streamed_export' ]
    def __init__(self, *, method, address, db_host, webserver, export_arguments, drms_client_type=None, drms_client=None, requestor=None, db_port=None, db_name=None, db_user=None):
        self._method = getattr(self, method)
        self._address = address
        self._db_host = db_host
        self._webserver = webserver # dict
        self._export_arguments = export_arguments
        self._options = {}
        self._options['drms_client_type'] = drms_client_type
        self._options['drms_client'] = drms_client
        self._options['requestor'] = requestor
        self._options['db_port'] = db_port
        self._options['db_name'] = db_name
        self._options['db_user'] = db_user

    def start_premium_export(self):
        '''
        export_arguments:
        {
          specification : <DRMS record-set specification>,
          processing : <dict of processing argument>,
          package : <type and file name>,
          access : <http, ftp>,
          file_format <exported file type>,
          file_name_format : <format string for name of exported file>,
          size_ratio : <multiplier to obtain image size from original image size>,
          number_records : <maximum number of records exported>
        }
        '''
        response = perform_action(is_program=False, export_type='premium', address=self._address, db_host=self._db_host, webserver=self._webserver, export_arguments=self._export_arguments, options=self._options)
        return response

    def start_mini_export(self):
        '''
        export_arguments:
        {
          specification : <DRMS record-set specification>,
          file_name_format : <format string for name of exported file>,
          number_records : <maximum number of records exported>
        }
        '''
        response = perform_action(is_program=False, export_type='mini', address=self._address, db_host=self._db_host, webserver=self._webserver, export_arguments=self._export_arguments, options=self._options)
        return response

    def start_streamed_export(self):
        # returns dict
        response = perform_action(is_program=False, export_type='streamed', address=self._address, db_host=self._db_host, webserver=self._webserver, export_arguments=self._export_arguments, options=self._options)
        return response


if __name__ == "__main__":
    try:
        response = perform_action(is_program=True)
    except ExportError as exc:
        response = exc.response
    except Exception as exc:
        error = ExportActionError(exc_info=sys_exc_info())
        response = error.response

    print(response.generate_json())

    # Always return 0. If there was an error, an error code (the 'status' property) and message (the 'statusMsg' property) goes in the returned HTML.
    sys_exit(0)
else:
    pass
