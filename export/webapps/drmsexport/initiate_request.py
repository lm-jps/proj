#!/usr/bin/env python3

from argparse import Action as ArgsAction
import asyncio
from copy import deepcopy
from json import loads as json_loads, dumps as json_dumps
from os.path import join as path_join
import re
from socket import SHUT_RDWR
from sys import exc_info as sys_exc_info, exit as sys_exit, stdout as sys_stdout
from time import sleep
from threading import Timer
from urllib.parse import urlunsplit
from werkzeug.datastructures import Headers

from drms_export import Connection, create_generator, ExpServerBaseError, Error as ExportError, ErrorCode as ExportErrorCode, ErrorResponse, Response, get_arguments as ss_get_arguments, get_message, send_message
from drms_parameters import DRMSParams, DPMissingParameterError
from drms_utils import Arguments as Args, ArgumentsError as ArgsError, Choices, CmdlParser, Formatter as DrmsLogFormatter, Log as DrmsLog, LogLevel as DrmsLogLevel, LogLevelAction as DrmsLogLevelAction, MakeObject, StatusCode as SC
from utils import extract_program_and_module_args, get_db_host

DEFAULT_LOG_FILE = 'ir_log.txt'

FILE_NAME_SIZE = 256
DUMMY_FITS_FILE_NAME = 'data.FITS'

class StatusCode(SC):
    # fetch success status codes
    REQUEST_COMPLETE = (0, 'export has been completely processed')
    REQUEST_PROCESSING = (1, 'processing export')
    REQUEST_QUEUED = (2, 'export has been queued for processing') # status == 2 exists only in jsoc.export_new; jsoc_fetch op=exp_request creates a record in this table and sets the status value to 2 (or 12)
    # this should go away - if in jsoc.export_new, but not in jsoc.export, then this is the same as REQUEST_QUEUED
    REQUEST_NOT_QUEUED = (6, 'export has not been queued yet') # can no longer happen with jsoc_fetch op=exp_status call
    REQUEST_QUEUED_DEBUG = (12, 'export has been queued for processing')
    REQUEST_GENERATOR_READY = (20, 'export is now generator accessible; (`destination`, `generator`) provided with this response')

class ErrorCode(ExportErrorCode):
    # fetch error codes
    REQUEST_TOO_LARGE = (3, 'requested payload is too large') # can only happen with jsoc_fetch op=exp_request call
    REQUEST_FATAL_ERROR = (4, 'a fatal error occurred during processing')
    REQUEST_NOT_ONLINE = (5, 'exported files are no longer online')
    REQUEST_TOO_MANY = (7, 'too many simulatneous exports') # can only happen with jsoc_fetch op=exp_request call
    REQUEST_SEGMENT_OFFLINE = (25, 'one or more segment files are offline; contact the JSOC for assistance')

    # other IR errors
    PARAMETERS = (101, 'failure locating DRMS parameters')
    ARGUMENTS = (102, 'bad arguments')
    LOGGING = (103, 'failure logging messages')
    STREAM_FORMAT = (104, 'format error in downloaded payload')
    EXPORT_ACTION = (105, 'failure calling export action')
    RESPONSE = (106, 'unable to generate valid response')
    EXPORT_SERVER = (107, 'export-server communication error')
    UNHANDLED_EXCEPTION = (108, 'unhandled exception')

class IrBaseError(ExportError):
    def __init__(self, *, exc_info=None, error_message=None):
        if exc_info is not None:
            import traceback

            # for use with some exception handlers
            self.exc_info = exc_info
            e_type, e_obj, e_tb = exc_info
            file_info = traceback.extract_tb(e_tb)[0]
            file_name = file_info.filename if hasattr(file_info, 'filename') else ''
            line_number = str(file_info.lineno) if hasattr(file_info, 'lineno') else ''

            if error_message is None:
                error_message = f'{file_name}:{line_number}: {e_type.__name__}: {str(e_obj)}'
            else:
                error_message = f'{error_message} [ {file_name}:{line_number}: {e_type.__name__}: {str(e_obj)} ]'

        super().__init__(error_message=error_message)

class ParametersError(IrBaseError):
    _error_code = ErrorCode.PARAMETERS

class ArgumentsError(IrBaseError):
    _error_code = ErrorCode.ARGUMENTS

class LoggingError(IrBaseError):
    _error_code = ErrorCode.LOGGING

class StreamFormatError(IrBaseError):
    _error_code = ErrorCode.STREAM_FORMAT

class ExportActionError(IrBaseError):
    _error_code = ErrorCode.EXPORT_ACTION

class ResponseError(IrBaseError):
    _error_code = ErrorCode(ErrorCode.RESPONSE)

class ExportServerError(IrBaseError):
    _error_code = ErrorCode.EXPORT_SERVER

class UnhandledExceptionError(IrBaseError):
    _error_code = ErrorCode.UNHANDLED_EXCEPTION

class ExportArgumentsAction(ArgsAction):
    def __call__(self, parser, namespace, value, option_string=None):
        # convert json arguments to dict
        export_arguments = self.json_to_dict(value)
        if 'specification' not in export_arguments:
            raise ArgumentsError(error_message=f'missing required export argument `specification`')
        setattr(namespace, self.dest, export_arguments)

    @classmethod
    def json_to_dict(cls, json_text):
        return json_loads(json_text)

def name_to_ws_obj(name, drms_params):
    webserver_dict = {}
    webserver_dict['host'] = name
    if name is None or len(name) == 0 or name.strip().lower() == 'none':
        # assume public
        webserver_dict['public'] = True
    else:
        webserver_dict['public'] = True if name.lower() != drms_params.get_required('WEB_DOMAIN_PRIVATE') else False

    return MakeObject(name='webserver', data=webserver_dict)()

def webserver_action_constructor(self, drms_params):
    self._drms_params = drms_params

def webserver_action(self, parser, namespace, value, option_string=None):
    webserver_obj = name_to_ws_obj(value, self._drms_params)
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
    def get_arguments(cls, *, is_program, program_name=None, program_args=None, module_args=None, drms_params, refresh=True):
        if cls._arguments is None or refresh:
            try:
                private_db_host = drms_params.get_required('SERVER')
                db_port = int(drms_params.get_required('DRMSPGPORT'))
                db_name = drms_params.get_required('DBNAME')
                db_user = drms_params.get_required('WEB_DBUSER')
                download_web_domain = drms_params.get_required('WEB_DOMAIN_PUBLIC')
            except DPMissingParameterError as exc:
                raise ParametersError(exc_info=sys_exc_info(), error_message=str(exc))

            if is_program:
                try:
                    log_file = path_join(drms_params.get_required('EXPORT_LOG_DIR'), DEFAULT_LOG_FILE)
                except DPMissingParameterError as exc:
                    raise ParametersError(exc_info=sys_exc_info(), error_message=str(exc))

                args = None

                if program_args is not None and len(program_args) > 0:
                    args = program_args

                parser_args = { 'usage' : '%(prog)s address=<registered email address> export-type=<premium/mini/streamed> dbhost=<db host> arguments=<export arguments specific to export type> [ -l/--log-file=<path to log file> ] [ -L/--logging-level=<critical/error/warning/info/debug> ] [ -N/--dbname=<db name> ] [ -P/--dbport=<db port> ] [ -r/--requestor=<> ] [ -U/--dbuser=<db user>] [ -w/--webserver=<host> ]' }

                if program_name is not None and len(program_name) > 0:
                    parser_args['prog'] = program_name

                parser = CmdlParser(**parser_args)

                # required
                parser.add_argument('address', help='the email addressed registered for export', metavar='<email address>', dest='address', required=True)
                parser.add_argument('dbhost', help='the machine hosting the database that contains export requests from this site', metavar='<db host>', dest='db_host', required=True)
                parser.add_argument('export-type', help='the export type: premium, mini, or streamed', metavar='<export type>', choices=Choices(['premium', 'mini', 'streamed']), dest='export_type', required=True)
                parser.add_argument('arguments', help='export arguments', action=ExportArgumentsAction, dest='export_arguments', required=True)

                # optional
                parser.add_argument('-l', '--log-file', help='the path to the log file', metavar='<log file>', dest='log_file', default=log_file)
                parser.add_argument('-L', '--logging-level', help='the amount of logging to perform; in order of increasing verbosity: critical, error, warning, info, debug', metavar='<logging level>', dest='logging_level', action=DrmsLogLevelAction, default=DrmsLogLevel.ERROR)
                parser.add_argument('-N', '--dbname', help='the name of the database that contains export requests', metavar='<db name>', dest='db_name', default=db_name)
                parser.add_argument('-P', '--dbport', help='the port on the host machine that is accepting connections for the database', metavar='<db host port>', dest='db_port', type=int, default=db_port)
                parser.add_argument('-r', '--requestor', help='the name of the export user', metavar='<requestor>', dest='requestor', default=None)
                parser.add_argument('-U', '--dbuser', help='the name of the database user account', metavar='<db user>', dest='db_user', default=db_user)
                parser.add_argument('-w', '--webserver', help='the webserver invoking this script', metavar='<webserver>', action=create_webserver_action(drms_params), dest='webserver', default=name_to_ws_obj(None, drms_params))

                arguments = Arguments(parser=parser, args=args)
            else:
                # `program_args` has all `arguments` values, in final form; validate them
                def extract_module_args(*, address, export_type, db_host, export_arguments, file_specification=None, log=None, db_name=db_name, db_port=db_port, requestor=None, db_user=db_user, webserver=None):
                    arguments = {}

                    arguments['address'] = address
                    arguments['export_type'] = export_type
                    arguments['db_host'] = db_host

                    # replace the '*file*' specification with `file_specification`, if it is provided
                    export_arguments_dict = ExportArgumentsAction.json_to_dict(export_arguments)

                    if file_specification is not None:
                        export_arguments_dict['specification'] = file_specification

                    arguments['export_arguments'] = export_arguments_dict
                    arguments['db_name'] = db_name
                    arguments['db_port'] = db_port
                    arguments['requestor'] = requestor
                    arguments['db_user'] = db_user
                    arguments['webserver'] = name_to_ws_obj(webserver, drms_params) # sets webserver.public = True if webserver is None

                    InitiateRequestAction.set_log(log)

                    return arguments

                # dict
                module_args_dict = extract_module_args(**module_args)
                arguments = Arguments(parser=None, args=module_args_dict)

            # if the caller has specified a public webserver, make sure that the db specified is not the private one
            if arguments.webserver.public and arguments.db_host == private_db_host:
                raise ArgumentsError(error_message=f'cannot specify private db server to handle public webserver requests')

            arguments.private_db_host = private_db_host
            arguments.download_web_domain = download_web_domain

            cls.validate_export_arguments(arguments=arguments)
            cls._arguments = cls.make_objects(arguments=arguments)

        return cls._arguments

    @classmethod
    def validate_export_arguments(cls, *, arguments):
        if arguments.export_type == 'premium':
            required_args = ('access', 'package', 'specification',)
            optional_args = ('file-format', 'file-format-args', 'file-name-format', 'number-records', 'processing',)
        elif arguments.export_type == 'mini':
            required_args = ('specification',)
            optional_args = ('file-name-format', 'number-records',)
        elif arguments.export_type == 'streamed':
            required_args = ('specification',)
            optional_args = ('file-name-format', 'download-directory')

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

def get_response_dict(export_status_dict, log):
    log.write_debug([ f'[ get_response_dict ]' ])
    error_code, status_code = get_request_status(export_status_dict) # 2/12 ==> asynchronous processing, 0 ==> synchronous processing
    data = None
    export_directory = None
    keywords_file= None
    tar_file = None
    method = None
    request_url = None
    access = None
    scheme = None
    package = None
    file_format = None
    record_count = None
    file_count = None
    mb_exported = None
    contact = None
    fetch_error = None

    # request ID is not consistently provided in response JSON
    if 'requestid' in export_status_dict:
        # no error, and non-image file formats
        request_id = export_status_dict['requestid']
    elif 'reqid' in export_status_dict:
        # no error, image protocols
        request_id = export_status_dict['reqid']
    elif 'REQUEST_ID' in export_status_dict:
        # some errors, no request ID in response, so use the one stored by calling code
        request_id = export_status_dict['REQUEST_ID']
    else:
        request_id = None

    # on error, there is no request ID
    request_id = request_id if request_id and (len(request_id) > 0) else None

    if 'protocol' in export_status_dict:
        file_format = export_status_dict['protocol'].lower()
    elif 'FILE_FORMAT' in export_status_dict:
        file_format = export_status_dict['FILE_FORMAT'].lower()
    else:
        file_format = None

    file_count = export_status_dict.get('count')
    record_count = export_status_dict.get('rcount')
    mb_exported = export_status_dict.get('size')

    if error_code is not None:
        log.write_debug([ f'[ get_response_dict ] error calling fetch `{error_code.description()}` for request `{request_id}`'])
        # a fetch error occurred (`fetch_error` has the error message returned by fetch, `status_description` has the IR error message)
        fetch_error = export_status_dict.get('error') # None, unless an error occurred
        status_code = error_code
    else:
        # no fetch error occurred
        log.write_debug([ f'[ get_response_dict ] NO error calling fetch, status `{status_code.description()}` for request `{request_id}`'])

        if status_code == StatusCode.REQUEST_COMPLETE:
            log.write_debug([ f'[ get_response_dict ] exported data were generated for request `{request_id}`'])

            # not all of these attributes are provides in all responses, so use `get()`
            export_directory = export_status_dict.get('dir')
            keywords_file = export_status_dict.get('keywords')
            tar_file = export_status_dict.get('tarfile')
            method = export_status_dict.get('method')

            download_web_domain = export_status_dict.get('DOWNLOAD_WEB_DOMAIN')

            if download_web_domain is not None and export_directory is not None:
                # URL of the export SU
                request_url = urlunsplit(('https', download_web_domain, export_directory, None, None))

            if method is not None:
                try:
                    access = method[:3]
                except ValueError:
                    access = 'url'

            scheme = 'https' if access == 'url' else 'ftp'

            package = { 'type' : None if tar_file is None else 'tar', 'file_name' : None if tar_file is None else tar_file }

            file_information = export_status_dict.get('data')
            if file_information is not None:
                # none, unless synchronous export processing occurred
                if package['type'] == 'tar':
                    data = [ (record['record'], record['filename']) for record in file_information ]
                else:
                    # make URLs
                    url_information = { 'record' : [], 'filename' : [], 'url' : []} # ( 'record' : [ record_specs ], 'filename' : [ filenames ], 'url' : [ urls ])
                    file_information_resolved = None

                    if file_format in ['mpg', 'mp4']:
                        file_information_adjusted = deepcopy(file_information)

                        # look at the first record only (there is only one)
                        if file_information_adjusted[0]['record'].startswith('movie'):
                            file_information_adjusted[0]['record'] = None

                        file_information_resolved = file_information_adjusted
                    else:
                        file_information_resolved = file_information

                    if package['type'] == 'tar':
                        # record.filename is full path
                        for record in file_information_resolved:
                            url_information['record'].append(record['record'])
                            url_information['filename'].append(path_basename(record['filename']))
                            url_information['url'].append(urlunsplit((scheme, download_web_domain, record['filename'], None, None)))
                    else:
                        # record.filename is base file name
                        for record in file_information_resolved:
                            url_information['record'].append(record['record'])
                            url_information['filename'].append(record['filename'])
                            url_information['url'].append(urlunsplit((scheme, download_web_domain, path_join(export_directory, record['filename']), None, None)))

                    data = list(zip(url_information['record'], url_information['url']))

    response_dict = deepcopy(export_status_dict)

    if isinstance(status_code, ErrorCode):
        contact = export_status_dict.get('contact')
        response_dict.update({ 'error_code' : error_code, 'error_message' : fetch_error, 'contact' : contact })
    else:
        response_dict.update({ 'status_code' : status_code, 'fetch_error' : fetch_error, 'request_id' : request_id, 'file_format' : file_format, 'package' : package, 'access' : access, 'number_records' : record_count, 'number_files' : file_count, 'mb_exported' : mb_exported, 'sums_directory' : export_directory, 'keywords_file' : keywords_file, 'tar_file' : tar_file, 'request_url' : request_url, 'export_data' : data })

    log.write_debug([ f'[ get_response_dict] response dictionary `{str(response_dict)}`'])

    return (error_code is not None, response_dict)

def get_response(export_status_dict, log):
    # was there an error?
    error_code = None

    error_occurred, response_dict = get_response_dict(export_status_dict, log)

    if error_occurred:
        # an error occurred in jsoc_fetch; this is not the same thing as an error happening in MR so we have to manually
        # create an error response (the status_code will not be a MR error code, but the fetch_status will be a fetch error code)
        response = ErrorResponse.generate_response(**response_dict)
    else:
        response = Response.generate_response(**response_dict)

    return response

def get_request_status(export_status_dict):
    error_code = None
    status_code = None

    try:
        error_code = ErrorCode(int(export_status_dict['status']))
    except KeyError:
        pass

    if error_code is None:
        try:
            status_code = StatusCode(int(export_status_dict['status']))
        except KeyError:
            raise InvalidMessageError(exc_info=sys_exc_info(), error_message=f'unexpected fetch status returned {str(export_status_dict["status"])}')

    return (error_code, status_code)

def get_request_status_code(export_status_dict):
    error_code, status_code = get_request_status(export_status_dict)
    code = error_code if error_code is not None else status_code

    return code

def request_is_complete(export_status_dict):
    error_code, status_code = get_request_status(export_status_dict)
    return status_code == StatusCode.REQUEST_COMPLETE

def request_is_pending(export_status_dict):
    error_code, status_code = get_request_status(export_status_dict)
    return status_code == StatusCode.REQUEST_PROCESSING or status_code == StatusCode.REQUEST_QUEUED or status_code == StatusCode.REQUEST_NOT_QUEUED or status_code == StatusCode.REQUEST_QUEUED_DEBUG

def stop_loop(do_loop):
    do_loop = False

def export_premium(*, address, requestor=None, db_host, download_web_domain, log, export_arguments):
    '''
    export_arguments:
    {
        access : <http, ftp>, # required
        file-format <exported file type>, # optional
        file-format-args <exported file-type arguments>, # optional
        file-name-format : <format string for name of exported file>, # optional
        number-records : <maximum number of records exported>, # optional
        processing : <dict of processing argument>, # optional
        package : <type, file_name>, # required
        specification : <DRMS record-set specification> # required
    }
    '''
    response = None
    status_code = None

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

    try:
        # use socket server to call jsoc_fetch
        nested_arguments = ss_get_arguments(is_program=False, module_args={})

        with Connection(server=nested_arguments.server, listen_port=nested_arguments.listen_port, timeout=nested_arguments.message_timeout, log=log) as connection:
            try:
                skip_quit = False

                message = { 'request_type' : 'premium_export', 'address' : address, 'requestor' : requestor, 'specification' : export_arguments['specification'], 'processing' : export_arguments['processing'], 'method' : method, 'file_format' : export_arguments['file-format'], 'file_format_args' : export_arguments['file-format-args'], 'file_name_format' : export_arguments['file-name-format'], 'number_records' : export_arguments['number-records'], 'db_host' : db_host }

                response = send_request(message, connection, log)
                export_status_dict = json_loads(response)

                if export_status_dict.get('export_server_status') == 'export_server_error':
                    raise ExportServerError(error_message=f'{export_status_dict["error_message"]}')

                status_code = get_request_status_code(export_status_dict)

                # to make URLs for the downloadables
                export_status_dict['DOWNLOAD_WEB_DOMAIN'] = download_web_domain

                if request_is_complete(export_status_dict):
                    log.write_info([ f'[ export_premium ] request {str(export_status_dict["requestid"])} was processed synchronously and is complete (description `{status_code.description()}`)' ])
                elif request_is_pending(export_status_dict):
                    # jsoc_fetch executed the `url` branch of code, so the request is now asynchronous
                    log.write_info([ f'[ export_premium ] request {str(export_status_dict["requestid"])} is being processed asynchronously (description `{status_code.description()}`)' ])

                    # this loop is to ensure that the status code is not REQUEST_NOT_QUEUED; originally, this meant
                    # in jsoc.export_new, but not jsoc.export - but now jsoc_fetch will not return this code; jsoc_fetch
                    # will return REQUEST_PROCESSING if the export is in jsoc.export_new, but not jsoc.export
                    # ** so this loop is a no-op
                    do_loop = True
                    timer = Timer(8, stop_loop, args=(do_loop,))
                    timer.start()
                    while do_loop:
                        message = { 'request_type' : 'export_status', 'address' : address, 'request_id' : export_status_dict['requestid'], 'db_host' : db_host }

                        response = send_request(message, connection, log)
                        export_status_dict = json_loads(response)

                        if export_status_dict.get('export_server_status') == 'export_server_error':
                            raise ExportServerError(error_message=f'{export_status_dict["error_message"]}')

                        status_code = get_request_status_code(export_status_dict)

                        if isinstance(status_code, ErrorCode):
                            # GRRRRR! if an error occurred, then export_status_dict contains almost no info, not even the request ID; add
                            # property that does contain the request id value
                            export_status_dict[REQUEST_ID] = export_status_dict['requestid']

                        export_status_dict['DOWNLOAD_WEB_DOMAIN'] = download_web_domain

                        log.write_debug([ f'[ export_premium ] request {str(export_status_dict["requestid"])} status is `{status_code.description()}`' ])

                        if status_code != StatusCode.REQUEST_NOT_QUEUED:
                            timer.cancel()
                            break

                        sleep(1)
                else:
                    # must be an error - get_response() will pass the error code backt to the client
                    log.write_debug([ f'[ export_premium ] export request failed (description `{status_code.description()}`)' ])
            except ExportServerError:
                skip_quit = True
                raise
            finally:
                if not skip_quit:
                    message = { 'request_type' : 'quit' }
                    send_request(message, connection, log)
    except ExpServerBaseError as exc:
        raise ExportServerError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')
    except Exception as exc:
        raise ExportServerError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

    # if jsoc_fetch processed the request synchronously, then it will have returned status == 0 to the socket server - there is no request id created, so `request` will not have an `id` property; if jsoc_fetch processed the request asynchronously (because data were offline), then `request` will have an `id` property that contains the request ID needed by exportdata.html so it can poll for export-processing completion

    try:
        response = get_response(export_status_dict, log)
    except Exception as exc:
        raise ResponseError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

    return response

def export_mini(*, address, requestor=None, db_host, download_web_domain, log, export_arguments):
    '''
    export_arguments:
    {
      specification : <DRMS record-set specification>, # required
      file-name-format : <format string for name of exported file>, # optional
      number-records : <maximum number of records exported> # optional
    }
    '''
    response = None

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

    try:
        # use socket server to call jsoc_fetch
        nested_arguments = ss_get_arguments(is_program=False, module_args={})

        with Connection(server=nested_arguments.server, listen_port=nested_arguments.listen_port, timeout=nested_arguments.message_timeout, log=log) as connection:
            try:
                skip_quit = False

                message = { 'request_type' : 'mini_export', 'address' : address, 'requestor' : requestor, 'specification' : export_arguments['specification'], 'file_name_format' : export_arguments['file-name-format'], 'number_records' : export_arguments['number-records'], 'db_host' : db_host }

                response = send_request(message, connection, log)
                export_status_dict = json_loads(response)

                if export_status_dict.get('export_server_status') == 'export_server_error':
                    raise ExportServerError(error_message=f'{export_status_dict["error_message"]}')

                status_code = get_request_status_code(export_status_dict)

                # to make URLs for the downloadables
                export_status_dict['DOWNLOAD_WEB_DOMAIN'] = download_web_domain

                log.write_debug([ f'[ export_mini ] request {str(export_status_dict["requestid"])} status is `{status_code.description()}`' ])

                if request_is_complete(export_status_dict):
                    # if jsoc_fetch processed the request synchronously, then it will have returned status == 0 and
                    # there will be no request ID created
                    log.write_info([ f'[ export_mini ] request for {str(export_arguments["specification"])} was processed synchronously and is complete' ])
                elif request_is_pending(export_status_dict):
                    # jsoc_fetch executed the `url` branch of code, because data were offline, so the request is
                    # now asynchronous; if jsoc_fetch processed the request asynchronously, then it will have returned
                    # status == 2 and there will be a request ID created
                    log.write_info([ f'[ export_mini ] request {str(export_status_dict["requestid"])} is being processed asynchronously' ])

                    log.write_info([ f'[ export_mini ] data offline for request `{str(export_status_dict["requestid"])}`; must perform premium export - export request is being processed asynchronously' ])

                    # this loop is to ensure that the status code is not REQUEST_NOT_QUEUED; originally, this meant
                    # in jsoc.export_new, but not jsoc.export - but now jsoc_fetch will not return this code; jsoc_fetch
                    # will return REQUEST_PROCESSING if the export is in jsoc.export_new, but not jsoc.export
                    # ** so this loop is a no-op
                    do_loop = True
                    timer = Timer(8, stop_loop, args=(do_loop,))
                    timer.start()
                    while do_loop:
                        message = { 'request_type' : 'export_status', 'address' : address, 'request_id' : export_status_dict['requestid'], 'db_host' : db_host }

                        response = send_request(message, connection, log)
                        export_status_dict = json_loads(response)

                        if export_status_dict.get('export_server_status') == 'export_server_error':
                            raise ExportServerError(error_message=f'{export_status_dict["error_message"]}')

                        status_code = get_request_status_code(export_status_dict)

                        if isinstance(status_code, ErrorCode):
                            # GRRRRR! if an error occurred, then export_status_dict contains almost no info, not even the request ID; add
                            # property that does contain the request id value
                            export_status_dict[REQUEST_ID] = export_status_dict['requestid']

                        # to make URLs for the downloadables
                        export_status_dict['DOWNLOAD_WEB_DOMAIN'] = download_web_domain

                        log.write_debug([ f'[ export_premium ] request {str(export_status_dict["requestid"])} status is `{status_code.description()}`' ])

                        if status_code != StatusCode.REQUEST_NOT_QUEUED:
                            timer.cancel()
                            break

                        sleep(1)
                else:
                    # must be an error - get_response() will pass the error code backt to the client
                    pass
            except ExportServerError:
                skip_quit = True
                raise
            finally:
                if not skip_quit:
                    message = { 'request_type' : 'quit' }
                    send_request(message, connection, log)
    except ExpServerBaseError as exc:
        raise ExportServerError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')
    except Exception as exc:
        raise ExportServerError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

    try:
        response = get_response(export_status_dict, log)
    except Exception as exc:
        raise ResponseError(exc_info=sys_exc_info())

    return response

# extract name of file to be downloaded from socket
# returns tuple (<boolean>, <data>) where the first element indicates if there
# are more data to download
def read_from_connection(*, destination):
    if destination['has_header'] and not destination['header_extracted']:
        # this is the generator's first iteration
        # a FITS file prepended with filename; read the filename first
        number_bytes_read = 0
        file_name_buffer = []
        while True:
            bytes_read = destination['connection'].recv(FILE_NAME_SIZE) # blocking
            file_name_buffer.append(bytes_read)
            number_bytes_read += len(bytes_read)

            if number_bytes_read == FILE_NAME_SIZE:
                break

        # truncate padding (0 bytes)
        file_name = b''.join(file_name_buffer).rstrip(b'\x00').decode()

        # force alphanumeric (plus '_'), preserving the file extension
        matches = re.search(r'(^.+)[.]([^.]+$)', file_name)
        if matches is not None:
            base = matches.group(1)
            extension = matches.group(2)
            file_name = f'{re.sub(r"[^a-zA-Z0-9_]", "_", base)}.{extension}'

        destination['file_name'] = file_name # for use in `perform_action()` caller
        destination['header_extracted'] = True

        return (True, '')
    else:
        # dump the remainder of the FITS file (stdout is what the child process is dumping to); partial
        # reads are fine
        bytes_read = destination['connection'].recv(4096)
        if not bytes_read:
            # server done writing payload data (server shutdown connection) - close client end of connection
            destination['connection'].shutdown(SHUT_RDWR)
            destination['connection'].close()

            return (False, b'')
        else:
            return (True, bytes_read)

def export_streamed(*, request_action, address, db_host, log, export_arguments):
    '''
    export_arguments:
    {
      specification : <DRMS record-set specification>, # required
      file-name-format : <format string for name of exported file> # optional
      download-directory : <directory to which exported files are downloaded> # optional
    }
    '''
    response = None
    destination = None
    generator = None

    # new socket-server stuff

    # `header_extracted` - set to True after generator has read file name from child process stdout
    # `file_name` - the name to give the file that contains the downloaded content
    # `stream_reader` - the function called by the generator; it reads content from the child process' stdout
    # `has_header` - if True, then `stream_reader` will receive the header (file_name) BEFORE the generator
    #   is created
    destination = { 'header_extracted' : False, 'file_name' : None, 'reader' : read_from_connection, 'has_header' : True }

    connection = None
    close_connection = False
    try:
        # use socket server to call jsoc_fetch
        nested_arguments = ss_get_arguments(is_program=False, module_args={})

        # over connection...
        # 1. client sends 'streamed_export' message
        # 2. server sends response telling client to start downloading payload over socket
        # 3. client reads payload from socket until server shuts down server end of socket
        # 4. client closes client end of socket

        # call connect() when instantiating Connection outside of a context manager; the CM returns
        # a socket object, but the constructor does not
        connection = Connection(server=nested_arguments.server, listen_port=nested_arguments.listen_port, timeout=nested_arguments.message_timeout, log=log).connect()
        destination['connection'] = connection

        # this message will start the child process on the export server; the child process will write to connection
        message = { 'request_type' : 'streamed_export', 'address' : address, 'specification' : export_arguments['specification'], 'file_name_format' : export_arguments['file-name-format'], 'db_host' : db_host }

        # response is pretty much empty; as long as it indicates success, the child process was started on the export server
        response = send_request(message, connection, log)
        export_status_dict = json_loads(response)

        if export_status_dict.get('export_server_status') == 'export_server_error':
            raise ExportServerError(error_message=f'{export_status_dict["error_message"]}')

        # no error, so it is ok for client to start downloading payload, which happens through generator

        # create_generator() makes a generator object that can be used to download the payload:
        #   - reads header from connection and stores in destination
        #   - returns generator object to download the rest of the payload over the connection
        generator = create_generator(destination=destination)

        if request_action is None:
            download_directory = export_arguments.get('download-directory')
            file_name = destination.get('file_name')

            if download_directory is not None and file_name is not None:
                destination['file_path'] = path_join(download_directory, file_name)
        else:
            # if this is being run in a module context, then the caller will iterate through generator;
            # the caller also needs access to the destination (which contains the file name for the
            # downloaded content) to generate HTTP headers
            request_action.destination = destination
            request_action.generator = generator

        response = Response.generate_response(status_code=StatusCode.REQUEST_GENERATOR_READY)
    except ExportServerError as exc:
        close_connection = True
        response = ErrorResponse.generate_response(error_code=ErrorCode.REQUEST_FATAL_ERROR, error_message=f'{exc.message}')
    except Exception as exc:
        close_connection = True
        response = ErrorResponse.generate_response(error_code=ErrorCode.REQUEST_FATAL_ERROR, error_message=f'{str(exc)}')

    if close_connection:
        if connection is not None:
            connection.shutdown(SHUT_RDWR)
            connection.close()

    else:
        # the socket connection needs to remain open until generation is complete; so the code that
        # iterates through the generator must close the connection (which is accessible through the
        # destination returned)
        pass

    # this never gets back to browser - browser just displays the file data
    # streamed back to it; the browser reads the headers sent back (which
    # should include a 200 if everything worked out)
    return (response, (destination, generator))

def send_request(request, connection, log):
    json_message = json_dumps(request)
    send_message(connection, json_message)
    message = get_message(connection)

    return message

def perform_action(*, action_obj, is_program, program_name=None, **kwargs):
    # catch all expections so we can always generate a response
    response = None
    destination = None
    generator = None
    log = None

    try:
        program_args, module_args = extract_program_and_module_args(is_program=is_program, **kwargs)

        try:
            drms_params = DRMSParams()

            if drms_params is None:
                raise ParametersError(error_message='unable to locate DRMS parameters package')

            arguments = Arguments.get_arguments(is_program=is_program, program_name=program_name, program_args=program_args, module_args=module_args, drms_params=drms_params)
        except Exception as exc:
            raise ArgumentsError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

        if is_program:
            try:
                formatter = DrmsLogFormatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
                log = DrmsLog(arguments.log_file, arguments.logging_level, formatter)
                InitiateRequestAction._log = log
            except Exception as exc:
                raise LoggingError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')
        else:
            log = action_obj.log

        if is_program:
            log.write_debug([ f'[ perform_action ] program invocation' ])
        else:
            log.write_debug([ f'[ perform_action ] module invocation' ])

        log.write_debug([ f'[ perform_action ] action arguments: {str(arguments)}' ])

        # parse specification to obtain series so we can check for pass-through series (relevant only if the user is on a public webserver);
        # parsing checks syntax, it does not check for the existence of series in the specification, so either a public or private drms client can be used
        if not InitiateRequestAction.is_valid_specification(arguments.export_arguments['specification'], arguments.db_host, arguments.webserver.host):
            raise ArgumentsError(error_message=f'invalid record-set specification')

        parsed_specification = InitiateRequestAction.get_parsed_specification(arguments.export_arguments['specification'], arguments.db_host, arguments.webserver.host)

        if not parsed_specification.attributes.hasfilts and arguments.number_records is None:
            raise ArgumentsError(error_message=f'must specify either a record-set filter, or a maximum number of records')

        series = []
        for subset in parsed_specification.attributes.subsets:
            series.append(subset.seriesname)

        resolved_db_host = get_db_host(webserver=arguments.webserver, series=series, private_db_host=arguments.private_db_host, db_host=arguments.db_host, db_port=arguments.db_port, db_name=arguments.db_name, db_user=arguments.db_user, exc=ExportActionError, log=log)

        if resolved_db_host is None:
            raise ArgumentsError(error_message=f'unable to determine db host suitable for series {", ".join(series)}')

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
                response = export_premium(address=arguments.address, requestor=arguments.requestor, db_host=resolved_db_host, download_web_domain=arguments.download_web_domain, log=log, export_arguments=export_arguments)
            elif arguments.export_type == 'mini':
                # supports no processing, no tar, http access, native file format attributes; use securedrms.SecureClient.export(method='url_quick')
                log.write_info([ f'[ perform_action ] servicing mini request for user `{arguments.address}`: {str(export_arguments)}' ])
                response = export_mini(address=arguments.address, requestor=arguments.requestor, db_host=resolved_db_host, download_web_domain=arguments.download_web_domain, log=log, export_arguments=export_arguments)
            elif arguments.export_type == 'streamed':
                # supports no-processing, no-tar, stream-access, native file format attributes; exports a single file only, all SUs must be online; use securedrms.SecureClient.export_package()
                log.write_info([ f'[ perform_action ] servicing streamed request for user `{arguments.address}`: {str(export_arguments)}' ])
                (response, (destination, generator)) = export_streamed(request_action=action_obj, address=arguments.address, db_host=resolved_db_host, log=log, export_arguments=export_arguments)
        except Exception as exc:
            raise ExportActionError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')
    except IrBaseError as exc:
        response = exc.response
        error_message = exc.message

        if log:
            log.write_error([ error_message ])
        elif is_program:
            print(error_message)
    except Exception as exc:
        response = UnhandledExceptionError(exc_info=sys_exc_info(), error_message=f'{str(exc)}').response
        error_message = str(exc)

        if log:
            log.write_error([ error_message ])
        elif is_program:
            print(error_message)

    if log is not None:
        log.write_info([ f'[ perform_action ] request complete; status {response.attributes.drms_export_status_description}' ])

    return (response, (destination, generator))

# for use in export web app
from action import Action
from parse_specification import ParseSpecificationAction

class InitiateRequestAction(Action):
    actions = [ 'start_premium_export', 'start_mini_export', 'start_streamed_export' ]

    _log = None

    def __init__(self, *, method, address, db_host, export_arguments, file_specification=None, log=None, db_name=None, db_port=None, requestor=None, db_user=None, webserver=None):
        self._method = getattr(self, method)
        self._address = address
        self._db_host = db_host
        self._export_arguments = export_arguments # json text
        self._options = {}
        self._options['file_specification'] = file_specification
        self._options['log'] = log
        self._options['db_name'] = db_name
        self._options['db_port'] = db_port
        self._options['requestor'] = requestor
        self._options['db_user'] = db_user
        self._options['webserver'] = webserver # host name - gets converted to object in `get_arguments()`

        self.destination = None
        self.generator = None

    def start_premium_export(self):
        '''
        export_arguments:
        {
          specification : <DRMS record-set specification>,
          processing : <dict of processing argument>,
          package : <type and file name>,
          access : <http, ftp>,
          file-format <exported file type>,
          file-name-format : <format string for name of exported file>,
          size-ratio : <multiplier to obtain image size from original image size>,
          number-records : <maximum number of records exported>
        }
        '''
        response = perform_action(action_obj=self, is_program=False, export_type='premium', address=self._address, db_host=self._db_host, export_arguments=self._export_arguments, options=self._options)
        return response

    def start_mini_export(self):
        '''
        export_arguments:
        {
          specification : <DRMS record-set specification>,
          file-name-format : <format string for name of exported file>,
          number-records : <maximum number of records exported>
        }
        '''
        response = perform_action(action_obj=self, is_program=False, export_type='mini', address=self._address, db_host=self._db_host, export_arguments=self._export_arguments, options=self._options)
        return response

    def start_streamed_export(self):
        '''
        export_arguments:
        {
          specification : <DRMS record-set specification>,
          file-name-format : <format string for name of exported file>
          download-directory : <directory to which exported files are downloaded> # optional
        }
        '''
        response = perform_action(action_obj=self, is_program=False, export_type='streamed', address=self._address, db_host=self._db_host, export_arguments=self._export_arguments, options=self._options)
        return response

    def generate_headers(self):
        if self.destination is not None and self.destination['header_extracted'] and self.destination['file_name'] is not None:
            headers = Headers()
            headers.add('Content-Disposition', 'attachment', filename=f'{self.destination["file_name"]}')
            headers.add('Content-Transfer-Encoding', 'BINARY')
            return headers

    @classmethod
    def is_valid_arguments(cls, arguments_json, log):
        cls.set_log(log)
        is_valid = None
        try:
            json_loads(arguments_json)
            is_valid = True
        except:
            is_valid = False

        return is_valid

    @classmethod
    def get_parsed_specification(cls, specification, db_host, webserver):
        action = Action.action(action_type='parse_specification', args={ 'log' : cls._log, 'specification' : specification, 'db_host' : db_host, 'webserver' : webserver })
        parsed_specification = action()

        if isinstance(parsed_specification, ErrorResponse) or parsed_specification is None:
            cls._log.write_error([ f'[ get_parsed_specification ] {parsed_specification.attributes.error_message}'])
            raise ArgumentsError(error_message=f'unable to parse specification `{specification}`')

        return parsed_specification

    @classmethod
    def is_valid_specification(cls, specification, db_host, webserver):
        is_valid = None

        try:
            if db_host is None:
                # if this method is called before URL arguments are parsed, then `db_host` is not known;
                # use default (public) export DB (specification parser does not use DB so it does not matter)
                try:
                    drms_params = DRMSParams()

                    if drms_params is None:
                        raise ParametersError(error_message='unable to locate DRMS parameters package')

                    db_host_resolved = drms_params.get_required('EXPORT_DB_HOST_DEFAULT')
                except DPMissingParameterError as exc:
                    raise ParametersError(exc_info=sys_exc_info(), error_message=str(exc))
            else:
                db_host_resolved = db_host

            # parse specification
            response = cls.get_parsed_specification(specification, db_host_resolved, webserver)

            is_valid = False if isinstance(response, ErrorResponse) else True
        except:
            is_valid = False

        return is_valid

    @property
    def log(self):
        return self.__class__._log

    @log.setter
    def log(self, log):
        self.__class__._log = log

    @classmethod
    def set_log(cls, log=None):
        if cls._log is None:
            cls._log = DrmsLog(None, None, None) if log is None else log

    @classmethod
    def get_log(cls):
        return cls._log

if __name__ == "__main__":
    response, (destination, generator) = perform_action(action_obj=None, is_program=True)

    if not isinstance(response, ErrorCode):
        if response.status_code == StatusCode.REQUEST_GENERATOR_READY:
            if destination is not None:
                file_path = destination.get('file_path')

            if generator is not None:
                # iterate through generator, which generates binary export content (e.g., a FITS file);
                # write output to stdout; the caller can save a file or redirect stdout to a file

                if file_path is not None:
                    with open(file_path, mode='wb') as file_out:
                        for data in generator:
                            file_out.write(data)
                else:
                    for data in generator:
                        sys_stdout.buffer.write(data)
        else:
            print(response.generate_json())

    # Always return 0. If there was an error, an error code (the 'status' property) and message (the 'statusMsg' property) goes in the returned HTML.
    sys_exit(0)
else:
    pass
