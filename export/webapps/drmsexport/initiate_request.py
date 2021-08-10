#!/usr/bin/env python3


from argparse import Action as ArgsAction
from check_dbserver import StatusCode as CdbStatusCode
from parse_specification import StatusCode as PsStatusCode
from drms_export import Error as ExportError, ErrorCode as ExportErrorCode, Response, securedrms
from drms_parameters import DRMSParams, DPMissingParameterError
from drms_utils import Arguments as Args, ArgumentsError as ArgsError, Choices, CmdlParser, Formatter as DrmsLogFormatter, Log as DrmsLog, LogLevel as DrmsLogLevel, LogLevelAction as DrmsLogLevelAction, MakeObject, StatusCode as SC
from json import dumps as json_dumps, loads as json_loads
from os.path import join as path_join
from sys import exc_info as sys_exc_info, exit as sys_exit
from threading import Event, Timer

DEFAULT_LOG_FILE = 'ir_log.txt'

FILE_NAME_SIZE = 256

class StatusCode(SC):
    # fetch success status codes
    REQUEST_COMPLETE = 0, 'request has been completely processed'
    REQUEST_PROCESSING = 1, 'processing request'
    REQUEST_QUEUED = 2, 'request has been queued for processing'
    # this should go away - if in jsoc.export_new, but not in jsoc.export, then this is the same as REQUEST_QUEUED
    REQUEST_NOT_QUEUED = 6, 'request has not been queued yet'

    @classmethod
    def get_status_code(cls, request_status):
        if request_status == 6:

            return StatusCode.REQUEST_NOT_FOUND
        if request_status == 2  or request_status == 12:
            return StatusCode.REQUEST_QUEUED
        if request_status == 1:
            return StatusCode.REQUEST_PROCESSING
        if request_status == 0:
            return StatusCode.REQUEST_COMPLETE

class ErrorCode(ExportErrorCode):
    # fetch error codes
    REQUEST_TOO_LARGE = 3, 'requested payload is too large'
    REQUEST_FATAL_ERROR = 4, 'a fatal error occurred during processing'
    REQUEST_NOT_ONLINE = 5, 'exported files are no longer online'
    REQUEST_TOO_MANY = 7, 'too many simulatneous exports'

    # other IR errors
    PARAMETERS = 101, 'failure locating DRMS parameters'
    ARGUMENTS = 102, 'bad arguments'
    LOGGING = 103, 'failure logging messages'
    DRMS_CLIENT = 104, 'drms client error'
    STREAM_FORMAT = 105, 'format error in downloaded payload'
    EXPORT_ACTION = 106, 'failure calling export action'

    @classmethod
    def get_error_code(cls, request_error):
        if request_status == 7:
            return ErrorCode.REQUEST_REJECTED
        if request_status == 5:
            return ErrorCode.REQUEST_NOT_ONLINE
        if request_status == 4:
            return ErrorCode.REQUEST_FATAL
        if request_status == 3:
            return ErrorCode.REQUEST_TOO_LARGE

class IrBaseError(ExportError):
    def __init__(self, *, exc_info=None, msg=None):
        if exc_info is not None:
            self.exc_info = exc_info
            e_type, e_obj, e_tb = exc_info

            if msg is None:
                msg = f'{e_type.__name__}: {str(e_obj)}'
            else:
                msg = f'{msg} [ {e_type.__name__}: {str(e_obj)} ]'

        super().__init__(msg=msg)

class ParametersError(IrBaseError):
    _error_code = ErrorCode(ErrorCode.PARAMETERS)

    def __init__(self, *, exc_info=None, msg=None):
        super().__init__(exc_info=exc_info, msg=msg)

class ArgumentsError(IrBaseError):
    _error_code = ErrorCode(ErrorCode.ARGUMENTS)

    def __init__(self, *, exc_info=None, msg=None):
        super().__init__(exc_info=exc_info, msg=msg)

class LoggingError(IrBaseError):
    _error_code = ErrorCode(ErrorCode.LOGGING)

    def __init__(self, *, exc_info=None, msg=None):
        super().__init__(exc_info=exc_info, msg=msg)

class DRMSClientError(IrBaseError):
    _error_code = ErrorCode(ErrorCode.DRMS_CLIENT)

    def __init__(self, *, exc_info=None, msg=None):
        super().__init__(exc_info=exc_info, msg=msg)

class StreamFormatError(IrBaseError):
    _error_code = ErrorCode(ErrorCode.STREAM_FORMAT)

    def __init__(self, *, exc_info=None, msg=None):
        super().__init__(exc_info=exc_info, msg=msg)

class ExportActionError(IrBaseError):
    _error_code = ErrorCode(ErrorCode.EXPORT_ACTION)

    def __init__(self, *, exc_info=None, msg=None):
        super().__init__(exc_info=exc_info, msg=msg)

class ExportTypeAction(ArgsAction):
    def __call__(self, parser, namespace, value, option_string=None):
        if value == 'premium':
            setattr(namespace, 'package', None)
        elif value == 'mini':
            setattr(namespace, 'processing', None)
            setattr(namespace, 'package', None)
            setattr(namespace, 'access', 'http')
            setattr(namespace, 'file_format', 'native')
        else:
            setattr(namespace, 'access', 'streamed')

        setattr(namespace, self.dest, value)

class ExportArgumentsAction(ArgsAction):
    def __call__(self, parser, namespace, value, option_string=None):
        # convert json arguments to dict
        export_arguments = json_loads(value)
        if 'specification' not in export_arguments:
            raise ArgumentsError(f'missing required export argument `specification`')
        setattr(namespace, self.dest, export_arguments)

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
                raise ParametersError(msg=str(exc))

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
                parser.add_argument('db_host', help='the machine hosting the database that contains export requests from this site', metavar='<db host>', dest='db_host', required=True)
                parser.add_argument('export_type', help='the export type: premium, mini, or streamed', metavar='<webserver>', choices=Choices(['premium', 'mini', 'streamed']), action=ExportTypeAction, dest='export_type', required=True)
                parser.add_argument('webserver', help='the webserver invoking this script', metavar='<webserver>', action=create_webserver_action(drms_params), dest='webserver', required=True)
                parser.add_argument('address', help='the email addressed registered for export', metavar='<email address>', dest='address', required=True)
                parser.add_argument('arguments', help='export arguments', action=ExportArgumentsAction, dest='export_arguments', required=True)

                arguments = Arguments(parser=parser, args=args)
            else:
                # `program_args` has all `arguments` values, in final form; validate them
                def extract_module_args(*, db_host, export_type, webserver, address, export_arguments, drms_client_type='ssh', requestor=None, log_file=log_file, logging_level=logging_level, db_port=db_port, db_name=db_name, db_user=db_user):
                    arguments = {}

                    arguments['db_host'] = db_host
                    arguments['export_type'] = export_type
                    arguments['webserver'] = webserver
                    arguments['address'] = address
                    arguments['export_arguments'] = export_arguments
                    arguments['drms_client_type'] = drms_client_type
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
                raise ArgumentsError(msg=f'cannot specify private db server to handle public webserver requests')

            arguments.private_db_host = private_db_host
            arguments.drms_client = None

            cls._arguments = arguments

        return cls._arguments

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

def get_response_dict(request, status_or_error_code):
    fetch_status = request.status # 2/12 ==> asynchronous processing, 0 ==> synchronous processing
    error = request.error # None, unless an error occurred
    requestid = request.request_id # None, unless asynchronous processing
    count = request.file_count # None, unless synchronous processing
    rcount = request.record_count
    size = request.size
    method = request.method
    protocol = request.file_format

    if status_or_error_code == StatusCode.REQUEST_COMPLETE:
        data = json_dumps(list(zip(request.urls.record.to_list(), request.urls.url.to_list()))) # None, unless synchronous processing
    else:
        data = None

    response_dict = { 'status_code' : status_or_error_code, 'message' : status_or_error_code.description(), 'fetch_status' : fetch_status, 'fetch_error' : error, 'requestid' : requestid, 'count' : count, 'rcount' : rcount, 'size' : size, 'method' : method, 'protocol' : protocol, 'data' : data }

    return response_dict

def get_response(request):
    # was there an error?
    error_code = None

    try:
        error_code = ErrorCode(int(request.status))
    except KeyError:
        pass

    print(f'error code is {str(error_code)}')

    if error_code is not None:
        response_dict = get_response_dict(request, error_code)
        response = ErrorResponse.generate_response(**response_dict)
    else:
        try:
            status_code = StatusCode(int(request.status))
        except ValueError:
            raise DRMSClientError(exc_info=sys_exc_info(), msg=f'unexpected fetch status returned {str(request.status)}')

        response_dict = get_response_dict(request, status_code)
        response = Response.generate_response(**response_dict)

    return response

def stop_loop(do_loop):
    do_loop = False

def export_premium(*, drms_client, address, requestor=None, log, **export_arguments):
    '''
    export_arguments:
    {
      specification : <DRMS record-set specification>, # required
      processing : <dict of processing argument>, # optional
      package : <type and file name>, # required
      access : <http, ftp>, # required
      file_format <exported file type>, # optional
      file_format_args <exported file-type arguments>, # optional
      file_name_format : <format string for name of exported file>, # optional
      number_records : <maximum number of records exported> # optional
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

        request = drms_client.export(email=address, requestor=requestor, ds=export_arguments['specification'], processing=export_arguments['processing'], method=method, protocol=export_arguments['file_format'], protocol_args=export_arguments['file_format_args'] if 'file_format_args' in export_arguments else None, filename_fmt=export_arguments['file_name_format'] if 'file_name_format' in export_arguments else None, n=export_arguments['number_records'] if 'number_records' in export_arguments else None, synchronous_export=False)

        # if jsoc_fetch processed the request synchronously, then it will have returned status == 0 to `drms_client` - there is no request id created, so `request` will not have an `id` property; if jsoc_fetch processed the request asynchronously (because data were offline), then `request` will have an `id` property that contains the request ID needed by exportdata.html so it can poll for export-processing completion
        if 'id' in request:
            # jsoc_fetch executed the `url` branch of code, so the request is now asynchronous
            log.write_info(f'[ export_premium ] request {str(request.id)} is being processed asynchronously')

            # we want to wait for the request to appear in the request db table (jsoc.export); this should happen relatively quickly, so have a small timeout and if the timeout event occurs, return an error to the caller; otherwise, examine the retured status code - if it is QUEUED or PROCESSING, then make a 'pending' response; if it is COMPLETE, then make a 'complete' response; otherwise make an 'error' response
            do_loop = True
            timer = Timer(8, stop_loop, args=(do_loop,))
            while do_loop:
                request = drms_client.export_from_id(request.id) # updates status with jsoc_fetch exp_status call
                status_code = StatusCode.get_status_code(request.status)
                if status_code != StatusCode.REQUEST_NOT_FOUND:
                    break
        else:
            raise DRMSClientError(msg=f'unexpected `export_premium` response; missing request ID')
    except securedrms.SecureDRMSError as exc:
        raise DRMSClientError(exc_info=sys_exc_info())

    response = get_response(request)

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
        method='url-quick'

        request = drms_client.export(email=address, requestor=requestor, ds=export_arguments['specification'], method='url_quick', protocol='as-is', protocol_args=None, filename_fmt=export_arguments['file_name_format'] if 'file_name_format' in export_arguments else None, n=export_arguments['number_records'] if 'number_records' in export_arguments else None, synchronous_export=False)

        # if jsoc_fetch processed the request synchronously, then it will have returned status == 0 to `drms_client` - there is no request id created, so `request` will not have an `id` property; if jsoc_fetch processed the request asynchronously (because data were offline), then `request` will have an `id` property that contains the request ID needed by exportdata.html so it can poll for export-processing completion
        if request.request_id is not None:
            # jsoc_fetch executed the `url` branch of code, so the request is now asynchronous
            log.write_info(f'[ export_mini ] request {str(request.id)} is being processed asynchronously')

            # we want to wait for the request to appear in the request db table (jsoc.export); this should happen relatively quickly, so have a small timeout and if the timeout event occurs, return an error to the caller; otherwise, examine the retured status code - if it is QUEUED or PROCESSING, then make a 'pending' response; if it is COMPLETE, then make a 'complete' response; otherwise make an 'error' response

            # set timeout timer for 8 seconds
            do_loop = True
            timer = Timer(8, stop_loop, args=(do_loop,))
            while do_loop:
                request = drms_cccccclient.export_from_id(request.id) # updates status with jsoc_fetch exp_status call
                export_status_code = StatusCode.get_status_code(request.status)
                if export_status_code != StatusCode.REQUEST_NOT_FOUND:
                    break
        else:
            # jsoc_fetch processed the request synchronously; export files exist in a temporary SU now
            log.write_info([ f'[ export_mini ] request for `{export_arguments["specification"]}` was processed synchronously and is now complete' ])
    except securedrms.SecureDRMSError as exc:
        raise DRMSClientError(exc_info=sys_exc_info())
    except Exception as exc:
        raise DRMSClientError(exc_info=sys_exc_info())

    response = get_response(request)

    return response

def export_streamed(*, drms_client, address, requestor=None, log, **export_arguments):
    '''
    export_arguments:
    {
      specification : <DRMS record-set specification>, # required
      file_name_format : <format string for name of exported file> # optional
    }
    '''
    try:
        method = 'url_direct'

        download_stream = io.BytesIO()

        # spawn a thread to download and write to `download_stream`; will return while streaming is occurring;
        # drms-export-to-stdout will verify email address
        drms_client.export_package(spec=export_arguments['specification'], download_directory=None, filename_fmt=export_arguments['file_name_format'] if 'file_name_format' in export_arguments else None)

        file_name = ''

        # a FITS file prepended with filename; read the filename first
        while len(file_name) < FILE_NAME_SIZE:
            pipe_bytes = download_stream.read(FILE_NAME_SIZE - len(file_name))

            if len(pipe_bytes) == 0:
                raise StreamFormatError('improper formatting for FITS-file-name header')
            file_name = file_name + pipe_bytes.decode('UTF8')

        # truncate padding (0 bytes)
        file_name = file_name.rstrip('\x00')

        if len(file_name) == 0:
            # FITS file, which may have a filename, or it may not
            file_name = DUMMY_FITS_FILE_NAME

        # force alphanumeric (plus '_'), preserving the file extension
        matches = re.search(r'(^.+)[.]([^.]+$)', file_name)
        if matches is not None:
            base = matches.group(1)
            extension = matches.group(2)
            file_name = re.sub(r'[^a-zA-Z0-9_]', '_', base) + '.' + extension

        # we are streaming either a single FITS file; write HTTP header
        sys.stdout.buffer.write(b'Content-type: application/octet-stream\n')
        sys.stdout.buffer.write(b'Content-Disposition: attachment; filename="' + file_name.encode() + b'"\n')
        sys.stdout.buffer.write(b'Content-transfer-encoding: binary\n\n')

        # dump the remainder of the FITS file
        sys.stdout.buffer.write(download_stream.read())
        sys.stdout.buffer.flush()
    except securedrms.SecureDRMSError as exc:
        raise DRMSClientError(exc_info=sys_exc_info())

    response = Response.generate_response(status_code=StatusCode.SUCCESS)
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
                raise ParametersError(msg='unable to locate DRMS parameters package')

            arguments = Arguments.get_arguments(is_program=is_program, program_name=program_name, program_args=args, module_args=module_args, drms_params=drms_params)
        except ArgsError as exc:
            raise ArgumentsError(exc_info=sys_exc_info(), msg=f'{str(exc)}')
        except ExportError as exc:
            exc.exc_info = sys_exc_info()
            raise exc
        except Exception as exc:
            raise ArgumentsError(exc_info=sys_exc_info(), msg=f'{str(exc)}')

        try:
            formatter = DrmsLogFormatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
            log = DrmsLog(arguments.log_file, arguments.logging_level, formatter)
        except Exception as exc:
            raise LoggingError(msg=f'{str(exc)}')

        # parse specification to obtain series
        try:
            log.write_debug([ f'[ perform_action ] parsing record-set specification' ])
            action_type = 'parse_specification'
            action = Action.action(action_type=action_type, args={ 'specification' : arguments.export_arguments['specification'] })
            response = action() # __call__ returns dict

            print(f'response is {str(response)}')
            if response['status_code'] != PsStatusCode.SUCCESS:
                print(f'here OK')
                raise ExportActionError(msg=f'failure calling `{action_type}` action; status: `{response.status.description()}`')

            # `subsets` exists if status == PsStatusCode.SUCCESS
            subsets = response['subsets']
        except ExportError as exc:
            exc.exc_info = sys_exc_info()
            raise exc
        except Exception as exc:
            raise ExportActionError(exc_info=sys_exc_info(), msg=f'unable to parse record-set specification')

        try:
            if arguments.webserver.public:
                log.write_debug([ f'[ perform_action ] determining DB server' ])
                # need to determine if pass-through series have been specified; if so, use securedrms client that uses private db
                series_dict = { 'series' : [] }

                for subset in subsets:
                    series_dict['series'].append(subset['seriesname'])

                action_type = 'determine_db_server'
                action = Action.action(action_type=action_type, args={ 'public_dbhost' : arguments.db_host, 'series' : series_dict })
                response = action()
                if response['status_code'] != CdbStatusCode.SUCCESS:
                    raise ExportActionError(msg=f'failure calling `{action_type}` action; status: `{response.status.description()}`')
                if response['server'] is None:
                    raise ArgumentsError(msg=f'cannot service any series in `{", ".join(arguments.series)}`')
                db_host = response['server']
        except ExportError as exc:
            exc.exc_info = sys_exc_info()
            raise exc
        except Exception as exc:
            raise ExportActionError(exc_info=sys_exc_info(), msg=f'unable to determine supporting db server')

        if arguments.drms_client is None:
            try:
                use_public_db_host = True if db_host != arguments.private_db_host else False
                debug = True if arguments.logging_level == DrmsLogLevel.DEBUG else False
                factory = securedrms.SecureClientFactory(debug=debug, email=arguments.address)

                drms_server = 'jsoc_external' if use_public_db_host else 'jsoc_internal'
                use_ssh = True if arguments.drms_client_type == 'ssh' else False
                use_internal = False if use_public_db_host else True
                connection_info = { 'dbhost' : db_host, 'dbport' : arguments.db_port, 'dbname' : arguments.db_name, 'dbuser' : arguments.db_user }
                log.write_debug([ f'[ perform_action ] creating securedrms client' ])
                drms_client = factory.create_client(server=drms_server, use_ssh=use_ssh, use_internal=use_internal, connection_info=connection_info)
            except securedrms.SecureDRMSError as exc:
                raise DRMSClientError(exc_info=sys_exc_info())
            except Exception as exc:
                raise DRMSClientError(exc_info=sys_exc_info())
        else:
            drms_client = arguments.drms_client

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
                print('MMM')
                response = export_mini(drms_client=drms_client, address=arguments.address, requestor=arguments.requestor, log=log, **export_arguments)
                print('GGG')
            elif arguments.export_type == 'streamed':
                # supports no-processing, no-tar, stream-access, native file format attributes; exports a single file only, all SUs must be online; use securedrms.SecureClient.export_package()
                log.write_info([ f'[ perform_action ] servicing streamed request for user `{arguments.address}`: {str(export_arguments)}' ])
                response = export_streamed(drms_client=drms_client, address=arguments.address, requestor=arguments.requestor, log=log, **export_arguments)
        except ExportError as exc:
            if not hasattr(exc, 'exc_info'):
                print('MMM')
                exc.exc_info = sys_exc_info()
            print('UUU')
            raise exc
        except Exception as exc:
            print('VVV')
            raise ExportActionError(exc_info=sys_exc_info(), msg=f'{str(exc)}')
    except ExportError as exc:
        response = exc.response
        if log:
            e_type, e_obj, e_tb = exc.exc_info
            log.write_error([ f'ERROR LINE {str(e_tb.tb_lineno)}: {exc.message}' ])

    return response

# for use in export web app
from action import Action
class InitiateRequestAction(Action):
    actions = [ 'start_premium_export', 'start_mini_export', 'start_streamed_export' ]
    def __init__(self, *, method, webserver, dbhost=None, dbport=None, dbname=None, dbuser=None, address, requestor=None, export_arguments):
        self._method = getattr(self, method)
        self._webserver = webserver # dict
        self._address = adderss
        self._export_arguments = export_arguments
        self._options = {}
        self._options['requestor'] = requestor
        self._options['dbhost'] = dbhost
        self._options['dbport'] = dbport
        self._options['dbname'] = dbname
        self._options['dbuser'] = dbuser

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
        response = perform_action(is_program=False, export_type='premium', webserver=self._webserver, address=self._address, export_arguments=self._export_arguments, options=self._options)
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
        print('calling mini')
        response = perform_action(is_program=False, export_type='mini', webserver=self._webserver, address=self._address, export_arguments=self._export_arguments, options=self._options)
        print('mini returned')
        return response

    def start_streamed_export(self):
        # returns dict
        response = perform_action(is_program=False, export_type='streamed', webserver=self._webserver, address=self._address, export_arguments=self._export_arguments, options=self._options)
        return response


if __name__ == "__main__":
    try:
        response = perform_action(is_program=True)
        print('at the end')
        print(response.generate_json())
    except ExportError as exc:
        response = exc.response
    except Exception as exc:
        error = DRMSClientError(exc_info=sys_exc_info())
        response = error.response

    # Always return 0. If there was an error, an error code (the 'status' property) and message (the 'statusMsg' property) goes in the returned HTML.
    sys_exit(0)
else:
    pass
