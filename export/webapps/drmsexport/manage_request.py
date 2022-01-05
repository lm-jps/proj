#!/usr/bin/env python3

# THIS IS A RENAMED VERSION OF manage-request.py. It turns out using an '-' for a python file name was a bad decision; a module name cannot have an hyphen in it, so it is not possible to import one that does; AND CVS does not permit file name changes, so, manage-request.py will be used for the legacy exportdata.html, and manage_request.py will be used for the new flask app

# The arguments to this script are parsed by cgi.FieldStorage(), which knows how to parse
# both HTTP GET and POST requests. A nice feature is that we can test the script as it runs in a CGI context
# by simply running on the command line with a single argument that is equivalent to an HTTP GET parameter string
# (e.g., address=gimli@mithril.com&addresstab=jsoc.export_addresses&domaintab=jsoc.export_addressdomains).

# Parameters:
#   address (required) - The email address to check or register.
#   addresstab (required) - The database table containing all registered (or registration-pending) email addresses.
#   domaintab (required) - The database table containing all email domains.
#   dbuser (optional) - The database account to be used when connecting to the database. The default is the value of the WEB_DBUSER parameter in DRMSParams.
#   checkonly (optional) - If set to 1, then no attept is made to register an unregistered email. In this case, if no error occurs then the possible return status codes are RV_REGISTEREDADDRESS, RV_REGISTRATIONPENDING, or RV_UNREGISTEREDADDRESS. The default is False (unknown addresses are registered).

# Testing:
#   to test providing CGI arguments, set `TEST_CGI` to True, and then provide a single QUERY_STRING-like argument, like:
#   manage-request.py 'address=person@domain&operation=check'

from argparse import Action as ArgsAction
from copy import deepcopy
from datetime import timedelta
from json import dumps as json_dumps, loads as json_loads
from os.path import join as path_join
import psycopg2
from re import compile as re_compile
from sys import exc_info as sys_exc_info, exit as sys_exit, stdout as sys_stdout
from urllib.parse import urlunsplit

from drms_export import Connection, ExpServerBaseError, Error as ExportError, ErrorCode as ExportErrorCode, ErrorResponse, Response, get_arguments as ss_get_arguments, get_message, send_message
from drms_parameters import DRMSParams, DPMissingParameterError
from drms_utils import Arguments as Args, ArgumentsError as ArgsError, CmdlParser, Formatter as DrmsLogFormatter, Log as DrmsLog, LogLevel as DrmsLogLevel, LogLevelAction as DrmsLogLevelAction, MakeObject, StatusCode as SC
from utils import extract_program_and_module_args

DEFAULT_LOG_FILE = 'mr_log.txt'

REQUEST_ID_PATTERN = r'^JSOC_\d\d\d\d\d\d\d\d_((000)|(00\d)|(0\d\d)|(\d\d\d+))((_X)?_IN)?$'

class StatusCode(SC):
    # fetch success status codes
    REQUEST_COMPLETE = (0, 'request has been completely processed')
    REQUEST_PROCESSING = (1, 'processing request')
    REQUEST_QUEUED = (2, 'request has been queued for processing') # status == 2 exists only in jsoc.export_new; jsoc_fetch op=exp_request creates a record in this table and sets the status value to 2 (or 12)
    # this should go away - if in jsoc.export_new, but not in jsoc.export, then this is the same as REQUEST_QUEUED
    REQUEST_NOT_QUEUED = (6, 'request has not been queued yet') # can no longer happen with jsoc_fetch op=exp_status call
    REQUEST_QUEUED_DEBUG = (12, 'request has been queued for processing')

    # mr status codes
    NOT_PENDING = (101, 'no request pending')
    PENDING = (102, 'request {request_id} is pending')
    REQUEST_CANCELED = (103, 'pending request {request_id} was canceled')

class ErrorCode(ExportErrorCode):
    # fetch error codes
    REQUEST_TOO_LARGE = (3, 'requested payload is too large') # can only happen with jsoc_fetch op=exp_request call, but a user can get this info from the jsoc.export table
    REQUEST_FATAL_ERROR = (4, 'a fatal error occurred during processing')
    REQUEST_NOT_ONLINE = (5, 'exported files are no longer online')

    # other mr error codes
    PARAMETERS = (201, 'failure locating DRMS parameters')
    ARGUMENTS = (202, 'bad arguments')
    LOGGING = (203, 'failure logging messages')
    DB = (204, 'failure executing database command')
    DB_CONNECTION = (205, 'failure connecting to database')
    EXPORT_ACTION = (206, 'failure calling export action')
    CHECK = (207, 'unable to check export request for export user {address}')
    CANCEL = (208, 'unable to cancel export request for export user {address}')
    STATUS = (209, 'unable to check export-request for export user {address}')
    EXPORT_SERVER = (210, 'export-server communication error')
    UNHANDLED_EXCEPTION = (211, 'unhandled exception')

# exceptions
class MrBaseError(ExportError):
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

class ParametersError(MrBaseError):
    _error_code = ErrorCode.PARAMETERS
    # _header = f'if present, then `[cls.header]` will appear at the beginning of the error message'

class ArgumentsError(MrBaseError):
    _error_code = ErrorCode.ARGUMENTS

class LoggingError(MrBaseError):
    _error_code = ErrorCode.LOGGING

class DBError(MrBaseError):
    _error_code = ErrorCode.DB

class DBConnectionError(MrBaseError):
    _error_code = ErrorCode.DB_CONNECTION

class ExportActionError(MrBaseError):
    _error_code = ErrorCode.EXPORT_ACTION

class CheckError(MrBaseError):
    _error_code = ErrorCode.CHECK

class CancelError(MrBaseError):
    _error_code = ErrorCode.CANCEL

class StatusError(MrBaseError):
    _error_code = ErrorCode.STATUS

class ExportServerError(MrBaseError):
    _error_code = ErrorCode.EXPORT_SERVER

class UnhandledExceptionError(MrBaseError):
    _error_code = ErrorCode.UNHANDLED_EXCEPTION


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

# classes
class Arguments(Args):
    _arguments = None

    @classmethod
    def get_arguments(cls, *, is_program, program_name=None, program_args=None, module_args=None, drms_params, refresh=True):
        if cls._arguments is None or refresh:
            try:
                private_db_host = drms_params.get_required('SERVER')
                db_port = drms_params.get_required('DRMSPGPORT')
                db_name = drms_params.get_required('DBNAME')
                db_user = drms_params.get_required('WEB_DBUSER')
                pending_requests_table = drms_params.get_required('EXPORT_PENDING_REQUESTS_TABLE')
                timeout = drms_params.get_required('EXPORT_PENDING_REQUESTS_TIME_OUT')
                download_web_domain = drms_params.get_required('WEB_DOMAIN_PUBLIC')
            except DPMissingParameterError as exc:
                raise ParametersError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

            if is_program:
                try:
                    log_file = path_join(drms_params.get_required('EXPORT_LOG_DIR'), DEFAULT_LOG_FILE)
                except DPMissingParameterError as exc:
                    raise ParametersError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

                args = None

                if program_args is not None and len(program_args) > 0:
                    args = program_args

                parser_args = { 'usage' : '%(prog)s A/address=<registered email address> O/operation=<check/cancel/status> dbhost=<db host> [ -c/--drms-client-type=<ssh/http> ] [ -i/--id=<request ID> ] [ -l/--log-file=<log file path> [ -L/--logging-level=<critical/error/warning/info/debug> ] [ -N/--dbname=<db name> ] [ -P/--dbport=<db port> ] [ -p/--pending-reqs-table=<pending requests db table> ] [ -U/--dbuser=<db user>] [ -t/--timeout=<pending-request timeout> ] [ -w/--webserver=<host> ]' }
                if program_name is not None and len(program_name) > 0:
                    parser_args['prog'] = program_name

                parser = CmdlParser(**parser_args)

                # all arguments are considered optional in argparse (see `prefix_chars`); we can therefore do this:
                # parser.add_argument('-H', 'H', '--dbhost', ...)

                # required
                parser.add_argument('A', 'address', help='the export-registered email address', metavar='<email address>', dest='address', required=True)
                parser.add_argument('O', 'operation', help='the export-request operation to perform (check, cancel, status)', choices=['check', 'cancel', 'status'], metavar='<operation>', dest='operation', default='check', required=True)
                parser.add_argument('dbhost', help='the machine hosting the database that contains export requests from this site', metavar='<db host>', dest='db_host', required=True)

                # optional
                parser.add_argument('-i', '--id', help='request ID; required if operation == status', metavar='<export request ID>', dest='request_id', default=None)
                parser.add_argument('-l', '--log-file', help='the path to the log file', metavar='<log file>', dest='log_file', default=log_file)
                parser.add_argument('-L', '--logging-level', help='the amount of logging to perform; in order of increasing verbosity: critical, error, warning, info, debug', metavar='<logging level>', dest='logging_level', action=DrmsLogLevelAction, default=DrmsLogLevel.ERROR)
                parser.add_argument('-N', '--dbname', help='the name of the database used to manage pending export requests', metavar='<db name>', dest='db_name', default=db_name)
                parser.add_argument('-P', '--dbport', help='The port on the host machine that is accepting connections for the database', metavar='<db host port>', dest='db_port', default=db_port)
                parser.add_argument('-p', '--pending-reqs-table', help='the db table of pending requests', metavar='<pending requests table>', dest='pending_requests_table', default=pending_requests_table)
                parser.add_argument('-t', '--timeout', help='after this number of minutes have elapsed, requests are no longer considered pending', metavar='<timeout>', dest='timeout', default=timeout)
                parser.add_argument('-U', '--dbuser', help='the name of the database user account', metavar='<db user>', dest='db_user', default=db_user)
                parser.add_argument('-w', '--webserver', help='the webserver invoking this script', metavar='<webserver>', action=create_webserver_action(drms_params), dest='webserver', default=name_to_ws_obj(None, drms_params))

                arguments = Arguments(parser=parser, args=args)
            else:
                # module invocation
                def extract_module_args(*, address, operation, db_host, request_id=None, log=None, db_name=db_name, db_port=db_port, pending_requests_table=pending_requests_table, timeout=timeout, db_user=db_user, webserver=None):
                    arguments = {}

                    arguments['address'] = address
                    arguments['operation'] = operation
                    arguments['db_host'] = db_host
                    arguments['request_id'] = request_id
                    arguments['db_name'] = db_name
                    arguments['db_port'] = db_port
                    arguments['pending_requests_table'] = pending_requests_table
                    arguments['timeout'] = timeout
                    arguments['db_user'] = db_user
                    arguments['webserver'] = name_to_ws_obj(webserver, drms_params) # sets webserver.public = True if webserver is None

                    PendingRequestAction.set_log(log)

                    return arguments

                module_args_dict = extract_module_args(**module_args)
                arguments = Arguments(parser=None, args=module_args_dict)

            # if the caller has specified a public webserver, make sure that the db specified is not the private one
            if arguments.webserver.public and arguments.db_host == private_db_host:
                raise ArgumentsError(error_message=f'cannot specify private db server to handle public webserver requests')

            # if operation == status, then id is required
            if arguments.operation == 'status' and arguments.request_id is None:
                raise ArgumentsError(error_message=f'must specify `id` argument if operation is `status`')

            arguments.private_db_host = private_db_host
            arguments.download_web_domain = download_web_domain

            cls._arguments = arguments

        return cls._arguments

class OperationFactory():
    def __new__(cls, *, operation_name, address, log=None):
        operation = None
        if operation_name.lower() == CheckOperation._name:
            operation = CheckOperation(address, log)
        elif operation_name.lower() == CancelOperation._name:
            operation = CancelOperation(address, log)
        elif operation_name.lower() == StatusOperation._name:
            operation = StatusOperation(address, log)
        else:
            raise ArgumentsError(error_message=f'invalid operation type {operation}')

        log.write_debug([ f'[ OperationFactory ] created `{operation._name}` operation' ])
        return operation

# operations
class Operation():
    def __init__(self, address, log):
        self._address = address
        self._log = log
        self._request_id = None
        self._start_time = None
        self._response = None

    def __str__(self):
        return self._name

    def __call__(self, cursor, pending_requests_table, timeout):
        cmd = f"SELECT request_id, start_time FROM {pending_requests_table} WHERE lower(address) = lower('{self._address}') AND CURRENT_TIMESTAMP - start_time < interval '{str(timeout)} minutes'"
        self._log.write_debug([ f'[ Operation.__call__ ] executing SQL `{cmd}`'])

        try:
            cursor.execute(cmd)
            rows = cursor.fetchall()
            if len(rows) > 1:
                # only one request per address allowed currently
                raise self._generate_exception(exc=DbError, error_message=f'unexpected number of rows returned ({cmd})')
        except psycopg2.Error as exc:
            # handle database-command errors
            raise self._generate_exception(exc=DbError, exc_info=sys_exc_info(), error_message=str(exc))

        if len(rows) != 0:
            self._request_id = rows[0][0]
            self._start_time = rows[0][1]
            self._log.write_debug([ f'[ Operation.__call__ ] PG returned one row {(str(self._request_id), self._start_time.strftime("%Y-%m-%d %T"))}'])
        else:
            self._log.write_debug([ f'[ Operation.__call__ ] no export requests for user {self._address}' ])

    def _generate_exception(self, exc, exc_info=None, error_message=None):
        description = self._exception._error_code.description(address=self._address)
        message = f'{description}: {error_message}' if error_message is not None and len(error_message) > 0 else description
        return exc(error_message=message)

    @property
    def response(self):
        return self._response

class CheckOperation(Operation):
    _name = 'check'
    _exception = CheckError

    def __call__(self, *, cursor, pending_requests_table, timeout=timedelta(minutes=60)):
        super().__call__(cursor, pending_requests_table, timeout)

        if not self._request_id:
            self._response = NotPendingResponse.generate_response(address=self._address)
        else:
            self._response = PendingResponse.generate_response(address=self._address, request_id=self._request_id, start_time=self._start_time.strftime('%Y-%m-%d %T'))

class CancelOperation(Operation):
    _name = 'cancel'
    _exception = CancelError

    def __call__(self, *, cursor, pending_requests_table, timeout=timedelta(minutes=60)):
        # first run the Operation.process() code to obtain the request_id
        super().__call__(cursor, pending_requests_table, timeout)

        # then run the code to delete the pending request
        if not self._request_id:
            self._response = NotPendingResponse.generate_response(address=self._address)
        else:
            cmd = f"DELETE FROM {pending_requests_table} WHERE lower(address) = lower('{self._address}') AND lower(request_id) = lower('{self._request_id}')"
            self._log.write_debug([ f'[ CancelOperation.__call__ ] executing SQL `{cmd}`'])

            try:
                cursor.execute(cmd)
            except psycopg2.Error as exc:
                # cannot delete the pending request from the pending-requests db table
                error_msg = f'cannot delete pending request with id={self._request_id} and start_time={self._start_time.strftime("%Y-%m-%d %T")} ({str(exc)})'
                raise self._generate_exception(exc=self._exception, exc_info=sys_exc_info(), error_message=error_msg)

            self._response = CancelResponse.generate_response(address=self._address, request_id=self._request_id, start_time=self._start_time.strftime('%Y-%m-%d %T'))

class StatusOperation(Operation):
    _name = 'status'
    _exception = StatusError

    def __call__(self, *, cursor, pending_requests_table, timeout=timedelta(minutes=60), connection, request_id, db_host, download_web_domain):
        delete_request = None
        message = { 'request_type' : 'export_status', 'address' : self._address, 'request_id' : request_id, 'db_host' : db_host }
        response = send_request(message, connection, self._log)
        export_status_dict = json_loads(response)

        if export_status_dict.get('export_server_status') == 'export_server_error':
            raise ExportServerError(error_message=f'{export_status_dict["error_message"]}')

        status_code = self.get_request_status_code(export_status_dict)

        self._log.write_debug([ f'[ StatusOperation.__call__ ] export request {request_id} status is `{status_code.description(address=self._address)}`' ])

        if not isinstance(status_code, ErrorCode):
            if status_code == StatusCode.REQUEST_COMPLETE:
                # request should NOT be in pending requests table
                delete_request = True
            else:
                # request should be in pending requests table
                delete_request = False
        else:
            # GRRRRR! if an error occurred, then export_status_dict contains almost no info, not even the request ID; add
            # property that does contain the request id value
            export_status_dict['REQUEST_ID'] = request_id

            # request should NOT be in pending requests table
            delete_request = True

        # to make URLs for the downloadables
        export_status_dict['DOWNLOAD_WEB_DOMAIN'] = download_web_domain

        super().__call__(cursor, pending_requests_table, timeout)
        if self._request_id is not None and self._request_id == request_id:
            # request IS in pending requests table
            if delete_request:
                # request should NOT have been in the pending-requests table; log as a warning
                self._log.write_warning([ f'[ StatusOperation.__call__ ] export request `{self._request_id}` has completed, but it was not removed from the pending-requests table; removing orphaned request now'])

                # and then call CancelOperation code to remove it from the pending-requests table
                try:
                    operation = OperationFactory(operation_name=CancelOperation._name, address=self._address, log=self._log)
                    operation(cursor=cursor, pending_requests_table=pending_requests_table, timeout=timeout)
                except MrBaseError as exc:
                    error_msg = f'unable to cancel orphaned pending request `{self._request_id}` for user `{self._address}`'
                    raise self._generate_exception(exc=self._exception, exc_info=sys_exc_info(), error_message=error_msg)
        else:
            # cannot locate the pending request in the pending-requests db table
            if not delete_request:
                # request should have been in the pending-requests table; log as a warning
                self._log.write_warning([ f'[ StatusOperation.__call__ ] export system is processing request {self._request_id}, but that request is not in the pending-requests table'])

        self._response = self.get_response(export_status_dict)

    def get_response_dict(self, export_status_dict):
        self._log.write_debug([ f'[ StatusOperation.get_response_dict ]' ])
        error_code, status_code = self.get_request_status(export_status_dict) # 2/12 ==> asynchronous processing, 0 ==> synchronous processing
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
            self._log.write_debug([ f'[ StatusOperation.get_response_dict ] error calling fetch `{error_code.description()}` for request `{str(request_id)}`'])
            # a fetch error occurred (`fetch_error` has the error message returned by fetch, `status_description` has the IR error message)
            fetch_error = export_status_dict.get('error') # None, unless an error occurred
            status_code = error_code
        else:
            # no fetch error occurred
            self._log.write_debug([ f'[ StatusOperation.get_response_dict ] NO error calling fetch, status `{status_code.description(address=self._address)}` for request `{str(request_id)}`'])

            if status_code == StatusCode.REQUEST_COMPLETE:
                self._log.write_debug([ f'[ StatusOperation.get_response_dict ] exported data were generated for request `{request_id}`'])

                # not all of these attributes are provides in all responses, so use `get()`
                export_directory = export_status_dict.get('dir')
                keywords_file = export_status_dict.get('keywords')
                tar_file = export_status_dict.get('tarfile')
                method = export_status_dict.get('method')

                download_web_domain = export_status_dict.get('DOWNLOAD_WEB_DOMAIN')
                if download_web_domain is not None and export_directory is not None:
                    # URL of the export SU
                    request_url = urlunsplit(('http', download_web_domain, export_directory, None, None))

                if method is not None:
                    try:
                        access = method[:3]
                    except ValueError:
                        access = 'url'

                scheme = 'https' if access == 'url' else 'ftp'

                package = { 'type' : None if tar_file is None else 'tar', 'file_name' : None if tar_file is None else tar_file }

                file_information = export_status_dict.get('data')
                if file_information is not None:
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
            response_dict.update({ 'error_code' : error_code, 'error_message' : fetch_error, 'request_id' : request_id, 'contact' : contact })
        else:
            response_dict.update({ 'status_code' : status_code, 'fetch_error' : fetch_error, 'request_id' : request_id, 'file_format' : file_format, 'package' : package, 'access' : access, 'number_records' : record_count, 'number_files' : file_count, 'mb_exported' : mb_exported, 'sums_directory' : export_directory, 'keywords_file' : keywords_file, 'tar_file' : tar_file, 'request_url' : request_url, 'export_data' : data })

        self._log.write_debug([ f'[ StatusOperation.get_response_dict] response dictionary `{str(response_dict)}`'])

        return (error_code is not None, response_dict)

    def get_response(self, export_status_dict):
        self._log.write_debug([ f'[ StatusOperation.get_response ]' ])
        # was there an error?
        error_code = None

        error_occurred, response_dict = self.get_response_dict(export_status_dict)

        if error_occurred:
            # an error occurred in jsoc_fetch; this is not the same thing as an error happening in MR so we have to manually
            # create an error response (the status_code will not be a MR error code, but the fetch_status will be a fetch error code)
            response = ErrorResponse.generate_response(**response_dict)
        else:
            response = StatusResponse.generate_response(**response_dict)

        return response

    def get_request_status(self, export_status_dict):
        self._log.write_debug([ f'[ StatusOperation.get_status ]' ])
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

    def get_request_status_code(self, export_status_dict):
        error_code, status_code = self.get_request_status(export_status_dict)
        code = error_code if error_code is not None else status_code

        return code

# responses
class ManageRequestResponse(Response):
    _status_code = None
    _comment = None

class NotPendingResponse(ManageRequestResponse):
    _status_code = StatusCode.NOT_PENDING
    _comment = 'no existing export request for export user {address}'

    @classmethod
    def generate_response(cls, *, status_code=None, address, **response_dict):
        response_dict['comment'] = cls._comment.format(address=address)
        return super().generate_response(status_code=status_code, address=address, **response_dict)

class PendingResponse(ManageRequestResponse):
    _status_code = StatusCode.PENDING
    _comment = 'existing export request for export user {address} [ request_id={request_id}, start_time={start_time} ]'

    @classmethod
    def generate_response(cls, *, status_code=None, address, request_id, start_time, **response_dict):
        response_dict['comment'] = cls._comment.format(address=address, request_id=request_id, start_time=start_time)
        # request_id needed when export.py converts ManageRequestResponse into dict
        return super().generate_response(status_code=status_code, address=address, request_id=request_id, start_time=start_time, **response_dict)

class CancelResponse(ManageRequestResponse):
    _status_code = StatusCode.REQUEST_CANCELED
    _comment = 'existing export request for export user {address} [ request_id={request_id}, start_time={start_time} ] was canceled'

    @classmethod
    def generate_response(cls, *, status_code=None, address, request_id, start_time, **response_dict):
        response_dict['comment'] = cls._comment.format(address=address, request_id=request_id, start_time=start_time)
        # request_id needed when export.py converts ManageRequestResponse into dict
        return super().generate_response(status_code=status_code, address=address, request_id=request_id, start_time=start_time, **response_dict)

class StatusResponse(ManageRequestResponse):
    _status_code = None
    _comment = f'`status_description` describes fetch status'

    @classmethod
    def generate_response(cls, *, status_code=None, **response_dict):
        response_dict['comment'] = cls._comment
        return super().generate_response(status_code=status_code, **response_dict)


# for use in export web app
from action import Action
class PendingRequestAction(Action):
    actions = [ 'check_pending_request', 'cancel_pending_request', 'get_export_status' ]

    _reg_ex = None
    _log = None

    def __init__(self, *, method, address, db_host, request_id=None, log=None, db_name=None, db_port=None, pending_requests_table=None, timeout=None, db_user=None, webserver=None):
        self._method = getattr(self, method)
        self._address = address
        self._db_host = db_host
        self._options = {}
        self._options['request_id'] = request_id
        self._options['log'] = log
        self._options['db_name'] = db_name
        self._options['db_port'] = db_port
        self._options['pending_requests_table'] = pending_requests_table
        self._options['timeout'] = timeout
        self._options['db_user'] = db_user
        self._options['webserver'] = webserver # host name - gets converted to object in `get_arguments()`

    def check_pending_request(self):
        response = perform_action(action_obj=self, is_program=False, operation='check', address=self._address, db_host=self._db_host, options=self._options)
        return response

    def cancel_pending_request(self):
        response = perform_action(action_obj=self, is_program=False, operation='cancel', address=self._address, db_host=self._db_host, options=self._options)
        return response

    def get_export_status(self):
        response = perform_action(action_obj=self, is_program=False, operation='status', address=self._address, db_host=self._db_host, options=self._options)
        return response

    @classmethod
    def get_reg_ex(cls, logging_level=None):
        if cls._reg_ex is None:
            cls._reg_ex = re_compile(REQUEST_ID_PATTERN)

        return cls._reg_ex

    @classmethod
    def is_valid_request_id(cls, address, log):
        cls.set_log(log)
        reg_ex = cls.get_reg_ex()
        return reg_ex.match(address) is not None

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

def requires_private_db(request_id):
    reg_ex = PendingRequestAction.get_reg_ex()
    match = reg_ex.match(request_id)

    if match is None:
        raise ArgumentsError(error_message=f'invalid request ID {request_id}')

    # private if ends in '_' ['X' '_'] 'I' 'N'
    return True if match.group(6) is not None else False

def send_request(request, connection, log):
    json_message = json_dumps(request)
    send_message(connection, json_message)
    message = get_message(connection)

    return message

def perform_action(*, action_obj, is_program, program_name=None, **kwargs):
    response = None
    log = None

    try:
        program_args, module_args = extract_program_and_module_args(is_program=is_program, **kwargs)

        try:
            drms_params = DRMSParams()

            if drms_params is None:
                raise ParametersError(error_message='unable to locate DRMS parameters package')

            arguments = Arguments.get_arguments(is_program=is_program, program_name=program_name, program_args=program_args, module_args=module_args, drms_params=drms_params)
        except ArgsError as exc:
            raise ArgumentsError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')
        except Exception as exc:
            raise ArgumentsError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

        if is_program:
            try:
                formatter = DrmsLogFormatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
                log = DrmsLog(arguments.log_file, arguments.logging_level, formatter)
                PendingRequestAction.set_log(log)
            except Exception as exc:
                raise LoggingError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')
        else:
            log = action_obj.log

        if is_program:
            log.write_debug([ f'[ perform_action ] program invocation' ])
        else:
            log.write_debug([ f'[ perform_action ] module invocation' ])

        log.write_debug([ f'[ perform_action ] action arguments: {str(arguments)}' ])

        try:
            operation = OperationFactory(operation_name=arguments.operation, address=arguments.address, log=log)
            response = None

            try:
                # `arguments.db_host` contains the pending requests table - the public webserver uses a public db host,
                # and the private webserver uses the private db host
                with psycopg2.connect(host=arguments.db_host, port=str(arguments.db_port), database=arguments.db_name, user=arguments.db_user) as conn:
                    with conn.cursor() as cursor:
                        if isinstance(operation, StatusOperation):
                            # there are two types of request IDs that the public webserver supports:
                            #   JSOC_20210821_050 - external DB handled the request
                            #   JSOC_20210821_051_X_IN - internal DB handled the request
                            try:
                                private_db_needed = requires_private_db(arguments.request_id)
                            except Exception as exc:
                                raise ExportActionError(exc_info=sys_exc_info(), error_message=str(exc))

                            if not private_db_needed and arguments.db_host == arguments.private_db_host:
                                raise ArgumentsError(error_message=f'request ID {arguments.request_id} requires public database access, but client is the private webserver (uses public database host `{arguments.db_host}`)')

                            if private_db_needed:
                                resolved_db_host = arguments.private_db_host
                                log.write_debug([ f'[ perform_action ] accessing private database host {arguments.db_host}' ])
                            else:
                                # must be public db host - a private host would have caused ArgumentsError exception
                                resolved_db_host = arguments.db_host
                                log.write_debug([ f'[ perform_action ] access public database host {arguments.db_host}' ])

                            try:
                                # use socket server to call jsoc_fetch
                                nested_arguments = ss_get_arguments(is_program=False, module_args={})

                                with Connection(server=nested_arguments.server, listen_port=nested_arguments.listen_port, timeout=nested_arguments.message_timeout, log=log) as connection:
                                    operation(cursor=cursor, pending_requests_table=arguments.pending_requests_table, timeout=arguments.timeout, connection=connection, request_id=arguments.request_id, db_host=resolved_db_host, download_web_domain=arguments.download_web_domain)

                                    message = { 'request_type' : 'quit' }
                                    send_request(message, connection, log)

                            except ExpServerBaseError as exc:
                                raise ExportServerError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')
                            except Exception as exc:
                                raise ExportServerError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')
                        else:
                            operation(cursor=cursor, pending_requests_table=arguments.pending_requests_table, timeout=arguments.timeout)

                        response = operation.response
            except psycopg2.OperationalError as exc:
                # closes the cursor and connection
                raise DBConnectionError(exc_info=sys_exc_info(), error_message=f'unable to connect to the database: {str(exc)}')
        except Exception as exc:
            raise ExportActionError(exc_info=sys_exc_info(), error_message=str(exc))
    except ExpServerBaseError as exc:
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

    return response

if __name__ == "__main__":
    response = perform_action(action_obj=None, is_program=True)
    print(response.generate_json())

    # Always return 0. If there was an error, an error code (the 'status' property) and message (the 'statusMsg' property) goes in the returned HTML.
    sys_exit(0)
else:
    pass

# add py: base/export/scripts/checkexpdbserver, base/export/scripts/jsocextfetch, base/export/scripts/jsocextinfo, base/export/scripts/showextinfo, base/export/scripts/showextseries, base/export/scripts/seriesinfo, base/export/scripts/checkAddress
# add C programs (in drms/base): drms-export-as-fits, jsoc_info, show_info, drms_parserecset; AAAH! use securedrms.py

# dir structure - use links in code tree to create this
# proj/export/webapps/drmsexport
#   drmsexport
#     __init__.py
#     action.py
#     arguments.py
#     check_address.py
#     check_dbserver.py
#     drmsparams.py
#     error.py
#     export.py
#     export_wsgi.py
#     seriesinfo.py (like show_info/jsoc_info/showextinfo/jsocextinfo)
#     initiate_request.py (like jsoc_fetch/jsocextfetch)
#     manage_request.py
#     securedrms.py(copied from base/libs/py)
#     statuscode.py
#  setup.py
