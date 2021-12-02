#!/usr/bin/env python3

from argparse import Action as ArgsAction
import inspect
from json import loads as json_loads, dumps as json_dumps
from os import environ
from os.path import join as path_join
from sys import exc_info as sys_exc_info, exit as sys_exit
from drms_export import Connection, Error as ExportError, ErrorCode as ExportErrorCode, ExpServerBaseError, get_arguments as ss_get_arguments, get_message, Response, send_message
from drms_parameters import DRMSParams, DPMissingParameterError
from drms_utils import Arguments as Args, Choices, CmdlParser, Formatter as DrmsLogFormatter, Log as DrmsLog, LogLevel as DrmsLogLevel, LogLevelAction as DrmsLogLevelAction, MakeObject, StatusCode as SC
from utils import extract_program_and_module_args

DEFAULT_LOG_FILE = 'ps_log.txt'

class StatusCode(SC):
    SUCCESS = (0, 'success')

class ErrorCode(ExportErrorCode):
    PARAMETERS = (1, 'failure locating DRMS parameters')
    ARGUMENTS = (2, 'bad arguments')
    LOGGING = (3, 'failure logging messages')
    EXPORT_SERVER = (4, 'export-server communication error')
    UNHANDLED_EXCEPTION = (5, 'unhandled exception')

class PsBaseError(ExportError):
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

class ParametersError(PsBaseError):
    _error_code = ErrorCode.PARAMETERS

class ArgumentsError(PsBaseError):
    _error_code = ErrorCode.ARGUMENTS

class LoggingError(PsBaseError):
    _error_code = ErrorCode.LOGGING

class ExportServerError(PsBaseError):
    _error_code = ErrorCode.EXPORT_SERVER

class UnhandledExceptionError(PsBaseError):
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

class Arguments(Args):
    _arguments = None

    @classmethod
    def get_arguments(cls, *, is_program, program_name=None, program_args=None, module_args=None, drms_params, refresh=True):
        if cls._arguments is None or refresh:
            try:
                log_file = path_join(drms_params.get_required('EXPORT_LOG_DIR'), DEFAULT_LOG_FILE)
                private_db_host = drms_params.get_required('SERVER')
                db_port = int(drms_params.get_required('DRMSPGPORT'))
                db_name = drms_params.get_required('DBNAME')
                db_user = drms_params.get_required('WEB_DBUSER')
            except DPMissingParameterError as exc:
                raise ParametersError(error_message=str(exc))

            if is_program:
                args = None

                if program_args is not None and len(program_args) > 0:
                    args = program_args

                parser_args = { 'usage' : f'%(prog)s specification=<specification> dbhost=<db host> [ -l/--log-file=<log file path> [ -L/--logging-level=<critical/error/warning/info/debug> ] [ -H/--dbhost=<db host> ] [ -N/--dbname=<db name> ] [ -P/--dbport=<db port> ] [ -U/--dbuser=<db user>] [ --webserver=<host> ]' }

                if program_name is not None and len(program_name) > 0:
                    parser_args['prog'] = program_name

                parser = CmdlParser(**parser_args)

                # required
                parser.add_argument('specification', help='the DRMS record-set specification', metavar='<specification>', dest='specification', required=True)
                parser.add_argument('dbhost', help='the machine hosting the database that contains export requests from this site', metavar='<db host>', dest='db_host', required=True)

                # optional
                parser.add_argument('-l', '--log-file', help='the path to the log file', metavar='<log file>', dest='log_file', default=log_file)
                parser.add_argument('-L', '--logging-level', help='the amount of logging to perform; in order of increasing verbosity: critical, error, warning, info, debug', metavar='<logging level>', dest='logging_level', action=DrmsLogLevelAction, default=DrmsLogLevel.ERROR)
                parser.add_argument('-N', '--dbname', help='the name of the database that contains export requests', metavar='<db name>', dest='db_name', default=db_name)
                parser.add_argument('-P', '--dbport', help='the port on the host machine that is accepting connections for the database', metavar='<db host port>', dest='db_port', type=int, default=db_port)
                parser.add_argument('-U', '--dbuser', help='the name of the database user account', metavar='<db user>', dest='db_user', default=db_user)
                parser.add_argument('-w', '--webserver', help='the webserver invoking this script', metavar='<webserver>', action=create_webserver_action(drms_params), dest='webserver', default=name_to_ws_obj(None, drms_params))

                arguments = Arguments(parser=parser, args=args)
            else:
                def extract_module_args(*, specification, db_host, log_file=log_file, logging_level='error', db_name=db_name, db_port=db_port, db_user=db_user, webserver=None):
                    arguments = {}

                    arguments['specification'] = specification
                    arguments['db_host'] = db_host
                    arguments['log_file'] = log_file
                    arguments['logging_level'] = DrmsLogLevelAction.string_to_level(logging_level)
                    arguments['db_name'] = db_name
                    arguments['db_port'] = db_port
                    arguments['db_user'] = db_user
                    arguments['webserver'] = name_to_ws_obj(webserver, drms_params) # sets webserver.public = True if webserver is None

                    return arguments

                # dict
                module_args_dict = extract_module_args(**module_args)
                arguments = Arguments(parser=None, args=module_args_dict)

            # if the caller has specified a public webserver, make sure that the db specified is not the private one
            if arguments.webserver.public and arguments.db_host == private_db_host:
                raise ArgumentsError(error_message=f'cannot specify private db server to handle public webserver requests')

            arguments.private_db_host = private_db_host

            cls._arguments = arguments

        return cls._arguments

# for use in export web app
from action import Action
class ParseSpecificationAction(Action):
    actions = [ 'parse_specification' ]

    _log = None

    def __init__(self, *, method, specification, db_host, logging_level=None, db_name=None, db_port=None, db_user=None, webserver=None):
        self._method = getattr(self, method)
        self._specification = specification
        self._db_host = db_host # host webserver uses (private webserver uses private db host)
        self._options = {}
        self._options['logging_level'] = logging_level
        self._options['db_name'] = db_name
        self._options['db_port'] = db_port
        self._options['db_user'] = db_user
        self._options['webserver'] = webserver

    def parse_specification(self):
        # returns dict
        prog_name = inspect.currentframe().f_code.co_name
        response = perform_action(action_obj=self, is_program=False, program_name=prog_name, specification=self._specification, db_host=self._db_host, options=self._options)
        return response

    @property
    def log(self):
        return self.__class__._log

    @log.setter
    def log(self, log):
        self.__class__._log = log

def send_request(request, connection, log):
    json_message = json_dumps(request)

    log.write_debug([ f'[ send_request ] sending message to server:' ])
    log.write_debug([ f'[ send_request ] {json_message}' ])
    send_message(connection, json_message)

    message = get_message(connection)
    log.write_debug([ f'[ send_request ] server response:' ])
    log.write_debug([ f'[ send_request ] {message}' ])

    return message

def perform_action(*, action_obj, is_program, program_name=None, **kwargs):
    # catch all expections so we can always generate a response
    response = None
    log = None

    try:
        # removes each options with a value of None
        program_args, module_args = extract_program_and_module_args(is_program=is_program, **kwargs)

        try:
            drms_params = DRMSParams()

            if drms_params is None:
                raise ParametersError(error_message=f'unable to locate DRMS parameters package')

            arguments = Arguments.get_arguments(is_program=is_program, program_name=program_name, program_args=program_args, module_args=module_args, drms_params=drms_params)
        except Exception as exc:
            raise ArgumentsError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

        if action_obj is None or action_obj.log is None:
            try:
                formatter = DrmsLogFormatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
                log = DrmsLog(arguments.log_file, arguments.logging_level, formatter)
                if action_obj is not None:
                    action_obj.log = log
            except Exception as exc:
                raise LoggingError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')
        else:
            log = action_obj.log

        log.write_debug([ f'[ perform_action ] action arguments: {str(arguments)}' ])

        try:
            # use socket server to call drms_parserecset
            nested_arguments = ss_get_arguments(is_program=False, module_args={})

            with Connection(server=nested_arguments.server, listen_port=nested_arguments.listen_port, timeout=nested_arguments.message_timeout, log=log) as connection:
                message = { 'request_type' : 'parse_specification', 'specification' : arguments.specification }
                response = send_request(message, connection, log)

                # message is raw JSON from drms_parserecset
                response_dict = json_loads(response)

                if 'status' in response_dict and response_dict['status'] == 'export_server_error':
                    raise ExportServerError(error_message=f'{response_dict["error_message"]}')
                message = { 'request_type' : 'quit' }
                send_request(message, connection, log)

            if response_dict['errMsg'] is not None:
                raise ExportServerError(error_message=f'failure parsing record-set specification {arguments.export_arguments["specification"]}: {response_dict["errMsg"]}')
        except ExpServerBaseError as exc:
            raise ExportServerError(exc_info=sys_exc_info(), error_message=f'{exc.message}')
        except PsBaseError:
            raise
        except Exception as exc:
            raise ExportServerError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

        response = Response.generate_response(status_code=StatusCode.SUCCESS, **response_dict)
    except PsBaseError as exc:
        response = exc.response
        error_message = exc.message

        if log:
            log.write_error([ error_message ])
        else:
            print(error_message)
    except Exception as exc:
        response = UnhandledExceptionError(exc_info=sys_exc_info(), error_message=f'{str(exc)}').response
        error_message = str(exc)

        if log:
            log.write_error([ error_message ])
        else:
            print(error_message)

    return response

if __name__ == '__main__':
    response = perform_action(action_obj=None, is_program=True)
    print(response.generate_json())

    sys_exit(0)
else:
    pass
