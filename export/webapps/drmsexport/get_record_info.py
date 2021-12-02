#!/usr/bin/env python3

from argparse import Action as ArgsAction
from copy import deepcopy
from functools import lru_cache
from json import loads as json_loads, dumps as json_dumps
from os.path import join as path_join
from sys import exc_info as sys_exc_info, exit as sys_exit

from drms_export import Connection, ExpServerBaseError, Error as ExportError, ErrorCode as ExportErrorCode, ErrorResponse, Response, get_arguments as ss_get_arguments, get_message, send_message
from drms_parameters import DRMSParams, DPMissingParameterError
from drms_utils import Arguments as Args, ArgumentsError as ArgsError, Choices, CmdlParser, Formatter as DrmsLogFormatter, ListAction, Log as DrmsLog, LogLevel as DrmsLogLevel, LogLevelAction as DrmsLogLevelAction, MakeObject, StatusCode as SC
from utils import extract_program_and_module_args, get_db_host

DEFAULT_LOG_FILE = 'gri_log.txt'

class StatusCode(SC):
    # info success status codes
    SUCCESS = (0, 'success')

class ErrorCode(ExportErrorCode):
    # info failure error codes
    FAILURE = (1, 'info failure')

    PARAMETERS = (101, 'failure locating DRMS parameters')
    ARGUMENTS = (102, 'bad arguments')
    LOGGING = (103, 'failure logging messages')
    DRMS_CLIENT = (104, 'drms client error')
    EXPORT_ACTION = (105, 'failure calling export action')
    RESPONSE = (106, 'unable to generate valid response')

class GriBaseError(ExportError):
    def __init__(self, *, exc_info=None, error_message=None):
        if exc_info is not None:
            self.exc_info = exc_info
            e_type, e_obj, e_tb = exc_info

            if error_message is None:
                error_message = f'{e_type.__name__}: {str(e_obj)}'
            else:
                error_message = f'{error_message} [ {e_type.__name__}: {str(e_obj)} ]'

        super().__init__(error_message=error_message)

class ParametersError(GriBaseError):
    _error_code = ErrorCode.PARAMETERS

class ArgumentsError(GriBaseError):
    _error_code = ErrorCode.ARGUMENTS

class LoggingError(GriBaseError):
    _error_code = ErrorCode.LOGGING

class DRMSClientError(GriBaseError):
    _error_code = ErrorCode.DRMS_CLIENT

class ExportActionError(GriBaseError):
    _error_code = ErrorCode.EXPORT_ACTION

class ResponseError(GriBaseError):
    _error_code = ErrorCode(ErrorCode.RESPONSE)

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
                raise ParametersError(exc_info=sys_exc_info(), error_message=str(exc))

            if is_program:
                args = None

                if program_args is not None and len(program_args) > 0:
                    args = program_args

                parser_args = { 'usage' : '%(prog)s specification=<DRMS record-set specification> dbhost=<db host> [ -c/--drms-client-type=<ssh/http>] [ -k/--keywords=<keywords> ] [ -K/--links=<links> ] [ -l/--log-file=<log file path> [ -L/--logging-level=<critical/error/warning/info/debug> ] [ -N/--dbname=<db name> ] [ -P/--dbport=<db port> ] [ -s/--segments=<segments> ] [ -U/--dbuser=<db user>] [ -w/--webserver=<host> ] ' }
                if program_name is not None and len(program_name) > 0:
                    parser_args['prog'] = program_name

                parser = CmdlParser(**parser_args)

                # required
                parser.add_argument('specification', help='the DRMS record-set specification identifying the records for which information is to be obtained', metavar='<DRMS record-set specification>', dest='specification', required=True)
                parser.add_argument('dbhost', help='the machine hosting the database that contains export requests from this site', metavar='<db host>', dest='db_host', required=True)

                # optional
                parser.add_argument('-c', '--drms-client-type', help='securedrms client type (ssh, http)', choices=[ 'ssh', 'http' ], dest='drms_client_type', default='ssh')
                parser.add_argument('-k', '--keywords', help='list of keywords for which information is returned', action=ListAction, dest='keywords', default=None)
                parser.add_argument('-K', '--links', help='list of links for which information is returned', action=ListAction, dest='links', default=None)
                parser.add_argument('-l', '--log-file', help='the path to the log file', metavar='<log file>', dest='log_file', default=log_file)
                parser.add_argument('-L', '--logging-level', help='the amount of logging to perform; in order of increasing verbosity: critical, error, warning, info, debug', metavar='<logging level>', dest='logging_level', action=DrmsLogLevelAction, default=DrmsLogLevel.ERROR)
                parser.add_argument('-N', '--dbname', help='the name of the database that contains export requests', metavar='<db name>', dest='db_name', default=db_name)
                parser.add_argument('-n', '--number-records', help='the maximum number of records for which information is returned', metavar='<maximum number of records>', dest='number_records', type=int, default=None)
                parser.add_argument('-P', '--dbport', help='the port on the host machine that is accepting connections for the database', metavar='<db host port>', dest='db_port', type=int, default=db_port)
                parser.add_argument('-s', '--segments', help='list of segments for which information is returned', action=ListAction, dest='segments', default=None)
                parser.add_argument('-U', '--dbuser', help='the name of the database user account', metavar='<db user>', dest='db_user', default=db_user)
                parser.add_argument('-w', '--webserver', help='the webserver invoking this script', metavar='<webserver>', action=create_webserver_action(drms_params), dest='webserver', default=name_to_ws_obj(None, drms_params))

                arguments = Arguments(parser=parser, args=args)
                arguments.drms_client = None
            else:
                # `program_args` has all `arguments` values, in final form; validate them
                def extract_module_args(*, specification, db_host, drms_client_type='ssh', drms_client=None, keywords=None, links=None, log_file=log_file, logging_level='error', db_name=db_name, number_records=None, db_port=db_port, segments=None, db_user=db_user, webserver=None):
                    arguments = {}

                    arguments['specification'] = specification
                    arguments['db_host'] = db_host
                    arguments['drms_client_type'] = drms_client_type
                    arguments['drms_client'] = drms_client
                    arguments['keywords'] = keywords # list
                    arguments['links'] = links # list
                    arguments['log_file'] = log_file
                    arguments['logging_level'] = DrmsLogLevelAction.string_to_level(logging_level)
                    arguments['db_name'] = db_name
                    arguments['number_records'] = number_records
                    arguments['db_port'] = db_port
                    arguments['segments'] = segments # list
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

            if not GetRecordInfoAction.is_valid_specification(arguments.specification, arguments.db_host, arguments.webserver.host, arguments.logging_level._fullname):
                raise ArgumentsError(error_message=f'invalid record-set specification')

            parsed_specification = get_parsed_specification(arguments.specification, arguments.db_host, arguments.webserver.host)

            if not parsed_specification.attributes.hasfilts and arguments.number_records is None:
                raise ArgumentsError(error_message=f'must specify either a record-set filter, or a maximum number of records')

            cls._arguments = arguments

        return cls._arguments

def get_response(client_response_dict, log):
    log.write_debug([ f'[ get_response ]' ])

    response_dict = deepcopy(client_response_dict)
    info_status = response_dict['status']

    if info_status == 0:
        response_dict['status_code'] = StatusCode.SUCCESS
    elif info_status == 1:
        response_dict['status_code'] = StatusCode.FAILURE
    else:
        # there should be no other status possible
        response_dict['status_code'] = StatusCode.FAILURE

    return Response.generate_response(**response_dict)

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
                raise ParametersError(error_message='unable to locate DRMS parameters package')

            arguments = Arguments.get_arguments(is_program=is_program, program_name=program_name, program_args=program_args, module_args=module_args, drms_params=drms_params)
        except ArgsError as exc:
            raise ArgumentsError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')
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

        if is_program:
            log.write_debug([ f'[ perform_action ] program invocation' ])
        else:
            log.write_debug([ f'[ perform_action ] module invocation' ])

        log.write_debug([ f'[ perform_action ] action arguments: {str(arguments)}' ])

        parsed_specification = GetRecordInfoAction.get_parsed_specification(arguments.specification, arguments.db_host, arguments.webserver.host)
        series = []
        for subset in parsed_specification.attributes.subsets:
            series.append(subset.seriesname)

        resolved_db_host = get_db_host(webserver=arguments.webserver, series=series, private_db_host=arguments.private_db_host, db_host=arguments.db_host, db_port=arguments.db_port, db_name=arguments.db_name, db_user=arguments.db_user, exc=ExportActionError, log=log)

        if resolved_db_host is None:
            raise ArgumentsError(error_message=f'unable to determine db host suitable for series {", ".join(series)}')

        # make request to socket server
        try:
            # use socket server to call drms_parserecset
            nested_arguments = ss_get_arguments(is_program=False, module_args={})

            with Connection(server=nested_arguments.server, listen_port=nested_arguments.listen_port, timeout=nested_arguments.message_timeout, log=log) as connection:
                message = { 'request_type' : 'record_info', 'specification' : arguments.specification, 'keywords' : arguments.keywords, 'segments' : arguments.segments, 'links' : arguments.links, 'number_records' : arguments.number_records }
                response = send_request(message, connection, log)
                message = { 'request_type' : 'quit' }
                send_request(message, connection, log)

            # message is raw JSON from jsoc_info
            record_set_info_dict = json_loads(response)
            response = get_response(record_set_info_dict, log)
        except ExpServerBaseError as exc:
            raise ExportServerError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')
        except Exception as exc:
            raise ExportServerError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')
    except ExpServerBaseError as exc:
        response = exc.response
        error_message = exc.message

        if log:
            log.write_error([ error_message ])
        else:
            print(error_message)
    except Exception as exc:
        error_message = str(exc)

        if log:
            log.write_error([ error_message ])
        else:
            print(error_message)

    if log is not None:
        log.write_info([ f'[ perform_action ] request complete; status {response.status_code.description()}' ])

    return response

# for use in export web app
from action import Action
from parse_specification import ParseSpecificationAction

@lru_cache
def get_parsed_specification(specification, db_host, webserver):
    action = Action.action(action_type='parse_specification', args={ 'specification' : specification, 'db_host' : db_host, 'webserver' : webserver, 'logging_level' : GetRecordInfoAction.logging_level })
    action.log = GetRecordInfoAction._log
    parsed_specification = action()

    return parsed_specification

class GetRecordInfoAction(Action):
    actions = [ 'get_record_set_info' ]

    _log = None
    logging_level = None

    def __init__(self, *, method, specification, db_host, drms_client_type=None, drms_client=None, keywords=None, links=None, logging_level=None, db_name=None, number_records=None, db_port=None, segments=None, db_user=None, webserver=None):
        self._method = getattr(self, method)
        self._specification = specification
        self._db_host = db_host # the host `webserver` uses (private webserver uses private db host)
        self._options = {}
        self._options['drms_client_type'] = drms_client_type
        self._options['drms_client'] = drms_client
        self._options['keywords'] = keywords # list
        self._options['links'] = links # list
        self._options['logging_level'] = logging_level
        self._options['db_name'] = db_name
        self._options['number_records'] = number_records # int
        self._options['db_port'] = db_port # int
        self._options['segments'] = segments # list
        self._options['db_user'] = db_user
        self._options['webserver'] = webserver # host name - gets converted to object in `get_arguments()`

    def get_record_set_info(self):
        response = perform_action(action_obj=self, is_program=False, specification=self._specification, db_host=self._db_host, options=self._options)
        return response

    @property
    def log(self):
        return self.__class__._log

    @log.setter
    def log(self, log):
        self.__class__._log = log

    @classmethod
    def get_parsed_specification(cls, specification, db_host, webserver):
        parsed_specification = get_parsed_specification(specification, db_host, webserver)

        if isinstance(parsed_specification, ErrorResponse) or parsed_specification is None:
            raise ArgumentsError(error_message=f'unable to parse specification `{specification}`')

        return parsed_specification

    @classmethod
    def is_valid_specification(cls, specification, db_host, webserver, logging_level=None):
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
            cls._logging_level = logging_level
            response = cls.get_parsed_specification(specification, db_host_resolved, webserver)

            is_valid = False if isinstance(response, ErrorResponse) else True
        except:
            is_valid = False

        return is_valid

if __name__ == "__main__":
    try:
        response = perform_action(action_obj=None, is_program=True)
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
