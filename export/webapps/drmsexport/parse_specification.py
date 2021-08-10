#!/usr/bin/env python3

from drms_export import Error as ExportError, ErrorCode as ExportErrorCode, Response, securedrms
from drms_parameters import DRMSParams, DPMissingParameterError
from drms_utils import Arguments as Args, Choices, CmdlParser, Formatter as DrmsLogFormatter, Log as DrmsLog, LogLevel as DrmsLogLevel, LogLevelAction as DrmsLogLevelAction, MakeObject, StatusCode as SC
import inspect
from json import loads as json_loads, decoder as json_decoder
from os import environ
from os.path import join as path_join
from sys import exc_info as sys_exc_info, exit as sys_exit

DEFAULT_LOG_FILE = 'ps_log.txt'
PARSE_SPEC_BIN = 'drms_parserecset'

class StatusCode(SC):
    SUCCESS = 0, 'success'

class ErrorCode(ExportErrorCode):
    PARAMETERS = 1, 'failure locating DRMS parameters'
    ARGUMENTS = 2, 'bad arguments'
    LOGGING = 3, 'failure logging messages'
    DRMS_CLIENT = 4, 'drms client error'

class PsBaseError(ExportError):
    def __init__(self, *, exc_info=None, error_message=None):
        if exc_info is not None:
            self.exc_info = exc_info
            e_type, e_obj, e_tb = exc_info

            if error_message is None:
                error_message = f'{e_type.__name__}: {str(e_obj)}'

        super().__init__(error_message=error_message)

class ParametersError(PsBaseError):
    _error_code = ErrorCode(ErrorCode.PARAMETERS)

    def __init__(self, *, exc_info=None, error_message=None):
        super().__init__(exc_info=exc_info, error_message=error_message)

class ArgumentsError(PsBaseError):
    _error_code = ErrorCode(ErrorCode.ARGUMENTS)

    def __init__(self, *, exc_info=None, error_message=None):
        super().__init__(exc_info=exc_info, error_message=error_message)

class LoggingError(PsBaseError):
    _error_code = ErrorCode(ErrorCode.LOGGING)

    def __init__(self, *, exc_info=None, error_message=None):
        super().__init__(exc_info=exc_info, error_message=error_message)

class DRMSClientError(PsBaseError):
    _error_code = ErrorCode(ErrorCode.DRMS_CLIENT)

    def __init__(self, *, exc_info=None, error_message=None):
        super().__init__(exc_info=exc_info, error_message=error_message)

class Arguments(Args):
    _arguments = None

    @classmethod
    def get_arguments(cls, *, is_program, program_name=None, program_args, drms_params):
        if cls._arguments is None:
            try:
                log_file = path_join(drms_params.get_required('EXPORT_LOG_DIR'), DEFAULT_LOG_FILE)
                export_bin = drms_params.get_required('BIN_EXPORT')
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

                if program_name is not None and len(program_name) > 0:
                    parser = CmdlParser(prog=program_name, usage=f'%(prog)s spec=<specification> [ --dbhost=<db host> ] [ --dbport=<db port> ] [ --dbname=<db name> ] [ --dbuser=<db user>]')
                else:
                    parser = CmdlParser(usage=f'%(prog)s spec=<specification> [ --dbhost=<db host> ] [ --dbport=<db port> ] [ --dbname=<db name> ] [ --dbuser=<db user>]')

                # required
                parser.add_argument('spec', help='the DRMS record-set specification', metavar='<specification>', dest='specification', required=True)
                parser.add_argument('address', help='the email addressed registered for export', metavar='<email address>', dest='address', required=True)

                # optional
                parser.add_argument('-c', '--drms-client-type', help='securedrms client type (ssh, http)', choices=[ 'ssh', 'http' ], dest='drms_client_type', default='ssh')
                parser.add_argument('-l', '--log_file', help='the path to the log file', metavar='<log file>', dest='log_file', default=log_file)
                parser.add_argument('-L', '--logging_level', help='the amount of logging to perform; in order of increasing verbosity: critical, error, warning, info, debug', metavar='<logging level>', dest='logging_level', action=DrmsLogLevelAction, default=DrmsLogLevel.ERROR)
                parser.add_argument('-H', '--db_host', help='the machine hosting the database that contains export requests from this site', metavar='<db host>', dest='db_host', default=private_db_host)
                parser.add_argument('-P', '--db_port', help='the port on the host machine that is accepting connections for the database', metavar='<db host port>', dest='db_port', type=int, default=db_port)
                parser.add_argument('-N', '--db_name', help='the name of the database that contains export requests', metavar='<db name>', dest='db_name', default=db_name)
                parser.add_argument('-U', '--db_user', help='the name of the database user account', metavar='<db user>', dest='db_user', default=db_user)

                arguments = Arguments(parser=parser, args=args)
                arguments.drms_client = None
            else:
                def extract_module_args(*, spec, address, drms_client=None, drms_client_type='ssh', log_file=log_file, logging_level=logging_level, db_host=db_host, db_port=db_port, db_name=db_name, db_user=db_user):
                    arguments = {}

                    arguments['address'] = address
                    arguments['drms_client'] = drms_client
                    arguments['drms_client_type'] = drms_client_type
                    arguments['log_file'] = log_file
                    arguments['logging_level'] = logging_level
                    arguments['db_host'] = db_host
                    arguments['db_port'] = db_port
                    arguments['db_name'] = db_name
                    arguments['db_user'] = db_user

                    return arguments

                # dict
                module_args_dict = extract_module_args(**module_args)
                arguments = Arguments(parser=None, args=module_args_dict)

            arguments.private_db_host = private_db_host
            arguments.arch_dir = environ['JSOC_MACHINE']
            arguments.export_bin = export_bin

            cls._arguments = arguments

        return cls._arguments

# for use in export web app
from action import Action
class ParseSpecificationAction(Action):
    actions = [ 'parse_specification' ]
    def __init__(self, *, method, specification, db_host=None, db_port=None, db_name=None, db_user=None):
        self._method = getattr(self, method)
        self._specification = specification
        self._options = {}

        if db_host is not None:
            self._options['db_host'] = db_host

        if db_port is not None:
            self._options['db_port'] = db_port

        if db_name is not None:
            self._options['db_name'] = db_name

        if db_user is not None:
            self._options['db_user'] = db_user

    def parse_specification(self):
        # returns dict
        prog_name = inspect.currentframe().f_code.co_name
        response = perform_action(is_program=False, program_name=prog_name, spec=self._specification, options=self._options)
        return response.generate_dict()

def perform_action(*, is_program, program_name=None, **kwargs):
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
                raise ParametersError(error_message=f'unable to locate DRMS parameters package')

            arguments = Arguments.get_arguments(is_program=is_program, program_name=program_name, program_args=args, drms_params=drms_params)
        except ExportError as exc:
            exc.exc_info = sys_exc_info()
            raise exc
        except Exception as exc:
            raise ArgumentsError(exc_info=sys_exc_info())

        try:
            formatter = DrmsLogFormatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
            log = DrmsLog(arguments.log_file, arguments.logging_level, formatter)
        except Exception as exc:
            raise LoggingError(error_message=f'{str(exc)}')

        # run PARSE_SPEC_BIN
        try:
            if arguments.drms_client is None:
                # make a public client, since either public or private can parse record-set specifications
                debug = True if arguments.logging_level == DrmsLogLevel.DEBUG else False
                factory = securedrms.SecureClientFactory(debug=debug, email=arguments.address)
                use_ssh = True if arguments.drms_client_type == 'ssh' else False
                connection_info = { 'dbhost' : arguments.db_host, 'dbport' : arguments.db_port, 'dbname' : arguments.db_name, 'dbuser' : arguments.db_user }

                log.write_debug([ f'[ perform_action ] creating public securedrms client' ])
                public_drms_client = factory.create_client(server='jsoc_external', use_ssh=use_ssh, use_internal=False, connection_info=connection_info)
                drms_client = public_drms_client
            else:
                # could be either public or private client
                drms_client = arguments.drms_client

            response_dict = drms_client.parse_spec(arguments.specification)

            if response_dict['errMsg'] is not None:
                raise DRMSClientError(error_message=f'failure parsing record-set specification {arguments.export_arguments["specification"]}: {response_dict["errMsg"]}')
        except ExportError as exc:
            exc.exc_info = sys_exc_info()
            raise exc
        except Exception as exc:
            raise DRMSClientError(exc_info=sys_exc_info())

        response = Response.generate_response(status_code=StatusCode.SUCCESS, **response_dict)

    except ExportError as exc:
        response = exc.response
        if log:
            e_type, e_obj, e_tb = exc.exc_info
            log.write_error([ f'{e_type.__name__} (LINE {str(e_tb.tb_lineno)}): {str(e_obj)}' ])

    return response

if __name__ == '__main__':
    try:
        response = perform_action(is_program=True)
    except ExportError as exc:
        response = exc.response

    print(response.generate_json())

    sys_exit(0)
else:
    pass
