#!/usr/bin/env python3

# This is a version of checkExpDbServer.py modified to run inside a Flask context.

# check_dbserver.py determines if series are accessible publicly;

# returns a JSON object where the `server` attribute identifies the DB host that can serve data-series information for all
# series provided in the `series` program argument - if the server identified by the `public_db_host` program argument
# can serve all series, then the value of `server` is the value of `public_db_host`; otherwise the value of `server` is
# the internal DB host; in the returned JSON object, the `series` attribute contains a list of objects, one for each series
# provided in the `series` program argument; each object contains a single attribute - the name of the series; the value of
# this attribute is another object which contains a single attribute, `server`, whose value is the value of the
# `public_db_host` program argument if the series is DIRECTLY accessible from the public server, or whose value is the
# private server if the series is accessible from the private server only; `server` is None if the series does
# not exist, or if the series is accessible from the private server only BUT not on the whitelist; return examples:
#   { "server" : "hmidb", "series" : [{ "hmi.M_45s" : { "server" : "hmidb2" } }, { "hmi.on_white_list" : { "server" : "hmidb" }}], "status" : 0 }
#   { "server" : "hmidb2", "series" : [{ "hmi.M_45s" : { "server" : "hmidb2" } }, { "hmi.not_on_white_list" : { "server" : None }}], "status" : 0 }
#   { "server" : "hmidb2", "series" : [{ "hmi.M_45s" : { "server" : "hmidb2" } }, { "hmi.does_not_exist" : { "server" : None }}], "status" : 0 }

from argparse import Action as ArgsAction
from json import loads as json_loads
from sys import exc_info as sys_exc_info, exit as sys_exit

from drms_parameters import DRMSParams, DPMissingParameterError
from drms_utils import Arguments as Args, ArgumentsError as ArgsError, CmdlParser, MakeObject, StatusCode as ExportStatusCode
from drms_export import Response, Error as ExportError, ErrorCode as ExportErrorCode
from drms_export import securedrms
from utils import extract_program_and_module_args

class StatusCode(ExportStatusCode):
    SUCCESS = (0, 'success')

class ErrorCode(ExportErrorCode):
    PARAMETERS = (1, 'failure locating DRMS parameters')
    ARGUMENTS = (2, 'bad arguments')
    WHITELIST = (3, 'whitelists are unsupported')
    SERIES_INFO = (4, 'unable to obtain series information')
    SECURE_DRMS = (5, 'failure with securedrms interface')

class CdbBaseError(ExportError):
    def __init__(self, *, exc_info=None, error_message=None):
        if exc_info is not None:
            self.exc_info = exc_info
            e_type, e_obj, e_tb = exc_info

            if error_message is None:
                error_message = f'{e_type.__name__}: {str(e_obj)}'
            else:
                error_message = f'{error_message} [ {e_type.__name__}: {str(e_obj)} ]'

        super().__init__(error_message=error_message)

class ParametersError(CdbBaseError):
    _error_code = ErrorCode.PARAMETERS

class ArgumentsError(CdbBaseError):
    _error_code = ErrorCode.ARGUMENTS

class WhitelistError(CdbBaseError):
    _error_code = ErrorCode.WHITELIST

class SeriesInfoError(CdbBaseError):
    _error_code = ErrorCode.SERIES_INFO

class SecureDRMSError(CdbBaseError):
    _error_code = ErrorCode.SECURE_DRMS

class ValidateArgumentAction(ArgsAction):
    def __call__(self, parser, namespace, value, option_string=None):
        # the server specified must not be the internal server
        if self.dest == 'wl_file' and not parser.drms_params.WL_HASWL:
            raise ArgumentsError(error_message=f'{option_string} specified a white-list file, but this DRMS does not support series whitelists')

        if value.lower() == parser.drms_params.SERVER.lower():
            raise ArgumentsError(error_message=f'{option_string} specified the internal server, but you must specify an external server')

        setattr(namespace, self.dest, value)

class Arguments(Args):
    _arguments = None

    @classmethod
    def get_arguments(cls, *, is_program, program_name=None, program_args=None, module_args=None, drms_params, refresh=True):
        if cls._arguments is None or refresh:
            try:
                db_port = int(drms_params.get_required('DRMSPGPORT'))
                db_name = drms_params.get_required('DBNAME')
                db_user = drms_params.get_required('WEB_DBUSER')
                wl_file = drms_params.get_required('WL_FILE')

                private_db_host = drms_params.get_required('SERVER')
                has_wl = drms_params.get_required('WL_HASWL')
            except DPMissingParameterError as exc:
                raise ParametersError(error_message=str(exc))

            if is_program:
                args = None

                if program_args is not None and len(program_args) > 0:
                    args = program_args

                parser_args = { 'usage' : '%(prog)s public-db-host=<db host> series=<DRMS series list> [ -c/--drms-client-type=<ssh/http> ] [ -d/--debug ] [ -/P--dbport=<db port> ] [ -N/--dbname=<db name> ] [ -U/--dbuser=<db user> ] [ -w/--wlfile=<white-list text file> ]' }

                if program_name is not None and len(program_name) > 0:
                    parser_args['prog'] = program_name

                parser = CmdlParser(**parser_args)

                # to give the parser access to a few drms parameters inside ValidateArgumentAction
                parser.drms_params = drms_params

                # Required
                parser.add_argument('public-db-host', help='the machine hosting the EXTERNAL database that serves DRMS data series names.', metavar='<db host>', action=ValidateArgumentAction, dest='public_db_host', required=True)
                parser.add_argument('s', 'series', help='a comma-separated ist of series to be checked', metavar='<series>', action=ListAction, dest='series', required=True)

                # Optional
                parser.add_argument('-c', '--drms-client-type', help='securedrms client type (ssh, http)', choices=[ 'ssh', 'http' ], dest='drms_client_type', default='ssh')
                parser.add_argument('-d', '--debug', help='print helpful securedrms diagnostics', dest='debug', action='store_true', default=False)
                parser.add_argument('-P', '--dbport', help='the port on the machine hosting DRMS data series names', metavar='<db host port>', dest='db_port', default=db_port)
                parser.add_argument('-N', '--dbname', help='the name of the database serving DRMS series names', metavar='<db name>', dest='db_name', default=db_name)
                parser.add_argument('-U', '--dbuser', help='the user to log-in to the serving database as', metavar='<db user>', dest='db_user', default=db_user)
                parser.add_argument('-w', '--wlfile', help='the text file containing the definitive list of internal series accessible via the external web site', metavar='<white-list file>', dest='wl_file', action=ValidateArgumentAction, default=wl_file)

                arguments = Arguments(parser=parser, args=args)

                arguments.drms_client = None
            else:
                # `program_args` has all `arguments` values, in final form; validate them
                # `series` is py list of DRMS data series
                def extract_module_args(*, public_db_host, series, drms_client_type='ssh', drms_client=None, debug=False, db_port=db_port, db_name=db_name, db_user=db_user, wl_file=wl_file):
                    arguments = {}

                    arguments['public_db_host'] = public_db_host
                    arguments['series'] = series # list
                    arguments['drms_client_type'] = drms_client_type
                    arguments['drms_client'] = drms_client
                    arguments['debug'] = debug
                    arguments['db_port'] = db_port
                    arguments['db_name'] = db_name
                    arguments['db_user'] = db_user
                    arguments['wl_file'] = wl_file

                    return arguments

                # dict
                module_args_dict = extract_module_args(**module_args)
                arguments = Arguments(parser=None, args=module_args_dict)

            arguments.private_db_host = private_db_host
            arguments.has_wl = has_wl

            cls._arguments = arguments

        return cls._arguments


# for use in export web app
from action import Action
from get_series_info import GetSeriesInfoAction
class DetermineDbServerAction(Action):
    actions = [ 'determine_db_server' ]
    def __init__(self, *, method, public_db_host, series, drms_client_type=None, drms_client=None, db_port=None, db_name=None, db_user=None):
        self._method = getattr(self, method)
        self._public_db_host = public_db_host
        self._series = series # py list
        self._options = {}
        self._options['drms_client_type'] = drms_client_type
        self._options['drms_client'] = drms_client
        self._options['db_port'] = db_port
        self._options['db_name'] = db_name
        self._options['db_user'] = db_user

    def determine_db_server(self):
        # returns dict
        response = perform_action(is_program=False, public_db_host=self._public_db_host, series=self._series, options=self._options)
        return response

    # `series` is comma-separated list of series
    @classmethod
    def is_valid_series_set(cls, series, db_host, webserver, logging_level=None):
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

            return GetSeriesInfoAction.is_valid_series_set(series, db_host_resolved, webserver, logging_level)
        except:
            return False

def get_whitelist(wl_file):
    white_list = set()

    with open(wl_file, 'r') as f_white_list:
        # NOTE: This script does not attempt to validate the series in the whitelist - there could be invalid entries in that
        # file. Series from the whitelist that match a series in internal-DB series are returned to the caller.
        for series in f_white_list:
            white_list.add(series.strip().lower())

    return white_list

def perform_action(is_program, program_name=None, **kwargs):
    response = None

    try:
        program_args, module_args = extract_program_and_module_args(is_program=is_program, **kwargs)

        try:
            drms_params = DRMSParams()

            if drms_params is None:
                raise ParameterError(error_message='unable to locate DRMS parameters file (drmsparams.py)')

            arguments = Arguments.get_arguments(is_program=is_program, program_name=program_name, program_args=program_args, module_args=module_args, drms_params=drms_params)
        except ArgsError as exc:
            raise ArgumentsError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')
        except Exception as exc:
            raise ArgumentsError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

        try:
            factory = None

            # call show_series on the external host; if all series are on the external host, return the external server `arguments.db_host`; show_series does not provide a way to search for a list of series, so make a hash of all series, then assess each series for membsership in the list
            #
            # There is no JSOC_DBPORT DRMS parameter, much to my surprise. We need to append ':' + str(optD['dbport']) to optD['db_host'].

            # XXX use securedrms.py to call, on external server, show_series -qz f'JSOC_DBHOST={arguments.db_host}:str({arguments.dbport})' JSOC_DBNAME=arguments.dbname JSOC_DBUSER=arguments.dbuser; obtain json_response
            if arguments.drms_client is None:
                connection_info = { 'dbhost' : arguments.public_db_host, 'dbport' : arguments.db_port, 'dbname' : arguments.db_name, 'dbuser' : arguments.db_user }

                factory = securedrms.SecureClientFactory(debug=arguments.debug)
                if arguments.drms_client_type == 'ssh':
                    sshclient_external = factory.create_client(server='jsoc_external', use_ssh=True, use_internal=False, connection_info=connection_info)
                else:
                    sshclient_external = factory.create_client(server='jsoc_external', use_ssh=False, use_internal=False, connection_info=connection_info)
            else:
                # use client passed in from flask app - must be a public client
                if arguments.drms_client.use_internal:
                    raise SecureDRMSError(error_message=f'must provide a securedrms client with public access')

                sshclient_external = arguments.drms_client

            # we also need an internal client for pass-through series, if the servers support pass-through series
            if arguments.has_wl:
                if factory is None:
                    factory = securedrms.SecureClientFactory(debug=arguments.debug)

                connection_info = { 'dbhost' : arguments.private_db_host, 'dbport' : arguments.db_port, 'dbname' : arguments.db_name, 'dbuser' : arguments.db_user }
                if arguments.drms_client_type == 'ssh':
                    sshclient_internal = factory.create_client(server='jsoc_internal', use_ssh=True, use_internal=True, connection_info=connection_info)
                else:
                    sshclient_internal = factory.create_client(server='jsoc_internal', use_ssh=False, use_internal=True, connection_info=connection_info)

            # securedrms configuration does not include ssh_show_series_wrapper, so the results will include only series that are implemented in the public database (i.e., no "pass-through" series)
            response_dict = {}

            if len(arguments.series) > 0:
                response_dict['series'] = []
                use_public_server = True
                white_list = None
                supporting_server = False

                for series in arguments.series:
                    series_regex = series.lower().replace('.', '[.]')
                    public_series = sshclient_external.series(regex=series_regex, full=False)
                    if len(public_series) == 0:
                        if arguments.has_wl:
                            # try private server
                            private_series = sshclient_internal.series(regex=series_regex, full=False)

                            if white_list is None:
                                white_list = get_whitelist(arguments.wl_file) # if no whitelist exists, then this is the empty set

                            if series.lower() in white_list and len(private_series) != 0:
                                response_dict['series'].append({ series : { 'server' : arguments.private_db_host } })
                                supporting_server = True
                                use_public_server = False
                            else:
                                response_dict['series'].append({ series : { 'server' : None } })
                        else:
                            response_dict['series'].append({ series : { 'server' : None } })
                    else:
                        response_dict['series'].append({ series : { 'server' : arguments.public_db_host } })
                        supporting_server = True

                if supporting_server:
                    response_dict['server'] = arguments.public_db_host if use_public_server else arguments.private_db_host
                else:
                    response_dict['server'] = None
            else:
                response_dict['server'] = None
                response_dict['series'] = []

            # response examples:
            #   { "server" : "hmidb", "series" : [{ "hmi.M_45s" : { "server" : "hmidb2" } }, { "hmi.on_white_list" : { "server" : "hmidb" }}], "status" : 0 }
            #   { "server" : "hmidb2", "series" : [{ "hmi.M_45s" : { "server" : "hmidb2" } }, { "hmi.not_on_white_list" : { "server" : None }}], "status" : 0 }
            #   { "server" : "hmidb2", "series" : [{ "hmi.M_45s" : { "server" : "hmidb2" } }, { "hmi.does_not_exist" : { "server" : None }}], "status" : 0 }
            response = Response.generate_response(status_code=StatusCode.SUCCESS, **response_dict)
        except securedrms.SecureDRMSError as exc:
            response = SecureDRMSError(error_message=str(exc)).response
    except ExportError as exc:
        response = exc.response
        e_type, e_obj, e_tb = exc.exc_info
        error_msg = f'ERROR LINE {str(e_tb.tb_lineno)}: {exc.message}'

        print(f'{error_msg}')

    return response


# Parse arguments
if __name__ == "__main__":
    try:
        response = perform_action(is_program=True)
    except ExportError as exc:
        response = exc.response

    print(response.generate_json())

    # Always return 0. If there was an error, an error code (the 'status' property) and message (the 'statusMsg' property) goes in the returned HTML.
    sys_exit(0)
else:
    pass
