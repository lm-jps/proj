#!/usr/bin/env python3

# This is a version of checkExpDbServer.py modified to run inside a Flask context.

# check_dbserver.py determines if series are accessible publicly;

# returns a JSON object where the `server` attribute identifies the DB host that can serve data-series information for all
# series provided in the `series` program argument - if the server identified by the `public_dbhost` program argument
# can serve all series, then the value of `server` is the value of `public_dbhost`; otherwise the value of `server` is
# the internal DB host; in the returned JSON object, the `series` attribute contains a list of objects, one for each series
# provided in the `series` program argument; each object contains a single attribute - the name of the series; the value of
# this attribute is another object which contains a single attribute, `server`, whose value is the value of the
# `public_dbhost` program argument if the series is DIRECTLY accessible from the public server, or whose value is the
# private server if the series is accessible from the private server only; `server` is None if the series does
# not exist, or if the series is accessible from the private server only BUT not on the whitelist; return examples:
#   { "server" : "hmidb", "series" : [{ "hmi.M_45s" : { "server" : "hmidb2" } }, { "hmi.on_white_list" : { "server" : "hmidb" }}], "status" : 0 }
#   { "server" : "hmidb2", "series" : [{ "hmi.M_45s" : { "server" : "hmidb2" } }, { "hmi.not_on_white_list" : { "server" : None }}], "status" : 0 }
#   { "server" : "hmidb2", "series" : [{ "hmi.M_45s" : { "server" : "hmidb2" } }, { "hmi.does_not_exist" : { "server" : None }}], "status" : 0 }

from argparse import Action as ArgsAction
from json import loads as json_loads
from sys import exc_info as sys_exc_info, exit as sys_exit

from drms_parameters import DRMSParams
from drms_utils import Arguments as Args, ArgumentsError as ArgsError, CmdlParser as ArgsParser, MakeObject, StatusCode as ExportStatusCode
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
        if self.dest == 'wlfile' and not parser.drms_params.WL_HASWL:
            raise ArgumentsError(error_message=f'{option_string} specified a white-list file, but this DRMS does not support series whitelists')

        if value.lower() == parser.drms_params.SERVER.lower():
            raise ArgumentsError(error_message=f'{option_string} specified the internal server, but you must specify an external server')

        setattr(namespace, self.dest, value)

class SeriesAction(ArgsAction):
    def __call__(self, parser, namespace, value, option_string=None):
        series_dict = json_loads(value)
        setattr(namespace, self.dest, series_dict)

class Arguments(Args):
    _arguments = None

    @classmethod
    def get_arguments(cls, *, is_program, program_name=None, program_args=None, module_args=None, drms_params):
        if cls._arguments is None:
            try:
                db_port = int(drms_params.get_required('DRMSPGPORT'))
                db_name = drms_params.get_required('DBNAME')
                db_user = drms_params.get_required('WEB_DBUSER')
                wl_file = drms_params.get_required('WL_FILE')

                private_dbhost = drms_params.get_required('SERVER')
                has_wl = drms_params.get_required('WL_HASWL')
            except DPMissingParameterError as exc:
                raise ParametersError(error_message=str(exc))

            if is_program:
                args = None

                if program_args is not None and len(program_args) > 0:
                    args = program_args

                if program_name is not None and len(program_name) > 0:
                    parser = ArgsParser(prog=program_name, usage='%(prog)s [ -dhn ] public_dbhost=<db host> series=<DRMS series list> [ -cd ] [ --dbport=<db port> ] [ --dbname=<db name> ] [ --dbuser=<db user> ] [--wlfile=<white-list text file> ]')
                else:
                    parser = ArgsParser(usage='%(prog)s [ -dhn ] public_dbhost=<db host> series=<DRMS series list> [ -cd ] [ --dbport=<db port> ] [ --dbname=<db name> ] [ --dbuser=<db user> ] [--wlfile=<white-list text file> ]')

                # to give the parser access to a few drms parameters inside ValidateArgumentAction
                parser.drms_params = drms_params

                # Required
                parser.add_argument('public_dbhost', help='the machine hosting the EXTERNAL database that serves DRMS data series names.', metavar='<db host>', action=ValidateArgumentAction, dest='public_dbhost', required=True)
                parser.add_argument('s', 'series', help='a json string containing a list of series to be checked (`{ "series" : [ "series1", "series2", ... ] }`)', metavar='<series>', action=SeriesAction, dest='series', required=True)

                # Optional
                parser.add_argument('-c', '--drms-client-type', help='securedrms client type (ssh, http)', choices=[ 'ssh', 'http' ], dest='drms_client_type', default='ssh')
                parser.add_argument('-d', '--debug', help='print helpful securedrms diagnostics', dest='debug', action='store_true', default=False)
                parser.add_argument('-P', '--dbport', help='the port on the machine hosting DRMS data series names', metavar='<db host port>', dest='db_port', default=db_port)
                parser.add_argument('-N', '--dbname', help='the name of the database serving DRMS series names', metavar='<db name>', dest='db_name', default=db_name)
                parser.add_argument('-U', '--dbuser', help='the user to log-in to the serving database as', metavar='<db user>', dest='db_user', default=db_user)
                parser.add_argument('-w', '--wlfile', help='the text file containing the definitive list of internal series accessible via the external web site', metavar='<white-list file>', dest='wlfile', action=ValidateArgumentAction, default=wl_file)

                arguments = Arguments(parser=parser, args=args)

                arguments.drms_client = None
            else:
                # `program_args` has all `arguments` values, in final form; validate them
                def extract_module_args(*, public_dbhost, series, drms_client=None, drms_client_type='ssh', debug=False, db_port=db_port, db_name=db_name, db_user=db_user, wl_file=wl_file):
                    arguments = {}

                    arguments['public_dbhost'] = public_dbhost
                    arguments['series'] = series
                    arguments['drms_client'] = drms_client
                    arguments['drms_client_type'] = drms_client_type
                    arguments['debug'] = debug
                    arguments['db_port'] = db_port
                    arguments['db_name'] = db_name
                    arguments['db_user'] = db_user
                    arguments['wl_file'] = wl_file

                    return arguments

                # dict
                module_args_dict = extract_module_args(**module_args)
                arguments = Arguments(parser=None, args=module_args_dict)

            arguments.private_dbhost = private_dbhost
            arguments.has_wl = has_wl

            cls._arguments = arguments

        return cls._arguments


# for use in export web app
from action import Action
class DetermineDbServerAction(Action):
    actions = [ 'determine_db_server' ]
    def __init__(self, *, method, public_dbhost, series, drms_client=None, db_port=None, db_name=None, db_user=None):
        self._method = getattr(self, method)
        self._public_dbhost = public_dbhost
        self._series = series # dict with "series" : [ ... ]
        self._options = {}
        self._options['drms_client'] = drms_client

        if db_port is not None:
            self._options['db_port'] = db_port

        if db_name is not None:
            self._options['db_name'] = db_name

        if db_user is not None:
            self._options['db_user'] = db_user

    def determine_db_server(self):
        # returns dict
        response = perform_action(is_program=False, public_dbhost=self._public_dbhost, series=self._series, options=self._options)
        return response.generate_dict()

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

            # call show_series on the external host; if all series are on the external host, return the external server `arguments.dbhost`; show_series
            # does not provide a way to search for a list of series, so make a hash of all series, then assess each series for membsership in the list
            #
            # There is no JSOC_DBPORT DRMS parameter, much to my surprise. We need to append ':' + str(optD['dbport']) to optD['dbhost'].

            # XXX use securedrms.py to call, on external server, show_series -qz f'JSOC_DBHOST={arguments.dbhost}:str({arguments.dbport})' JSOC_DBNAME=arguments.dbname JSOC_DBUSER=arguments.dbuser; obtain json_response
            if arguments.drms_client is None:
                connection_info = { 'dbhost' : arguments.public_dbhost, 'dbport' : arguments.db_port, 'dbname' : arguments.db_name, 'dbuser' : arguments.db_user }

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

                connection_info = { 'dbhost' : arguments.private_dbhost, 'dbport' : arguments.db_port, 'dbname' : arguments.db_name, 'dbuser' : arguments.db_user }
                if arguments.drms_client_type == 'ssh':
                    sshclient_internal = factory.create_client(server='jsoc_internal', use_ssh=True, use_internal=True, connection_info=connection_info)
                else:
                    sshclient_internal = factory.create_client(server='jsoc_internal', use_ssh=False, use_internal=True, connection_info=connection_info)

            # securedrms configuration does not include ssh_show_series_wrapper, so the results will include only series that are implemented in the public database (i.e., no "pass-through" series)
            response_dict = {}

            if len(arguments.series['series']) > 0:
                response_dict['series'] = []
                use_public_server = True
                white_list = None
                supporting_server = False

                for series in arguments.series['series']:
                    series_regex = series.lower().replace('.', '[.]')
                    public_series = sshclient_external.series(regex=series_regex, full=False)
                    if len(public_series) == 0:
                        if arguments.has_wl:
                            # try private server
                            private_series = sshclient_internal.series(regex=series_regex, full=False)

                            if white_list is None:
                                white_list = get_whitelist(arguments.wlfile) # if no whitelist exists, then this is the empty set

                            if series.lower() in white_list and len(private_series) != 0:
                                response_dict['series'].append({ series : { 'server' : arguments.private_dbhost } })
                                supporting_server = True
                                use_public_server = False
                            else:
                                response_dict['series'].append({ series : { 'server' : None } })
                        else:
                            response_dict['series'].append({ series : { 'server' : None } })
                    else:
                        response_dict['series'].append({ series : { 'server' : arguments.public_dbhost } })
                        supporting_server = True

                if supporting_server:
                    response_dict['server'] = arguments.public_dbhost if use_public_server else arguments.private_dbhost
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
