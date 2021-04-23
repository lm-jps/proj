#!/usr/bin/env python3

# This is a version of checkExpDbServer.py modified to run inside a Flask context.

# Returns { "server" : "hmidb2", "series" : [{ "hmi.M_45s" : { "server" : "hmidb" } }, { "hmi.internalonly" : { "server" : "hmidb2" }}]}

from __future__ import print_function

import sys
import os
import pwd
from subprocess import check_output, check_call, CalledProcessError, STDOUT
import json
from drmsCmdl import CmdlParser
from drmsparams import DRMSParams
from statuscode import StatusCode as SC
from arguments import Arguments as Args
from response import Response
from error import Error as ExportError
import securedrms

class StatusCode(SC):
    SUCCESS = 0, 'success'

class ErrorCode(SC):
    PARAMETERS = 1, 'failure locating DRMS parameters'
    ARGUMENTS = 2, 'bad arguments'
    WHITELIST = 3, 'whitelists are unsupported'
    SERIES_INFO = 4, 'unable to obtain series information'

class ParametersError(ExportError):
    _status_code = StatusCode(ErrorCode.PARAMETERS)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class ArgumentsError(ExportError):
    _status_code = StatusCode(ErrorCode.ARGUMENTS)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class WhitelistError(ExportError):
    _status_code = StatusCode(ErrorCode.WHITELIST)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class SeriesInfoError(ExportError):
    _status_code = StatusCode(ErrorCode.SERIES_INFO)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class ValidateArgumentAction(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        # the server specified must not be the internal server
        if self.dest == 'wlfile' && not parser.drmsparams.WL_HASWL:
            raise ArgumentsError(f'{option_string} specified a white-list file, but this DRMS does not support series whitelists')

        if value.lower() == parser.drmsparams.SERVER.lower():
            raise ArgumentsError(f'{option_string} specified the internal server, but you must specify an external server')

        setattr(namespace, self.dest, value)

class Arguments(Args):
    _arguments = None

    @classmethod
    def get_arguments(cls, *, program_args, drms_params):
        if cls._arguments is None:
            try:
                dbport = int(drms_params.get_required('DRMSPGPORT'))
                dbname = drms_params.get_required('DBNAME')
                dbuser = drms_params.get_required('WEB_DBUSER')
                wl_file = drms_params.get_required('WL_FILE')

                internal_dbhost = drms_params.get_required('SERVER')
                has_wl = drms_params.get_required('WL_HASWL')
            except DPMissingParameterError as exc:
                raise ParametersError(msg=str(exc))

            args = None

            if program_args is not None and len(program_args) > 0:
                args = program_args

            parser = CmdlParser(usage='%(prog)s [ -dhn ] dbhost=<db host> series=<DRMS series list> [ --dbport=<db port> ] [ --dbname=<db name> ] [ --dbuser=<db user> ] [--wlfile=<white-list text file> ]')

            parser.drmsparams = drmsparams

            # Required
            parser.add_argument('H', 'dbhost', help='The machine hosting the EXTERNAL database that serves DRMS data series names.', metavar='<db host>', dest='dbhost', required=True, action=ValidateArgumentAction)
            parser.add_argument('s', 'series', help='A comma-separated list of series to be checked.', metavar='<series>', dest='series', required=True)

            # Optional
            parser.add_argument('-d', '--debug', help='Run in CGI mode, and print helpful diagnostics.).', dest='debug', action='store_true', default=False)
            parser.add_argument('-P', '--dbport', help='The port on the machine hosting DRMS data series names.', metavar='<db host port>', dest='dbport', default=dbport)
            parser.add_argument('-N', '--dbname', help='The name of the database serving DRMS series names.', metavar='<db name>', dest='dbname', default=dbname)
            parser.add_argument('-U', '--dbuser', help='The user to log-in to the serving database as.', metavar='<db user>', dest='dbuser', default=dbuser)
            parser.add_argument('-w', '--wlfile', help='The text file containing the definitive list of internal series accessible via the external web site.', metavar='<white-list file>', dest='wlfile', action=ValidateArgumentAction, default=wl_file)

            cls._arguments = Arguments(parser=parser, args=args)
            cls._arguments.set_arg('internal_dbhost', internal_dbhost)
            cls._arguments.set_arg('has_wl', has_wl)

        return cls._arguments



# Return codes for cmd-line run.
RET_SUCCESS = 0
RET_BADARGS = 1
RET_DRMSPARAMS = 2
RET_WHITELIST = 4
RET_ARCH = 5
RET_SHOWSERIES = 6


rv = RET_SUCCESS

# Parse arguments
if __name__ == "__main__":
    try:
        arguments = get_arguments(kwargs)

        if arguments.operation == 'register':
            response = check_and_register()
        else:
            response = check()




    server = 'Unknown' # This is the database host to use when querying about all the series provided as the series-list argument.
    seriesObjs = [] # A list of series objects with info about the series (like the db server for that series

    try:

        # Call show_series on the external host. If all series are on the external host, return the external server optD['dbhost']. show_series
        # does not provide a way to search for a list of series, so we have to make a hash of all series, then assess each series for membsership in the list
        # of external series.
        #
        # There is no JSOC_DBPORT DRMS parameter, much to my surprise. We need to append ':' + str(optD['dbport']) to optD['dbhost'].

        # XXX user securedrms.py to call, on external server, show_series -qz f'JSOC_DBHOST={arguments.dbhost}:str({arguments.dbport})' JSOC_DBNAME=arguments.dbname JSOC_DBUSER=arguments.dbuser; obtain json_response
        if json_response is not None:
            dict_response = json.loads(json_response)
        else:
            raise SeriesInfoError(f'unexpected response from show_series')

        external_map = {}
        internal_map = {}
        not_external = []
        passthru_series = {}

        for series in dict_response['names']:
            external_map[str(series['name']).strip().lower()] = str(series['name']).strip()

        for series in arguments.series:
            if series.lower() not in external_map:
                not_external.append(series)

        # if all series are external, then we can return the external database server
        if len(not_external) == 0:
            server = arguments.dbhost

           # seriesObjs - [{ "hmi.tdpixlist" : { "server" : "hmidb" } }, { "hmi.internalonly" : { "server" : "hmidb2" }}, ...]
            for series in arguments.series:
                sobj = { external_map[series.lower()] : {} }
                sobj[external_map[series.lower()]]['server'] = arguments.dbhost
                series_objs.append(sobj)
        else:
            # some series may be on the internal db server

            # XXX user securedrms.py to call, on internal server, show_series -qz f'JSOC_DBHOST={arguments.dbhost}:str({arguments.dbport})' JSOC_DBNAME=arguments.dbname JSOC_DBUSER=arguments.dbuser; obtain json_response
            if json_response is not None:
                dict_response = json.loads(json_response)
            else:
                raise SeriesInfoError(f'unexpected response from show_series')

            # Hash all series in the internal DB.
            for series in dict_response['names']:
                # Remove various whitespace too.
                internal_map[str(series['name']).strip().lower()] = str(series['name']).strip()

            # Fetch whitelist.
            with open(arguments.wlfile', 'r') as whitelist:
                # NOTE: This script does not attempt to validate the series in the whitelist - there could be invalid entries in that
                # file. Series from the whitelist that match a series in internal-DB series are returned to the caller.
                for series in whitelist:
                    if series.strip().lower() in internal_map:
                        passthru_series[series.strip().lower()] = True

            # Finally, check to see if all the series in not_external are in passthru_series.
            for series in not_external:
                if series.lower() not in passthru_series:
                    raise Exception('getArgs', 'Series ' + series + ' is not a valid series accessible from ' + optD['dbhost'] + '.', RET_BADARGS)

            # series_objs - [{ "hmi.tdpixlist" : { "server" : "hmidb" } }, { "hmi.internalonly" : { "server" : "hmidb2" }}, ...]
            for series in arguments.series:
                if series.lower() in external_map:
                    sobj = { external_map[series.lower()] : {} }
                    sobj[external_map[series.lower()]]['server'] = arguments.dbhost
                else:
                    sobj = { internal_map[series.lower()] : {} }
                    sobj[internal_map[series.lower()]]['server'] = arguments.internal_dbhost
                series_objs.append(sobj)

        # status_code is 0
        # Returns - { "server" : "hmidb2", "series" : [{ "hmi.M_45s" : { "server" : "hmidb" } }, { "hmi.internalonly" : { "server" : "hmidb2" }}], "err" : 0}
        response = Response.generate_response(status_code=StatusCode.SUCCESS, msg=StatusCode.SUCCESS.fullname, server=server, series=series_objs)
    except ExportError as exc:
        response = exc.response

    print(response.generate_json())

    # Always return 0. If there was an error, an error code (the 'status' property) and message (the 'statusMsg' property) goes in the returned HTML.
    sys.exit(0)
else:
    pass
