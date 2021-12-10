#!/usr/bin/env python3

from argparse import Action as ArgsAction
from copy import deepcopy
from distutils.util import strtobool
from collections import OrderedDict
from os.path import join as path_join
from psycopg2 import Error as PGError, connect as pg_connect
from sys import exc_info as sys_exc_info, exit as sys_exit

from drms_export import Error as ExportError, ErrorCode as ExportErrorCode, ErrorResponse, Response
from drms_parameters import DRMSParams, DPMissingParameterError
from drms_utils import Arguments as Args, ArgumentsError as ArgsError, Choices, CmdlParser, Formatter as DrmsLogFormatter, ListAction, Log as DrmsLog, LogLevel as DrmsLogLevel, LogLevelAction as DrmsLogLevelAction, MakeObject, StatusCode as SC
from utils import extract_program_and_module_args, get_db_host

DEFAULT_LOG_FILE = 'gsi_log.txt'

class StatusCode(SC):
    # info success status codes
    SUCCESS = (0, 'success')

class ErrorCode(ExportErrorCode):
    # info failure error codes
    FAILURE = (1, 'info failure')

    PARAMETERS = (101, 'failure locating DRMS parameters')
    ARGUMENTS = (102, 'bad arguments')
    LOGGING = (103, 'failure logging messages')
    EXPORT_ACTION = (104, 'failure calling export action')
    RESPONSE = (105, 'unable to generate valid response')
    DB_CONNECTION = (106, 'failure connecting to database')
    DB_COMMAND = (107, 'failure executing database command')
    UNHANDLED_EXCEPTION = (108, 'unhandled exception')

class GsiBaseError(ExportError):
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

class ParametersError(GsiBaseError):
    _error_code = ErrorCode.PARAMETERS

class ArgumentsError(GsiBaseError):
    _error_code = ErrorCode.ARGUMENTS

class LoggingError(GsiBaseError):
    _error_code = ErrorCode.LOGGING

class ExportActionError(GsiBaseError):
    _error_code = ErrorCode.EXPORT_ACTION

class ResponseError(GsiBaseError):
    _error_code = ErrorCode(ErrorCode.RESPONSE)

class DBConnectionError(GsiBaseError):
    _error_code = ErrorCode.DB_CONNECTION

class DBCommandError(GsiBaseError):
    _error_code = ErrorCode.DB_COMMAND

class UnhandledExceptionError(GsiBaseError):
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

class AttributeListAction(ListAction):
    def __call__(self, parser, namespace, values, option_string=None):
        try:
            strtobool(values)
            self.dest = values
        except ValueError:
            # not a boolean; should be a list
            super.__call__(parser, namespace, values, option_string)

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

                parser_args = { 'usage' : '%(prog)s series=<DRMS series list> db_host=<db host> [ --log-file=<log file path> ] [ --logging-level=<critical/error/warning/info/debug> ] [ -N/--dbname=<db name> ] [ -P/--dbport=<db port> ] [ -U/--dbuser=<db user>] [ -w/--webserver=<host> ]' }
                if program_name is not None and len(program_name) > 0:
                    parser_args['prog'] = program_name

                parser = CmdlParser(**parser_args)

                # required
                parser.add_argument('series', help='a comma-separated list of series to be checked', metavar='<DRMS series>', action=ListAction, dest='series', required=True)
                parser.add_argument('dbhost', help='the machine hosting the database that contains export requests from this site', metavar='<db host>', dest='db_host', required=True)

                # optional
                parser.add_argument('--keywords', help='a list of DRMS keywords for which information is to be displayed', action=AttributeListAction, dest='keywords', default=True)
                parser.add_argument('--links', help='a list of DRMS links for which information is to be displayed', action=AttributeListAction, dest='links', default=True)
                parser.add_argument('-l', '--log-file', help='the path to the log file', metavar='<log file>', dest='log_file', default=log_file)
                parser.add_argument('-L', '--logging-level', help='the amount of logging to perform; in order of increasing verbosity: critical, error, warning, info, debug', metavar='<logging level>', dest='logging_level', action=DrmsLogLevelAction, default=DrmsLogLevel.ERROR)
                parser.add_argument('-N', '--dbname', help='the name of the database that contains export requests', metavar='<db name>', dest='db_name', default=db_name)
                parser.add_argument('-P', '--dbport', help='the port on the host machine that is accepting connections for the database', metavar='<db host port>', dest='db_port', type=int, default=db_port)
                parser.add_argument('-r', '--parse-record-sets', help='if set, `series` is a list of record-set specifications', action='store_true', dest='parse_record_sets')
                parser.add_argument('--segments', help='a list of DRMS segments for which information is to be displayed', action=AttributeListAction, dest='segments', default=True)
                parser.add_argument('-U', '--dbuser', help='the name of the database user account', metavar='<db user>', dest='db_user', default=db_user)
                parser.add_argument('-w', '--webserver', help='the webserver invoking this script', metavar='<webserver>', action=create_webserver_action(drms_params), dest='webserver', default=name_to_ws_obj(None, drms_params))

                arguments = Arguments(parser=parser, args=args)
            else:
                # `program_args` has all `arguments` values, in final form; validate them
                # `series` is a str
                def extract_module_args(*, series, db_host, parse_record_sets=False, log=None, db_name=db_name, db_port=db_port, db_user=db_user, keywords=True, links=True, segments=True, webserver=None):
                    arguments = {}

                    arguments['series'] = series # a list
                    arguments['db_host'] = db_host
                    arguments['parse_record_sets'] = parse_record_sets
                    arguments['db_name'] = db_name
                    arguments['db_port'] = db_port
                    arguments['db_user'] = db_user
                    arguments['keywords'] = False if (type(keywords) == list and len(keywords) == 0) else keywords
                    arguments['links'] = False if (type(links) == list and len(links) == 0) else links
                    arguments['segments'] = False if (type(segments) == list and len(segments) == 0) else segments
                    arguments['webserver'] = name_to_ws_obj(webserver, drms_params) # sets webserver.public = True if webserver is None

                    GetSeriesInfoAction.set_log(log)

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

def get_response(series_info, log):
    log.write_debug([ f'[ get_response ]' ])

    response_dict = deepcopy(series_info._d)
    info_status = response_dict['status']

    # move status so it does not conflict with response status value
    del response_dict['status']
    response_dict['info_status'] = info_status

    if info_status == 0:
        response_dict['status_code'] = StatusCode.SUCCESS
    elif info_status == 1:
        response_dict['status_code'] = StatusCode.FAILURE
    else:
        # there should be no other status possible
        response_dict['status_code'] = StatusCode.FAILURE

    return Response.generate_response(**response_dict)

def get_series_info_from_db(*, series, cursor, keywords=None, links=None, segments=None, log, series_info):
    # series attributes
    # only one series
    if series.lower() not in series_info:
        series_info[series.lower()] = {}

        # extract namespace - instead of having the server do this, which is what should happen, do it here since the
        # format is simple; this avoids the overhead of using securedrms ssh and then executing a process remotely
        try:
            name_space = series[:series.index('.')]
        except ValueError as exc:
            raise ArgumentsError(exc_info=sys_exc_info(), error_message=f'invalid series `{series}`: {str(exc)}')

        # PG stores all timestamp with timezone column values in UTC
        sql = f"SELECT seriesname AS series, author, owner, unitsize, archive, retention, tapegroup, primary_idx as drmsprimekey, dbidx, to_timestamp(created || ' UTC', 'YYYY-MM-DD HH24:MI:SS') AS created, description FROM {name_space}.drms_series WHERE lower(seriesname) = '{series.lower()}'"

        try:
            log.write_debug([ f'[ get_series_info_from_db ] executing sql: {sql}' ])
            cursor.execute(sql)
            rows = cursor.fetchall()
            if len(rows) == 0:
                raise DBCommandError(error_message=f'unknown DRMS data series `{series}`')

            # there is only one row
            row = rows[0]

            keys = ('series', 'author', 'owner', 'unitsize', 'archive', 'retention', 'tapegroup', 'drmsprimekey', 'dbindex', 'created', 'description')
            vals = list(row[0:7])
            vals.append([ key.strip() for key in row[7].split(',') ])
            vals.append([ key.strip() for key in row[8].split(',') ])
            vals.append(row[9].strftime('%Y-%m-%d %H:%M:%S UTC'))
            vals.append(row[10])

            series_info[series.lower()]['attributes'] = dict(zip(keys, vals))
        except PGError as exc:
            raise DBCommandError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

        if keywords is not None:
            # only one series
            sql = ''
            series_info[series.lower()]['keywords'] = {}

            if type(keywords) == bool:
                if keywords:
                    # all keywords
                    sql = f"SELECT series, keyword, datatype, keyworddefault AS constantvalue, unit, isconstant, (flags >> 16)::integer AS rank, description FROM drms_followkeywordlink('{series.lower()}', NULL)"
            else:
                # keywords specified in list
                sql_elements = []
                for keyword in keywords:
                    sql_elements.append(f"SELECT series, keyword, datatype, keyworddefault AS constantvalue, unit, isconstant, (flags >> 16)::integer AS rank, description FROM drms_followkeywordlink('{series.lower()}', '{keyword.lower()}')")

                sql = '\nUNION\n'.join(sql_elements)

            if len(sql) > 0:
                try:
                    log.write_debug([ f'[ get_series_info_from_db ] executing sql: {sql}' ])
                    cursor.execute(sql)
                    keywords_with_info = set()
                    rows = cursor.fetchall()

                    for row in rows:
                        if row[1].lower() not in keywords_with_info:
                            series_info[series.lower()]['keywords'][row[1].lower()] = { 'data-type' : row[2], 'constant-value': row[5] if row[5] else 'na', 'physical-unit' : row[4], 'rank' : row[6], 'description' : row[7] }

                            keywords_with_info.add(row[1].lower())

                    # add elements for keywords in list but not in series
                    try:
                        iterator = iter(keywords)
                        for keyword in iterator:
                            if keyword.lower() not in keywords_with_info:
                                series_info[series.lower()]['keywords'][keyword.lower()] = { 'data-type' : 'na', 'constant-value': 'na', 'physical-unit' : 'na', 'rank' : -1, 'description' : 'unknown keyword' }
                    except:
                        pass
                except PGError as exc:
                    raise DBCommandError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

        if links is not None:
            # only one series
            sql = ''
            series_info[series.lower()]['links'] = {}

            if type(links) == bool:
                if links:
                    # all links
                    sql = f"SELECT seriesname AS series, linkname AS link, target_seriesname AS tail_series, type, description FROM {name_space}.drms_link WHERE lower(seriesname) = '{series.lower()}'"
            else:
                # links specified in list
                links_list = ','.join([ f"'{link.lower()}'" for link in links ])
                sql = f"SELECT seriesname AS series, linkname AS link, target_seriesname AS tail_series, type, description FROM {name_space}.drms_link WHERE lower(seriesname) = '{series.lower()}' AND lower(linkname) IN ({links_list})"

            if len(sql) > 0:
                try:
                    log.write_debug([ f'[ get_series_info_from_db ] executing sql: {sql}' ])
                    cursor.execute(sql)
                    links_with_info = set()
                    rows = cursor.fetchall()

                    for row in rows:
                        if row[1].lower() not in links_with_info:
                            series_info[series.lower()]['links'][row[1].lower()] = { 'tail-series' : row[2], 'type': row[2], 'description' : row[3] }
                            links_with_info.add(row[1].lower())

                    # add elements for links in list but not in series
                    try:
                        iterator = iter(links)
                        for link in iterator:
                            if link.lower() not in links_with_info:
                                series_info[series.lower()]['links'][link.lower()] = { 'tail-series' : 'na', 'type': 'na', 'description' : 'unknown link' }
                    except:
                        pass
                except PGError as exc:
                    raise DBCommandError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

        if segments is not None:
            # only one series
            sql = ''
            series_info[series.lower()]['segments'] = {}

            if type(segments) == bool:
                if segments:
                    # all segments
                    sql = f"SELECT series, segment, datatype, segnum, scope, numaxes, dimensions, unit, protocol, description FROM drms_followsegmentlink('{series.lower()}', NULL)"
            else:
                # segments specified in list
                sql_elements = []
                for segment in segments:
                    sql_elements.append(f"SELECT series, segment, datatype, segnum, scope, numaxes, dimensions, unit, protocol, description FROM drms_followsegmentlink('{series.lower()}', '{segment.lower()}')")

                sql = '\nUNION\n'.join(sql_elements)

            if len(sql) > 0:
                try:
                    log.write_debug([ f'[ get_series_info_from_db ] executing sql: {sql}' ])
                    cursor.execute(sql)
                    segments_with_info = set()
                    rows = cursor.fetchall()

                    for row in rows:
                        if row[1].lower() not in segments_with_info:
                            series_info[series.lower()]['segments'][row[1].lower()] = { 'data-type' : row[2], 'segment-number': row[3], 'scope' : row[4], 'number-axes' : row[5], 'dimensions' : row[6], 'physical-unit' : row[7], 'protocol' : row[8],'description' : row[9] }

                            segments_with_info.add(row[1].lower())

                    # add elements for segments in list but not in series
                    try:
                        iterator = iter(segments)
                        for segment in iterator:
                            if segment.lower() not in segments_with_info:
                                series_info[series.lower()]['segments'][segment.lower()] = { 'data-type' : 'na', 'segment-number': -1, 'scope' : 'na', 'number-axes' : -1, 'dimensions' : 'na', 'physical-unit' : 'na', 'protocol' : 'na','description' : 'unknown segment' }
                    except:
                        pass
                except PGError as exc:
                    raise DBCommandError(exc_info=sys_exc_info(), error_message=f'{str(exc)}')

def perform_action(*, action_obj, is_program, program_name=None, **kwargs):
    # catch all expections so we can always generate a response
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

        # if parse_record_set is True, then arguments.series is a list of record-set specifications; must extract series
        if arguments.parse_record_sets:
            log.write_debug([ f'[ perform_action ] parsing record sets' ])
            # parse specifications
            series = []
            for one_specification in arguments.series:
                log.write_debug([ f'[ perform_action ] specification `{one_specification}`' ])
                action = Action.action(action_type='parse_specification', args={ 'specification' : one_specification, 'db_host' : arguments.db_host, 'webserver' : arguments.webserver.host, 'log' : log })
                response = action()

                log.write_debug([ f'[ perform_action ] parse action response `{str(response)}`' ])

                if not isinstance(response, ErrorResponse):
                    for subset in response.attributes.subsets:
                        log.write_debug([ f'[ perform_action ] series extracted: `{subset.seriesname}`' ])
                        series.append(subset.seriesname)
                else:
                    log.write_error([ f'[ perform_action ] {response.attributes.error_message}'])
        else:
            series = arguments.series

        # connect to DB directly (fastest method)
        log.write_debug([ f'[ perform_action ] direct DB connection' ])

        with pg_connect(database=arguments.db_name, host=arguments.db_host, port=str(arguments.db_port), user=arguments.db_user) as conn:
            with conn.cursor() as cursor:
                unique_series = list(OrderedDict.fromkeys(series))
                keywords = arguments.keywords if type(arguments.keywords) == bool else list(OrderedDict.fromkeys(arguments.keywords))
                links = arguments.links if type(arguments.links) == bool else list(OrderedDict.fromkeys(arguments.links))
                segments = arguments.segments if type(arguments.segments) == bool else list(OrderedDict.fromkeys(arguments.segments))

                series_info = {}

                for one_series in unique_series:
                    get_series_info_from_db(series=one_series, cursor=cursor, keywords=keywords, links=links, segments=segments, log=log, series_info=series_info)

        response_dict = series_info
        response_dict['status_code'] = StatusCode.SUCCESS

        response = Response.generate_response(**response_dict)
    except GsiBaseError as exc:
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
        log.write_info([ f'[ perform_action ] request complete; status {response.status_code.description()}' ])

    return response

# for use in export web app
from action import Action
from parse_specification import ParseSpecificationAction
class GetSeriesInfoAction(Action):
    actions = [ 'get_series_info' ]

    _log = None

    def __init__(self, *, method, series, db_host, parse_record_sets=False, log=None, db_name=None, db_port=None, db_user=None, keywords=None, links=None, segments=None, webserver=None):
        self._method = getattr(self, method)
        self._series = series # py list
        self._db_host = db_host # host webserver uses (private webserver uses private db host)
        self._webserver = webserver # dict
        self._options = {}
        self._options['parse_record_sets'] = parse_record_sets
        self._options['log'] = log
        self._options['db_name'] = db_name
        self._options['db_port'] = db_port
        self._options['db_user'] = db_user
        self._options['keywords'] = keywords # py list
        self._options['links'] = links # py list
        self._options['segments'] = segments # py list
        self._options['webserver'] = webserver # host name - gets converted to object in `get_arguments()`

    def get_series_info(self):
        response = perform_action(action_obj=self, is_program=False, series=self._series, db_host=self._db_host, options=self._options)
        return response

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

    # `series` is a py list of series OR a list of record-set specifications
    @classmethod
    def is_valid_series_set(cls, series, db_host, webserver, log):
        cls.set_log(log)
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
            action = Action.action(action_type='parse_specification', args={ 'specification' : ','.join(series), 'db_host' : db_host_resolved, 'webserver' : webserver, 'log' : cls._log })
            response = action()

            # record-set specifications are allowed in lieu of series names, so it is OK if filters exist
            is_valid = False if isinstance(response, ErrorResponse) else True
        except:
            is_valid = False

        return is_valid


if __name__ == "__main__":
    response = perform_action(action_obj=None, is_program=True)
    print(response.generate_json())

    # Always return 0. If there was an error, an error code (the 'status' property) and message (the 'statusMsg' property) goes in the returned HTML.
    sys_exit(0)
else:
    pass
