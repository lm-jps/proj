#!/usr/bin/env python3

import argparse
import collections
import copy
import distutils.util
import enum
import json
import os.path
import psycopg2
import sys

from drms_export import Connection, Error as ExportError, ErrorCode as ExportErrorCode, ErrorResponse, ExpServerBaseError, get_arguments as ss_get_arguments, get_message, Response, send_message
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
    EXPORT_SERVER = (109, 'export-server communication error')

class RecordScope(enum.Enum):
    VARIABLE = (0, 'variable')
    CONSTANT = (1, 'constant')
    INDEX = (100, 'index')
    TS_EQ = (1000, 'ts_eq')
    SLOT = (1001, 'slot')
    ENUM = (1002, 'enum')
    CARR = (1003, 'carr')
    TS_SLOT = (1004, 'ts_slot')

    def __new__(cls, int_value, str_value):
        member = object.__new__(cls)
        member._int_value = int_value
        member._str_value = str_value

        if not hasattr(cls, '_all_members'):
            cls._all_members = {}
        cls._all_members[str(int_value)] = member

        return member

    def __int__(self):
        return self._int_value

    def __str__(self):
        return self._str_value

    @classmethod
    def _missing_(cls, int_value):
        return cls._all_members[str(int_value)]

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

class ExportServerError(GsiBaseError):
    _error_code = ErrorCode.EXPORT_SERVER

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
    cls = type('WebserverAction', (argparse.Action,),
    {
        '_drms_params' : drms_params,
        '__call__' : webserver_action,
    })

    return cls

class AttributeListAction(ListAction):
    def __call__(self, parser, namespace, values, option_string=None):
        try:
            distutils.util.strtobool(values)
            self.dest = values
        except ValueError:
            # not a boolean; should be a list
            super().__call__(parser, namespace, values, option_string)

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
                raise ParametersError(exc_info=sys.exc_info(), error_message=str(exc))

            if is_program:
                try:
                    log_file = os.path.join(drms_params.get_required('EXPORT_LOG_DIR'), DEFAULT_LOG_FILE)
                except DPMissingParameterError as exc:
                    raise ParametersError(exc_info=sys.exc_info(), error_message=str(exc))

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

    response_dict = copy.deepcopy(series_info._d)
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
            raise ArgumentsError(exc_info=sys.exc_info(), error_message=f'invalid series `{series}`: {str(exc)}')

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

            # we really should first check to see if a keyword is an index keyword, but that would require another DB query to get
            # keyword information; but we know that if a prime-key keyword ends in '_index' that it is an index keyword that has
            # an external name that is the original name with the '_index' lopped off
            vals.append([ key.strip()[0:key.strip().rfind('_index')] if (key.strip().rfind('_index') + len('_index') == len(key.strip())) else key.strip() for key in row[7].split(',') ])
            vals.append([ key.strip()[0:key.strip().rfind('_index')] if (key.strip().rfind('_index') + len('_index') == len(key.strip())) else key.strip() for key in row[8].split(',') ])
            vals.append(row[9].strftime('%Y-%m-%d %H:%M:%S UTC'))
            vals.append(row[10])

            series_info[series.lower()]['attributes'] = dict(zip(keys, vals))
        except psycopg2.Error as exc:
            raise DBCommandError(exc_info=sys.exc_info(), error_message=f'{str(exc)}')

        if keywords is not None:
            # only one series
            sql = ''
            series_info[series.lower()]['keywords'] = {}

            if type(keywords) == bool:
                if keywords:
                    # all keywords
                    sql = f"SELECT series, keyword, link, datatype, keyworddefault AS constantvalue, unit, scope_type, (flags >> 16)::integer AS rank, description FROM drms_followkeywordlink2('{series.lower()}', NULL) ORDER BY rank ASC"

                    sql2 = f"SELECT keyword, linkedkeyword FROM drms_linkedkeyword('{series.lower()}', NULL)"
            else:
                # keywords specified in list
                sql_elements = []
                sql_elements2 = []
                for keyword in keywords:
                    sql_elements.append(f"SELECT series, keyword, link, datatype, keyworddefault AS constantvalue, unit, scope_type, (flags >> 16)::integer AS rank, description FROM drms_followkeywordlink2('{series.lower()}', '{keyword.lower()}')")

                    sql_elements2.append(f"SELECT keyword, linkedkeyword FROM drms_linkedkeyword('{series.lower()}', '{keyword.lower()}')")

                sql = '\nUNION\n'.join(sql_elements)
                sql2 = '\nUNION\n'.join(sql_elements2) + 'ORDER BY rank ASC\n'

            if len(sql) > 0 and len(sql2) > 0:
                try:
                    log.write_debug([ f'[ get_series_info_from_db ] executing sql: {sql}' ])
                    cursor.execute(sql)
                    keywords_with_info = set()
                    rows = cursor.fetchall()

                    for row in rows:
                        keyword = row[1]
                        link = row[2]
                        data_type = row[3]
                        default_value = row[4]
                        unit = row[5]
                        scope = RecordScope(int(row[6]))
                        rank = row[7]
                        description = row[8]

                        if keyword.lower() not in keywords_with_info:
                            series_info[series.lower()]['keywords'][keyword.lower()] = { 'name' : keyword, 'linked-keyword' : f'{str(link)}->' if link is not None and len(str(link)) > 0 else None, 'data-type' : data_type, 'scope' : str(scope), 'default-value' : default_value, 'constant-value': default_value if scope == RecordScope.CONSTANT else None, 'physical-unit' : unit, 'is-slotted' : (int(scope) >= int(RecordScope.TS_EQ)), 'rank' : rank, 'description' : description }

                            keywords_with_info.add(keyword.lower())

                    log.write_debug([ f'[ get_series_info_from_db ] executing sql: {sql2}' ])
                    cursor.execute(sql2)
                    keywords_with_info2 = set()
                    rows = cursor.fetchall()

                    for row in rows:
                        keyword = row[0]
                        linked_keyword = row[1]

                        if keyword.lower() not in keywords_with_info2 and keyword.lower() in keywords_with_info:
                            linked_keyword_str = series_info[series.lower()]['keywords'][keyword.lower()].get('linked-keyword', None)
                            if linked_keyword_str is not None:
                                series_info[series.lower()]['keywords'][keyword.lower()]['linked-keyword'] = f'{linked_keyword_str}{linked_keyword}'

                            keywords_with_info2.add(keyword.lower())

                    # add elements for keywords in list but not in series
                    try:
                        iterator = iter(keywords)
                        for keyword in iterator:
                            if keyword.lower() not in keywords_with_info:
                                series_info[series.lower()]['keywords'][keyword.lower()] = { 'name' : keyword, 'linked-keyword' : None, 'data-type' : None, 'scope' : None, 'default-value' : None, 'constant-value': None, 'physical-unit' : None, 'rank' : None, 'description' : None }
                    except:
                        pass
                except psycopg2.Error as exc:
                    raise DBCommandError(exc_info=sys.exc_info(), error_message=f'{str(exc)}')

        if links is not None:
            # only one series
            sql = ''
            series_info[series.lower()]['links'] = {}

            if type(links) == bool:
                if links:
                    # all links
                    sql = f"SELECT seriesname AS series, linkname AS link, target_seriesname AS child_series, type, description FROM {name_space}.drms_link WHERE lower(seriesname) = '{series.lower()}'"
            else:
                # links specified in list
                links_list = ','.join([ f"'{link.lower()}'" for link in links ])
                sql = f"SELECT seriesname AS series, linkname AS link, target_seriesname AS child_series, type, description FROM {name_space}.drms_link WHERE lower(seriesname) = '{series.lower()}' AND lower(linkname) IN ({links_list})"

            if len(sql) > 0:
                try:
                    log.write_debug([ f'[ get_series_info_from_db ] executing sql: {sql}' ])
                    cursor.execute(sql)
                    links_with_info = set()
                    rows = cursor.fetchall()

                    for row in rows:
                        link = row[1]
                        child_series = row[2]
                        link_type = row[3]
                        description = row[4]

                        if link.lower() not in links_with_info:
                            series_info[series.lower()]['links'][link.lower()] = { 'name' : link, 'child-series' : child_series, 'type': link_type, 'description' : description }
                            links_with_info.add(link.lower())

                    # add elements for links in list but not in series
                    try:
                        iterator = iter(links)
                        for link in iterator:
                            if link.lower() not in links_with_info:
                                series_info[series.lower()]['links'][link.lower()] = { 'name' : link, 'child-series' : None, 'type': None, 'description' : None }
                    except:
                        pass
                except psycopg2.Error as exc:
                    raise DBCommandError(exc_info=sys.exc_info(), error_message=f'{str(exc)}')

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
                        segment = row[1]
                        data_type = row[2]
                        segment_number = row[3]
                        scope = row[4]
                        number_axes = row[5]
                        dimensions = row[6]
                        unit = row[7]
                        protocol = row[8]
                        description = row[9]

                        if segment.lower() not in segments_with_info:
                            series_info[series.lower()]['segments'][segment.lower()] = { 'name' : segment, 'data-type' : data_type, 'segment-number': segment_number, 'scope' : scope, 'number-axes' : number_axes, 'dimensions' : dimensions, 'physical-unit' : unit, 'protocol' : protocol,'description' : description }

                            segments_with_info.add(segment.lower())

                    # add elements for segments in list but not in series
                    try:
                        iterator = iter(segments)
                        for segment in iterator:
                            if segment.lower() not in segments_with_info:
                                series_info[series.lower()]['segments'][segment.lower()] = { 'name' : segment, 'data-type' : None, 'segment-number': None, 'scope' : None, 'number-axes' : None, 'dimensions' : None, 'physical-unit' : None, 'protocol' : None,'description' : None }
                    except:
                        pass
                except psycopg2.Error as exc:
                    raise DBCommandError(exc_info=sys.exc_info(), error_message=f'{str(exc)}')

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
            raise ArgumentsError(exc_info=sys.exc_info(), error_message=f'{str(exc)}')
        except Exception as exc:
            raise ArgumentsError(exc_info=sys.exc_info(), error_message=f'{str(exc)}')

        if action_obj is None or action_obj.log is None:
            try:
                formatter = DrmsLogFormatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
                log = DrmsLog(arguments.log_file, arguments.logging_level, formatter)
                if action_obj is not None:
                    action_obj.log = log
            except Exception as exc:
                raise LoggingError(exc_info=sys.exc_info(), error_message=f'{str(exc)}')
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

        with psycopg2.pg_connect(database=arguments.db_name, host=arguments.db_host, port=str(arguments.db_port), user=arguments.db_user) as conn:
            with conn.cursor() as cursor:
                unique_series = list(collections.OrderedDict.fromkeys(series))
                keywords = arguments.keywords if type(arguments.keywords) == bool else list(collections.OrderedDict.fromkeys(arguments.keywords))
                links = arguments.links if type(arguments.links) == bool else list(collections.OrderedDict.fromkeys(arguments.links))
                segments = arguments.segments if type(arguments.segments) == bool else list(collections.OrderedDict.fromkeys(arguments.segments))

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
        response = UnhandledExceptionError(exc_info=sys.exc_info(), error_message=f'{str(exc)}').response
        error_message = str(exc)

        if log:
            log.write_error([ error_message ])
        elif is_program:
            print(error_message)

    if log is not None:
        log.write_info([ f'[ perform_action ] request complete; status {response.attributes.drms_export_status_description}' ])

    return response

def send_request(request, connection, log):
    json_message = json.dumps(request)
    send_message(connection, json_message)
    message = get_message(connection)

    return message

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
                    raise ParametersError(exc_info=sys.exc_info(), error_message=str(exc))
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

class ShowSeriesLegacyResponse(Response):
    def remove_extraneous(self, response_dict, keys):
        for key in keys:
            if key in response_dict:
                del response_dict[key]

        return response_dict

    # override parent to remove attributes not present in legacy API
    def generate_serializable_dict(self):
        serializable_dict = copy.deepcopy(super().generate_serializable_dict())
        sanitized_serializable_dict = self.remove_extraneous(serializable_dict, [ 'drms_export_status', 'drms_export_status_code', 'drms_export_status_description' ] )

        return sanitized_serializable_dict

class ShowSeriesLegacyAction(Action):
    actions = [ 'legacy_get_series_list' ]

    _log = None

    def __init__(self, *, method, series_regex, db_host=None, log=None, **kwargs):
        self._method = getattr(self, method)
        self._series_regex = series_regex # ds
        self._options = {}
        self._options['db_host'] = db_host
        self._options['log'] = log if log is not None else DrmsLog(None, None, None)
        self._legacy_arguments = self.convert_booleans(kwargs)

    def convert_booleans(self, args_dict):
        converted = {}
        for key, val in args_dict.items():
            if type(val) == bool:
                converted[key] = int(val)
            else:
                converted[key] = val

        return converted

    def perform_action(self):
        try:
            log = self._options['log']
            nested_arguments = self.get_nested_arguments()

            try:
                response = self.send_request(nested_arguments=nested_arguments)
            except ExpServerBaseError as exc:
                raise ExportServerError(exc_info=sys.exc_info(), error_message=f'{str(exc)}')
            except Exception as exc:
                raise ExportServerError(exc_info=sys.exc_info(), error_message=f'{str(exc)}')
        except GsiBaseError as exc:
            response = exc.response
            error_message = exc.message

            if log:
                log.write_error([ error_message ])
        except Exception as exc:
            response = UnhandledExceptionError(exc_info=sys.exc_info(), error_message=f'{str(exc)}').response
            error_message = str(exc)

            if log:
                log.write_error([ error_message ])

        return response

    def get_nested_arguments(self):
        log = self._options['log']

        if self._options['db_host'] is None:
            drms_params = DRMSParams()

            if drms_params is None:
                raise ParametersError(error_message='unable to locate DRMS parameters package')

            self._options['db_host'] = drms_params.get_required('SERVER')

        # use socket server to call jsoc_fetch
        return ss_get_arguments(is_program=False, module_args={})

    def get_request_status_code(self, response_dict):
        error_code, status_code = self.get_request_status(response_dict)
        code = error_code if error_code is not None else status_code

        return code

    def get_response(self, client_response_dict, log):
        log.write_debug([ f'[ {self.__class__.__name__}.get_response ]' ])

        response_dict = copy.deepcopy(client_response_dict)
        response_dict['status_code'] = self.get_request_status_code(response_dict)
        response = ShowSeriesLegacyResponse.generate_response(**response_dict)

        return response

    def legacy_get_series_list(self):
        return self.perform_action()

    def send_request(self, *, nested_arguments):
        log = self._options['log']

        with Connection(server=nested_arguments.server, listen_port=nested_arguments.listen_port, timeout=nested_arguments.message_timeout, log=log) as connection:
            message = { 'request_type' : 'legacy_show_series', 'series_regex' : self._series_regex, 'db_host' : self._options['db_host'] }
            message.update(self._legacy_arguments) # error if the user has provided non-expected arguments
            response = send_request(message, connection, log)

            # message is raw JSON from checkAddress.py
            legacy_show_series_dict = json.loads(response)

            if legacy_show_series_dict.get('export_server_status') == 'export_server_error':
                raise ExportServerError(error_message=f'{legacy_show_series_dict["error_message"]}')

            message = { 'request_type' : 'quit' }
            send_request(message, connection, log)

        response = self.get_response(legacy_show_series_dict, log)
        return response

    def get_request_status(self, response_dict):
        error_code = None
        status_code = None

        try:
            error_code = ErrorCode(int(response_dict['status']))
        except KeyError:
            pass

        if error_code is None:
            try:
                status_code = StatusCode(int(response_dict['status']))
            except KeyError:
                raise ExportServerError(exc_info=sys.exc_info(), error_message=f'unexpected show_series status returned {str(response_dict["status"])}')

        return (error_code, status_code)

class ShowSeriesPublicLegacyAction(ShowSeriesLegacyAction):
    actions = [ 'legacy_get_public_series_list' ]

    _log = None

    def __init__(self, *, method, series_regex, db_host=None, log=None, **kwargs):
        self._method = getattr(self, method)
        self._series_regex = series_regex
        self._options = {}
        self._options['db_host'] = db_host
        self._options['log'] = log if log is not None else DrmsLog(None, None, None)
        self._legacy_arguments = self.convert_booleans(kwargs)

    def legacy_get_public_series_list(self):
        return self.perform_action()

    def send_request(self, *, nested_arguments):
        log = self._options['log']

        with Connection(server=nested_arguments.server, listen_port=nested_arguments.listen_port, timeout=nested_arguments.message_timeout, log=log) as connection:
            message = { 'request_type' : 'legacy_show_ext_series', 'series_regex' : self._series_regex, 'db_host' : self._options['db_host'] }
            message.update(self._legacy_arguments) # error if the user has provided non-expected arguments
            response = send_request(message, connection, log)

            # message is raw JSON from checkAddress.py
            legacy_show_series_dict = json.loads(response)

            if legacy_show_series_dict.get('export_server_status') == 'export_server_error':
                raise ExportServerError(error_message=f'{legacy_show_series_dict["error_message"]}')

            message = { 'request_type' : 'quit' }
            send_request(message, connection, log)

        response = self.get_response(legacy_show_series_dict, log)
        return response

    def get_request_status(self, response_dict):
        error_code = None
        status_code = None

        error_str = response_dict.get('errMsg')
        error_val = int(error_str is not None and len(error_str.strip()) > 0)

        try:
            error_code = ErrorCode(error_val)
        except KeyError:
            pass

        if error_code is None:
            try:
                status_code = StatusCode(error_val)
            except KeyError:
                raise ExportServerError(exc_info=sys.exc_info(), error_message=f'unexpected show_series status returned {error_str}')

        return (error_code, status_code)

def run_tests():
    formatter = DrmsLogFormatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    log = DrmsLog(sys.stdout, DrmsLogLevelAction.string_to_level('debug'), formatter)

    log.write_debug([ 'testing `legacy_get_series_list`'])
    sspla = ShowSeriesLegacyAction(method='legacy_get_series_list', series_regex='^su_arta', db_host='hmidb', log=log)
    response = sspla()
    response_dict = response.generate_serializable_dict()
    log.write_debug([ str(response_dict), '' ])

    log.write_debug([ 'testing `legacy_get_public_series_list`'])
    sspla = ShowSeriesPublicLegacyAction(method='legacy_get_public_series_list', series_regex='^hmi[.]m', info=1, db_host='hmidb2', log=log)
    response = sspla()
    response_dict = response.generate_serializable_dict()
    log.write_debug([ str(response_dict), '' ])

if __name__ == "__main__":
    response = perform_action(action_obj=None, is_program=True)
    print(response.generate_json())

    # Always return 0. If there was an error, an error code (the 'status' property) and message (the 'statusMsg' property) goes in the returned HTML.
    sys.exit(0)
else:
    pass
