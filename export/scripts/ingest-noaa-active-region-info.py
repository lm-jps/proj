from datetime import datetime
from enum import Enum
from ftplib import FTP
from io import BytesIO
from json import loads as json_loads, decoder
from os.path import join as path_join
from pytz import timezone
import re
from subprocess import run
from sys import exit as sys_exit, stderr as sys_stderr

from drms_parameters import DRMSParams, DPMissingParameterError
from drms_utils import Arguments as Args, ArgumentsError as ArgsError, CmdlParser, Formatter as DrmsLogFormatter, Log as DrmsLog, LogLevel as DrmsLogLevel, LogLevelAction as DrmsLogLevelAction

DEFAULT_LOG_FILE = 'inar_log.txt'

SWPC_FTP_HOST = 'ftp.swpc.noaa.gov'
SWPC_FTP_DIRECTORY = 'pub/forecasts/SRS'
SWPC_FILE_NAME_PATTERN = r'^\s*\d\d\d\dSRS.txt\s*$'
DRMS_TIME_STRING_PATTERN = r'^\s*(\d+)[.](\d+)[.](\d+)(_((\d+)(:(\d+)(:((\d+)([.]\d+)?))?)?)?)?(_([a-z]+))?\s*$'

JSOC_INFO_BIN = '/home/jsoc/cvs/Development/JSOC/bin/linux_avx/jsoc_info'
SET_INFO_BIN = '/home/jsoc/cvs/Development/JSOC/bin/linux_avx/set_info'
TIME_CONVERT_BIN = '/home/jsoc/cvs/Development/JSOC/bin/linux_avx/time_convert'
ACTIVE_REGION_SERIES = 'jsoc.noaa_active_regions'
ACTIVE_REGION_SERIES_FILE = 'swpc_file'
ACTIVE_REGION_SERIES_FILE_MOD = 'swpc_file_mod_time'

PART_IA_DEFAULT_ZURICH_CLASSIFICATION = '?'
PART_IA_DEFAULT_MAGNETIC_CLASSIFICATION = '?'
PART_IA_DEFAULT_NUMBER_SPOTS = 0
PART_IA_DEFAULT_AREA = 0
PART_IA_DEFAULT_EXTENT = -1

HEADER_1_PATTERN = r'^\s*srs number\s+[0-9a-z ]+?\d+\s+([a-z]+)\s+(\d\d\d\d)\s*$'
HEADER_2_PATTERN = r'^\s*report compiled\s+[a-z ]+?\d+\s+([a-z]+)\s*$'
PART_I_START_PATTERN = r'^\s*I[.]\s+.+?(\d+)[/](\d\d)(\d\d)[zZ]\s*$'
PART_I_PATTERN = r'^\s*(\d+)\s+([ns])(\d+)([ew])(\d+)\s+(\d+)\s+(\d+)\s+([a-z]+)\s+(\d+)\s+(\d+)\s+([a-z]+)\s*$'
PART_IA_START_PATTERN = r'^\s*IA[.]\s+.+?\s*$'
PART_IA_PATTERN = r'\s*(\d+)\s+([ns])(\d+)([ew])(\d+)\s+(\d+)\s*$'
PART_II_START_PATTERN = r'^\s*II[.]\s+.+?\s*$'

class ErrorCode(Enum):
    NONE = (0, 'no error - success!')
    ARGUMENTS = (1, 'bad arguments')
    PARAMETERS = (2, 'failure obtaining DRMS parameters')
    LOGGING = (3, 'failure creating or writing to log')
    HEADER_FORMAT = (4, 'header content does not conform to expected format')
    PART_I_TITLE_FORMAT = (5, 'part I title does not conform to expected format')
    PART_I_DATA_FORMAT = (6, 'part I data line does not conform to expected format')
    PART_IA_DATA_FORMAT = (7, 'part IA data line does not conform to expected format')
    PROCESSING_STATE = (8, 'unknown processing state')
    CHILD_PROCESS = (9, 'failure running child process')
    TIME_STRING_FORMAT = (10, 'invalid DRMS time string')
    REGEX_API = (11, 'failure creating regular-expression API')

    def __new__(cls, code, description):
        member = object.__new__(cls)
        member._code = code
        member._description = description
        return member

    def __int__(self):
        return self._code

    @property
    def description(self, **kwargs):
        return self._description.format(**kwargs)

class ParseState(Enum):
    # info failure error codes
    START = (0, '')
    HEADER_1 = (1, 'header 1')
    HEADER_2 = (2, 'header 2')
    PART_1 = (2, 'part 1')
    PART_1A = (3, 'part 1a')
    PART_2 = (3, 'part 2')
    END = (4, 'end')

class NoaaBaseException(Exception):
    def __init__(self, message):
        self._message = message

    def __str__(self):
        return self._message

    @property
    def error_code(self):
        return self.__class__._error_code

class ArgumentsError(NoaaBaseException):
    _error_code = ErrorCode.ARGUMENTS

class ParametersError(NoaaBaseException):
    _error_code = ErrorCode.PARAMETERS

class HeaderFormatError(NoaaBaseException):
    _error_code = ErrorCode.HEADER_FORMAT

class PartITitleFormatError(NoaaBaseException):
    _error_code = ErrorCode.PART_I_TITLE_FORMAT

class PartIDataFormatError(NoaaBaseException):
    _error_code = ErrorCode.PART_I_DATA_FORMAT

class PartIADataFormatError(NoaaBaseException):
    _error_code = ErrorCode.PART_IA_DATA_FORMAT

class ProcessingStateError(NoaaBaseException):
    _error_code = ErrorCode.PROCESSING_STATE

class ChildProcessError(NoaaBaseException):
    _error_code = ErrorCode.CHILD_PROCESS

class DrmsTimeStringFormatError(NoaaBaseException):
    _error_code = ErrorCode.TIME_STRING_FORMAT

class RegExApiError(NoaaBaseException):
    _error_code = ErrorCode.REGEX_API

class RegExType(type):
    _regex_classes = {}
    _required_attributes = [ 'regex_pattern', 'regex' ]

    def __new__(cls, name, bases, attributes):
        for required_att in cls._required_attributes:
            if required_att not in attributes:
                raise RegExApiError(f'missing required attribute {required_att} in RegExType class {str(cls)}')

        if name in cls._regex_classes:
            raise RegExApiError(f'duplicate RegExType class {name}')

        return super(RegExType, cls).__new__(cls, name, bases, attributes)

    def __init__(cls, name, bases, attributes):
        RegExType._regex_classes[name] = cls

    @classmethod
    def get_regex_classes(cls):
        regex_classes = []
        for name, regex_class in cls._regex_classes.items():
            regex_classes.append(regex_class)

        return regex_classes

    @classmethod
    def get_regex_names(cls):
        regex_names = []
        for name, regex_class in cls._regex_classes.items():
            regex_names.append(regex_class.regex_name)

        return regex_names

def regex_search(self, text):
    return self.regex.search(text)

def regex_str(self):
    return self.__name__

def create_regexes(regex_patterns):
    for regex_name, regex_pattern, ignore_case in regex_patterns:
        if ignore_case:
            regex = re.compile(regex_pattern, re.IGNORECASE)
        else:
            regex = re.compile(regex_pattern)

        regex_cls = RegExType(f'{regex_name}Meta', (object,),
        {
            'regex_pattern' : regex_pattern,
            'regex' : regex,
            'search' : regex_search,
            'regex_name' : regex_name,
            '__str__' : regex_str,
        })

        globals()[regex_name] = regex_cls()

REGEX_PATTERNS = (
    ('SwpcFileNameRegex', SWPC_FILE_NAME_PATTERN, False),
    ('DrmsTimeStringRegex', DRMS_TIME_STRING_PATTERN, True),
    ('Header1Regex', HEADER_1_PATTERN, True),
    ('Header2Regex', HEADER_2_PATTERN, True),
    ('PartIStartRegex', PART_I_START_PATTERN, False),
    ('PartIRegex', PART_I_PATTERN, True),
    ('PartIAStartRegex', PART_IA_START_PATTERN, True),
    ('PartIARegex', PART_IA_PATTERN, True),
    ('PartIIStartRegex', PART_II_START_PATTERN, True),
)

class Arguments(Args):
    _arguments = None

    @classmethod
    def get_arguments(cls, drms_params):
        try:
            # db_user = drms_params.get_required('WEB_DBUSER')
            log_file = path_join(drms_params.get_required('DRMS_LOG_DIR'), DEFAULT_LOG_FILE)
        except DPMissingParameterError as exc:
            raise ParametersError(exc_info=sys_exc_info(), error_message=str(exc))

        parser_args = { 'usage' : '%(prog)s  [ -e/--end-day=<end observation day> ] [ -s/--start-day=<start observation day> ] [ -L/--logging-level=<critical/error/warning/info/debug> ]' }

        parser = CmdlParser(**parser_args)

        # optional
        parser.add_argument('-e', '--end-day', help='end observation day (%m%d)', metavar='<end day>', dest='end_day', default=None)
        parser.add_argument('-l', '--log-file', help='the path to the log file', metavar='<log file>', dest='log_file', default=log_file)
        parser.add_argument('-L', '--logging-level', help='the amount of logging to perform; in order of increasing verbosity: critical, error, warning, info, debug', metavar='<logging level>', dest='logging_level', action=DrmsLogLevelAction, default=DrmsLogLevel.ERROR)
        parser.add_argument('-s', '--start-day', help='start observation day (%m%d)', metavar='<start day>', dest='start_day', default=None)

        cls._arguments = Arguments(parser=parser, args=None)
        return cls._arguments

def uniquify_id(id, observation_time):
    dt = datetime(2000, 1, 1, tzinfo=observation_time.tzinfo)
    if observation_time > dt and int(id) < 5000:
        unique_id = str(int(id) + 10000)
    else:
        unique_id = id
    return unique_id

def check_time_interval(state_data):
    start_day = state_data['time_window'].get('start_day', None) # no year (1900)
    end_day = state_data['time_window'].get('end_day', None) # no year (1900)
    process_data = True

    # remove year
    observation_month_and_day = datetime.strptime(state_data['observation_time'].strftime('%m%d'), '%m%d')

    if (start_day is not None and observation_month_and_day < start_day) or (end_day is not None and observation_month_and_day > end_day):
        if start_day is not None and end_day is not None:
            time_interval = f't >= {start_day.strftime("%b %d")} and t <= {end_day.strftime("%b %d")}'
        elif start_day is not None:
            time_interval = f't >= {start_day.strftime("%b %d")}'
        elif end_day is not None:
            time_interval = f't <= {end_day.strftime("%b %d")}'

        log.write_info([ f'[ check_time_interval ] observation time {state_data["observation_time"].strftime("%b %d")} outside of requested time interval {time_interval}'])
        process_data = False

    return process_data

def set_time_window_endpoint(time_window, endpoint, endpoint_day):
    dt = None
    try:
        dt = datetime.strptime(f'{endpoint_day}', '%b%d')
    except:
        try:
            dt = datetime.strptime(f'{endpoint_day}', '%m%d')
        except:
            try:
                dt = datetime.strptime(f'{endpoint_day}', '%B%d')
            except:
                log.write_info([ f'[ set_time_window_endpoint ] WARNING: invalid time window endpoint {endpoint_day}; ignoring' ])

    if dt is not None:
        time_window[endpoint] = dt
        log.write_debug([ f'[ set_time_window_endpoint ] window endpoint time {time_window[endpoint].strftime("%b %d")}' ])

def process_content(*, state, line_content, state_data):
    return_state = state

    if state == ParseState.START:
        match_obj = Header1Regex.search(line_content)
        if match_obj is not None:
            matches = match_obj.groups()
            if len(matches) != 2:
                raise HeaderFormatError(f'unexpected header format: {line_content}')

            report_month, report_year = matches
            state_data['report_month'] = datetime.strptime(f'{report_month}{report_year}', '%b%Y').strftime('%m')
            state_data['report_year'] = report_year

            return_state = ParseState.HEADER_1
    elif state == ParseState.HEADER_1:
        match_obj = Header2Regex.search(line_content)
        if match_obj is not None:
            matches = match_obj.groups()
            if len(matches) != 1:
                raise HeaderFormatError(f'unexpected header format: {line_content}')

            observation_month, = matches
            observation_month = datetime.strptime(f'{observation_month}', '%b').strftime('%m')
            state_data['observation_month'] = observation_month

            return_state = ParseState.HEADER_2
    elif state == ParseState.HEADER_2:
        match_obj = PartIStartRegex.search(line_content)
        if match_obj is not None:
            matches = match_obj.groups()
            if match_obj is not None:
                if len(matches) != 3:
                    raise PartITitleFormatError(f'unexpected part-I title format: {line_content}')

            observation_day, observation_hour, observation_minute = matches

            if int(state_data['report_month']) == 1 and int(state_data['observation_month']) == 12:
                observation_year = str(int(state_data['report_year']) - 1)
            else:
                observation_year = state_data['report_year']

            swpc_observation_time_str = f'{observation_year}.{state_data["observation_month"]}.{observation_day}_{observation_hour}:{observation_minute}_UT'
            state_data['swpc_observation_time_str'] = swpc_observation_time_str
            log.write_debug([ f'[ process_content ] swpc observation time string {swpc_observation_time_str}' ])

            return_state = ParseState.PART_1
    elif state == ParseState.PART_1:
        process_data = check_time_interval(state_data)

        if process_data:
            match_obj = PartIRegex.search(line_content)
            if match_obj is not None:
                matches = match_obj.groups()
                if len(matches) != 11:
                    raise PartIDataFormatError(f'unexpected part-I data format: {line_content}')

                part_i_id = matches[0]
                part_i_location = []
                part_i_location.append(f'-{matches[2]}' if matches[1].lower() == 's' else f'{matches[2]}')
                part_i_location.append(f'-{matches[4]}' if matches[3].lower() == 'e' else f'{matches[4]}')
                part_i_carrington_lon, part_i_area = matches[5:7]
                part_i_zurich_classification = f'{matches[7][0].upper()}{matches[7][1:].replace("-", "").lower()}'
                part_i_extent, part_i_number_spots, = matches[8:10]
                part_i_magnetic_classification = f'{matches[10][0].upper()}{matches[10][1:].replace("-", "").lower()}'

                # call set_info with key=val pairs
                keyword_dict = {}
                keyword_dict['regionnumber'] = uniquify_id(part_i_id, state_data['observation_time'])
                keyword_dict['observationtime'] = state_data['observation_time'].strftime('%Y.%m.%d_%H:%M:%S_%Z')
                keyword_dict['area'] = part_i_area
                keyword_dict['latitudehg'] = part_i_location[0]
                keyword_dict['longitudehg'] = part_i_carrington_lon
                keyword_dict['longitudecm'] = part_i_location[1]
                keyword_dict['longitudinalextent'] = part_i_extent
                keyword_dict['zurichclass'] = part_i_zurich_classification
                keyword_dict['spotcount'] = part_i_number_spots
                keyword_dict['magnetictype'] = part_i_magnetic_classification
                keyword_dict['swpc_file'] = state_data['swpc_file']
                keyword_dict['swpc_file_mod_time'] = state_data['swpc_file_mod_time'].strftime('%Y%m%d_%H%M%S')

                state_data['keyword_dict'] = keyword_dict
                log.write_debug([ f'[ process_content ] keyword dict {str(keyword_dict)}' ])
            else:
                state_data['keyword_dict'] = None
                match_obj = PartIAStartRegex.search(line_content)
                if match_obj is not None:
                    return_state = ParseState.PART_1A
        else:
            return_state = ParseState.END
    elif state == ParseState.PART_1A:
        process_data = check_time_interval(state_data)

        if process_data:
            match_obj = PartIARegex.search(line_content)
            if match_obj is not None:
                matches = match_obj.groups()
                if len(matches) != 6:
                    raise PartIADataFormatError(f'unexpected part-IA data format: {line_content}')

                part_ia_id = matches[0]
                part_ia_location = []
                part_ia_location.append(f'-{matches[2]}' if matches[1].lower() == 's' else f'{matches[2]}')
                part_ia_location.append(f'-{matches[4]}' if matches[3].lower() == 'e' else f'{matches[4]}')
                part_ia_carrington_lon = matches[5]
                part_ia_area = PART_IA_DEFAULT_AREA
                part_ia_zurich_classification = f'{PART_IA_DEFAULT_ZURICH_CLASSIFICATION[0].upper()}{PART_IA_DEFAULT_ZURICH_CLASSIFICATION[1:].replace("-", "").lower()}'
                part_ia_extent = PART_IA_DEFAULT_EXTENT
                part_ia_number_spots = PART_IA_DEFAULT_NUMBER_SPOTS
                part_ia_magnetic_classification = f'{PART_IA_DEFAULT_MAGNETIC_CLASSIFICATION[0].upper()}{PART_IA_DEFAULT_MAGNETIC_CLASSIFICATION[1:].replace("-", "").lower()}'

                keyword_dict = {}
                keyword_dict['regionnumber'] = uniquify_id(part_ia_id, state_data['observation_time'])
                keyword_dict['observationtime'] = state_data['observation_time'].strftime('%Y.%m.%d_%H:%M:%S_%Z')
                keyword_dict['area'] = part_ia_area
                keyword_dict['latitudehg'] = part_ia_location[0]
                keyword_dict['longitudehg'] = part_ia_carrington_lon
                keyword_dict['longitudecm'] = part_ia_location[1]
                keyword_dict['longitudinalextent'] = part_ia_extent
                keyword_dict['zurichclass'] = part_ia_zurich_classification
                keyword_dict['spotcount'] = part_ia_number_spots
                keyword_dict['magnetictype'] = part_ia_magnetic_classification
                keyword_dict['swpc_file'] = state_data['swpc_file']
                keyword_dict['swpc_file_mod_time'] = state_data['swpc_file_mod_time'].strftime('%Y%m%d_%H%M%S')

                state_data['keyword_dict'] = keyword_dict
                log.write_debug([ f'[ process_content ] keyword dict {str(keyword_dict)}' ])
            else:
                state_data['keyword_dict'] = None
                match_obj = PartIIStartRegex.search(line_content)
                if match_obj is not None:
                    return_state = ParseState.PART_2
        else:
            return_state = ParseState.END
    elif state == ParseState.PART_2:
        return_state = ParseState.END
    elif state ==  ParseState.END:
        return_state = ParseState.END
    else:
        # bad state
        raise ProcessingStateError(f'unexpected file-processing state {str(state)}')

    return return_state

def fetch_mod_times():
    log.write_debug([ f'[ fetch_mod_times ] fetching data-file modification times' ])
    mod_times = None

    # check all files
    command = [ JSOC_INFO_BIN, f'ds={ACTIVE_REGION_SERIES}[][]', f'op=rs_list', f's=1', f'key={ACTIVE_REGION_SERIES_FILE},{ACTIVE_REGION_SERIES_FILE_MOD}', f'DRMS_DBUTF8CLIENTENCODING=1', f'JSOC_DBHOST=hmidb2' ]

    try:
        log.write_debug([ f'[ fetch_mod_times ] running {" ".join(command)}' ])
        completed_proc = run(command, encoding='utf8', capture_output=True)
        if completed_proc.returncode != 0:
            raise ChildProcessError(f'non-zero return code: {str(completed_proc.returncode)}')

        response = json_loads(completed_proc.stdout)
        if response.get('keywords', None) is not None:
            mod_times = dict(zip(response['keywords'][0]['values'], [ datetime.strptime(time_string, '%Y%m%d_%H%M%S') for time_string in response['keywords'][1]['values'] ]))
    except decoder.JSONDecodeError as exc:
        raise ChildProcessError(f'child process did not return valid JSON response (command was `{" ".join(command)}``)')
    except Exception as exc:
        raise ChildProcessError(f'failure obtaining modification times of previously ingested data files (command was `{" ".join(command)}``)')

    log.write_info([ f'[ fetch_mod_times ] successfully fetched modification time for {str(len(mod_times)) if mod_times is not None else "0"} files' ])
    return mod_times

def get_new_files(files, mod_times):
    log.write_debug([ f'[ get_new_files ] filtering out up-to-date files' ])
    files_to_process = {}
    for file in files:
        file_name = file[0]

        if SwpcFileNameRegex.search(file_name) == None:
            log.write_info([ f'[ get_new_files ] WARNING: invalid file-name format `{file_name}`; ignoring' ])
            continue

        file_mod_time = datetime.strptime(file[1]["modify"], "%Y%m%d%H%M%S") #YYYMMDDHHMMSS

        # exclude files that have already been ingested (using file name and timestamp)
        last_mod_time = mod_times.get(file_name, None) if mod_times is not None else None
        if last_mod_time is None or file_mod_time > last_mod_time:
            log.write_debug([ f'[ get_new_files ] * file `{file_name}` (modification time {file_mod_time.strftime("%Y-%m-%dT%H:%M:%S")}) has NOT been ingested' ])
            files_to_process[file_name] = file_mod_time
        else:
            log.write_debug([ f'[ get_new_files ] file `{file_name}` (modification time {file_mod_time.strftime("%Y-%m-%dT%H:%M:%S")}) has ALREADY been ingested' ])

    log.write_info([ f'[ get_new_files ] found {str(len(files_to_process))} files to process' ])
    return files_to_process

def set_observation_time():
    # call time_convert to get proper DRMS time string
    log.write_debug([ f'[ set_observation_time ] converting SWPC observation time `{state_data["swpc_observation_time_str"]}` to a DRMS time string'])
    command = [ TIME_CONVERT_BIN, f'time={state_data["swpc_observation_time_str"]}', f'o=cal']

    try:
        log.write_debug([ f'[ set_observation_time ] running {" ".join(command)}'])
        completed_proc = run(command, encoding='utf8', capture_output=True)
    except Exception as exc:
        raise ChildProcessError(f'{str(exc)}')

    if completed_proc.returncode != 0:
        raise ChildProcessError(f'non-zero return code: {str(completed_proc.returncode)}')

    # python does not have a way to parse time zone from a non-standardized time string (like a DRMS time string);
    # use a regular expression
    drms_time_string = completed_proc.stdout.strip()
    log.write_debug([ f'[ set_observation_time ] DRMS time string `{drms_time_string}`'])
    match_obj = DrmsTimeStringRegex.search(drms_time_string)
    if match_obj is not None:
        matches = match_obj.groups()

        if len(matches) == 14:
            year_str, month_str, day_str = matches[0:3]
            hour_str = matches[5] if matches[5] is not None else '0'
            minute_str = matches[7] if matches[7] is not None else '0'
            second_str = matches[10] if matches[10] is not None else '0'
            microsecond = float(matches[11]) * 1000000 if matches[11] is not None else 0
            time_zone_str = matches[13] if matches[13] else None
        else:
            raise DrmsTimeStringFormatError(f'invalid DRMS time string `{drms_time_string}`')
    else:
        raise DrmsTimeStringFormatError(f'invalid DRMS time string `{drms_time_string}`')

    state_data['observation_time'] = datetime(int(year_str), int(month_str), int(day_str), hour=int(hour_str), minute=int(minute_str), second=int(second_str), microsecond=microsecond, tzinfo=timezone(time_zone_str) if time_zone_str is not None else None)

    log.write_info([ f'[ set_observation_time ] observation time `{state_data["observation_time"].strftime("%Y-%m-%dT%H:%M:%S%Z")}`' ])

if __name__ == "__main__":
    exit_code = None

    try:
        drms_params = DRMSParams()
        arguments = Arguments.get_arguments(drms_params)

        try:
            formatter = DrmsLogFormatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
            log = DrmsLog(arguments.log_file, arguments.logging_level, formatter)
        except Exception as exc:
            raise LoggingError(f'{str(exc)}')

        log.write_info([ f'[ __main__ ] arguments: {str(arguments)}' ])

        create_regexes(REGEX_PATTERNS)
        log.write_debug([f'successfully created regular expression metaclasses {str(RegExType.get_regex_names())}'])

        time_window = {}
        if arguments.start_day is not None:
            set_time_window_endpoint(time_window, 'start_day', arguments.start_day)

        if arguments.end_day is not None:
            set_time_window_endpoint(time_window, 'end_day', arguments.end_day)

        if 'start_day' in time_window and 'end_day' in time_window and time_window['start_day'] > time_window['end_day']:
            raise ArgumentsError(f'start day argument `{time_window["start_day"].strftime("%b %d")}` is AFTER end day argument `{time_window["end_day"].strftime("%b %d")}`')

        log.write_info([ f'[ __main__ ] requested time interval {time_window["start_day"].strftime("%b %d") if time_window["start_day"] is not None else "*"} <-> {time_window["end_day"].strftime("%b %d") if time_window["end_day"] is not None else "*"}' ])

        with FTP(SWPC_FTP_HOST) as ftp_client:
            log.write_info([ f'[ __main__ ] connected to {SWPC_FTP_HOST}' ])
            ftp_client.login()
            log.write_debug([ f'[ __main__ ] logged in' ])
            ftp_client.cwd(SWPC_FTP_DIRECTORY)
            log.write_debug([ f'[ __main__ ] cwd to {SWPC_FTP_DIRECTORY}' ])
            files = ftp_client.mlsd()

            mod_times = fetch_mod_times()
            files_to_process = get_new_files(files, mod_times)

            for file_name, file_mod_time in files_to_process.items():
                state = ParseState.START
                with BytesIO() as stream:
                    state_data = { "swpc_file" : file_name, "swpc_file_mod_time" : file_mod_time, "time_window" : time_window }

                    log.write_debug([ f'[ __main__ ] downloading `{file_name}`' ])
                    ftp_client.retrbinary(f'RETR {file_name}', stream.write)

                    stream.seek(0)
                    while True:
                        line_content = stream.readline().decode().strip()
                        log.write_debug([ f'[ __main__ ] processing line {line_content}' ])

                        if line_content != b'':
                            log.write_debug([ f'[ __main__ ] state is {str(state)}' ])
                            next_state = process_content(state=state, line_content=line_content, state_data=state_data)

                            if state == ParseState.HEADER_2:
                                set_observation_time()
                            elif state == ParseState.PART_1 or state == ParseState.PART_1A:
                                # call set_info with key=val pairs
                                keyword_dict = state_data.get('keyword_dict', None)

                                if keyword_dict is not None:
                                    command = [ SET_INFO_BIN, '-c', f'ds={ACTIVE_REGION_SERIES}', f'DRMS_DBUTF8CLIENTENCODING=1', f'JSOC_DBHOST=hmidb2' ]
                                    command.extend([ f'{element[0]}={str(element[1])}' for element in keyword_dict.items() ])

                                    try:
                                        log.write_info([ f'running {" ".join(command)}' ])
                                        completed_proc = run(command, encoding='utf8', capture_output=True)
                                        if completed_proc.returncode != 0:
                                            raise ChildProcessError(f'non-zero return code: {str(completed_proc.returncode)}')
                                    except Exception as exc:
                                        raise ChildProcessError(f'failure ingesting active-region data for ({state_data["observation_time"].strftime("%Y-%m-%dT%H:%M:%S%Z")}, {str(keyword_dict["regionnumber"])})')

                            elif state == ParseState.END:
                                break

                            state = next_state
                        else:
                            break

        exit_code = ErrorCode.NONE
        log.write_info([ f'[ __main__] {exit_code.description}' ])
    except NoaaBaseException as exc:
        exit_code = exc.error_code
        log.write_error([ f'[ __main__] failure : "{exit_code.description}", : error message "{str(exc)}"' ])

    sys_exit(int(exit_code))

        # algorithm for obs time
        # - get 1-based month index and year from SRS Number XXX line (the report month and year)
        # - if no valid month index and year, return failure
        # - get month abbreviation from Report compiled line, unless line is '::::::::::', in which
        #   case get the month from the line after '::::::::::' (the observation month)
        # - if no valid month, return failure
        # - convert month to 1-based index
        # - if report month index is 1 and observation month is 12, then subtract one from report year
        # - get day from I.  Regions with Sunspots. line
