from datetime import datetime
from enum import Enum
from ftplib import FTP
from io import BytesIO
from json import loads as json_loads, decoder
import re
from subprocess import run
from sys import stderr as sys_stderr

# connect
SWPC_FTP_HOST = 'ftp.swpc.noaa.gov'
SWPC_FTP_DIRECTORY = 'pub/forecasts/SRS'
SWPC_FILE_NAME_PATTERN = r'^\s*\d\d\d\dSRS.txt\s*$'

JSOC_INFO_BIN = '/home/jsoc/cvs/Development/JSOC/bin/linux_avx/jsoc_info'
SET_INFO_BIN = '/home/jsoc/cvs/Development/JSOC/bin/linux_avx/set_info'
TIME_CONVERT_BIN = '/home/jsoc/cvs/Development/JSOC/bin/linux_avx/time_convert'
ACTIVE_REGION_SERIES = 'jsoc.noaa_active_regions'
ACTIVE_REGION_SERIES_FILE = 'swpc_file'
ACTIVE_REGION_SERIES_FILE_MOD = 'swpc_file_mod_time'

class ParseState(Enum):
    # info failure error codes
    START = (0, '')
    HEADER_1 = (1, 'header 1')
    HEADER_2 = (2, 'header 2')
    PART_1 = (2, 'part 1')
    PART_1A = (3, 'part 1a')
    PART_2 = (3, 'part 2')
    END = (4, 'end')


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

# discard other parts (part II)

header_1_regex = re.compile(HEADER_1_PATTERN, re.IGNORECASE)
header_2_regex = re.compile(HEADER_2_PATTERN, re.IGNORECASE)
part_1_start_regex = re.compile(PART_I_START_PATTERN)
part_i_regex = re.compile(PART_I_PATTERN, re.IGNORECASE)
part_ia_start_regex = re.compile(PART_IA_START_PATTERN, re.IGNORECASE)
part_ia_regex = re.compile(PART_IA_PATTERN, re.IGNORECASE)
part_2_start_regex = re.compile(PART_II_START_PATTERN, re.IGNORECASE)

def uniquify_id(id, observation_time):
    if observation_time > datetime.strptime('20000101', '%Y%m%d') and int(id) < 5000:
        unique_id = str(int(id) + 10000)
    else:
        unique_id = id
    return unique_id

def process_content(*, state, line_content, state_data):
    return_state = state

    if state == ParseState.START:
        match_obj = header_1_regex.search(line_content)
        if match_obj is not None:
            matches = match_obj.groups()
            if len(matches) != 2:
                raise

            report_month, report_year = matches
            state_data['report_month'] = datetime.strptime(f'{report_month}{report_year}', '%b%Y').strftime('%m')
            state_data['report_year'] = report_year

            return_state = ParseState.HEADER_1
    elif state == ParseState.HEADER_1:
        match_obj = header_2_regex.search(line_content)
        if match_obj is not None:
            matches = match_obj.groups()
            if len(matches) != 1:
                raise

            observation_month, = matches
            observation_month = datetime.strptime(f'{observation_month}', '%b').strftime('%m')
            state_data['observation_month'] = observation_month

            return_state = ParseState.HEADER_2
    elif state == ParseState.HEADER_2:
        match_obj = part_1_start_regex.search(line_content)
        if match_obj is not None:
            matches = match_obj.groups()
            if match_obj is not None:
                if len(matches) != 3:
                    raise

            observation_day, observation_hour, observation_minute = matches

            if int(state_data['report_month']) == 1 and int(state_data['observation_month']) == 12:
                observation_year = str(int(state_data['report_year']) - 1)
            else:
                observation_year = state_data['report_year']

            swpc_observation_time_str = f'{observation_year}.{state_data["observation_month"]}.{observation_day}_{observation_hour}:{observation_minute}_UT'
            state_data['swpc_observation_time_str'] = swpc_observation_time_str

            return_state = ParseState.PART_1
    elif state == ParseState.PART_1:
        match_obj = part_i_regex.search(line_content)
        if match_obj is not None:
            print('XX1')
            matches = match_obj.groups()
            if len(matches) != 11:
                raise

            part_i_id = matches[0]
            part_i_location = []
            part_i_location.append(f'-{matches[2]}' if matches[1].lower() == 'e' else f'{matches[2]}')
            part_i_location.append(f'-{matches[4]}' if matches[3].lower() == 's' else f'{matches[4]}')
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
            keyword_dict['longitudnalextent'] = part_i_extent
            keyword_dict['zurichclass'] = part_i_zurich_classification
            keyword_dict['magnetictype'] = part_i_magnetic_classification
            keyword_dict['swpc_file'] = state_data['swpc_file']
            keyword_dict['swpc_file_mod_time'] = state_data['swpc_file_mod_time'].strftime('%Y%m%d_%H%M%S')

            print(f'setting keyword dict')
            state_data['keyword_dict'] = keyword_dict
        else:
            print('XX2')
            state_data['keyword_dict'] = None
            match_obj = part_ia_start_regex.search(line_content)
            if match_obj is not None:
                return_state = ParseState.PART_1A
    elif state == ParseState.PART_1A:
        match_obj = part_ia_regex.search(line_content)
        if match_obj is not None:
            matches = match_obj.groups()
            if len(matches) != 6:
                raise

            part_ia_id = matches[0]
            part_ia_location = []
            part_ia_location.append(f'-{matches[2]}' if matches[1].lower() == 'e' else f'{matches[2]}')
            part_ia_location.append(f'-{matches[4]}' if matches[3].lower() == 's' else f'{matches[4]}')
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
            keyword_dict['longitudnalextent'] = part_ia_extent
            keyword_dict['zurichclass'] = part_ia_zurich_classification
            keyword_dict['magnetictype'] = part_ia_magnetic_classification
            keyword_dict['swpc_file'] = state_data['swpc_file']
            keyword_dict['swpc_file_mod_time'] = state_data['swpc_file_mod_time'].strftime('%Y%m%d_%H%M%S')

            print(f'setting keyword dict')
            state_data['keyword_dict'] = keyword_dict
        else:
            state_data['keyword_dict'] = None
            match_obj = part_2_start_regex.search(line_content)
            if match_obj is not None:
                return_state = ParseState.PART_2
    elif state == ParseState.PART_2:
        return_state = ParseState.END
    elif state ==  ParseState.END:
        return_state = ParseState.END
    else:
        # bad state
        raise ValueError(f'bad state {str(state)}')

    return return_state

if __name__ == "__main__":
    with FTP(SWPC_FTP_HOST) as ftp_client:
        ftp_client.login()
        ftp_client.cwd(SWPC_FTP_DIRECTORY)
        files = ftp_client.mlsd()

        mod_times = None

        # check all files
        # midnight_today = datetime.strptime(datetime.today().strftime('%Y%m%d'), '%Y%m%d')
        command = [ JSOC_INFO_BIN, f'ds={ACTIVE_REGION_SERIES}[][]', f'op=rs_list', f's=1', f'key={ACTIVE_REGION_SERIES_FILE},{ACTIVE_REGION_SERIES_FILE_MOD}', f'DRMS_DBUTF8CLIENTENCODING=1', f'JSOC_DBHOST=hmidb2' ]

        try:
            print(f'running {" ".join(command)}')
            completed_proc = run(command, encoding='utf8', capture_output=True)
            if completed_proc.returncode != 0:
                print(f'return code: {str(completed_proc.returncode)}', file=sys.stderr)
                raise

            response = json_loads(completed_proc.stdout)
            if response.get('keywords', None) is not None:
                mod_times = dict(zip(response['keywords'][0]['values'], response['keywords'][1]['values']))
        except decoder.JSONDecodeError as exc:
            print(f'{str(exc)}', file=sys_stderr)
        except Exception as exc:
            print(f'{str(exc)}', file=sys_stderr)

        file_name_regex = re.compile(SWPC_FILE_NAME_PATTERN)

        files_to_process = {}
        for file in files:
            file_name = file[0]

            if file_name_regex.search(file_name) == None:
                continue

            file_mod_time = datetime.strptime(file[1]["modify"], "%Y%m%d%H%M%S") #YYYMMDDHHMMSS
            #if file_mod_time >= midnight_today:
            # exclude files that have already been ingested (using file name and timestamp)
            last_mod_time = mod_times.get(file_name, None) if mod_times is not None else None
            if last_mod_time is None or file_mod_time > last_mod_time:
                files_to_process[file_name] = file_mod_time

        for file_name, file_mod_time in files_to_process.items():
            state = ParseState.START
            with BytesIO() as stream:
                state_data = { "swpc_file" : file_name, "swpc_file_mod_time" : file_mod_time }
                print(f'downloading {file_name}')
                ftp_client.retrbinary(f'RETR {file_name}', stream.write)

                stream.seek(0)
                while True:
                    line_content = stream.readline().decode().strip()
                    print(f'processing line {line_content}')

                    if line_content != b'':
                        print(f'state is {str(state)}')
                        next_state = process_content(state=state, line_content=line_content, state_data=state_data)
                        print(f'next state is {str(next_state)}')

                        if state == ParseState.HEADER_2:
                            # call time_convert to get proper DRMS time string
                            command = [ TIME_CONVERT_BIN, f'time={state_data["swpc_observation_time_str"]}', f'o=cal']

                            try:
                                print(f'running {" ".join(command)}')
                                completed_proc = run(command, encoding='utf8', capture_output=True)
                                if completed_proc.returncode != 0:
                                    print(f'return code: {str(completed_proc.returncode)}', file=sys.stderr)
                                    raise
                                state_data['observation_time'] = datetime.strptime(completed_proc.stdout.strip(), '%Y.%m.%d_%H:%M:%S_%Z')
                            except Exception as exc:
                                print(f'{str(exc)}', file=sys_stderr)
                        elif state == ParseState.PART_1 or state == ParseState.PART_1A:
                            # call set_info with key=val pairs
                            keyword_dict = state_data.get('keyword_dict', None)

                            if keyword_dict is not None:
                                command = [ SET_INFO_BIN, '-c', f'ds=ACTIVE_REGION_SERIES' ]
                                command.extend([ f'{element[0]}={str(element[1])}' for element in keyword_dict.items() ])
                                print(f'command type is {type(command)}')
                                print(f'running {" ".join(command)}')
                        elif state == ParseState.END:
                            break

                        state = next_state
                    else:
                        break

                break

        # algorithm for obs time
        # - get 1-based month index and year from SRS Number XXX line (the report month and year)
        # - if no valid month index and year, return failure
        # - get month abbreviation from Report compiled line, unless line is '::::::::::', in which
        #   case get the month from the line after '::::::::::' (the observation month)
        # - if no valid month, return failure
        # - convert month to 1-based index
        # - if report month index is 1 and observation month is 12, then subtract one from report year
        # - get day from I.  Regions with Sunspots. line
