import calendar
from datetime import datetime, timedelta
from json import load as json_load, loads as json_loads
import os
import pandas as pd
from subprocess import PIPE, Popen, run
from sys import stderr as sys_stderr

# date when no missing filtergrams
reference_date = datetime.strptime('20180101', '%Y%m%d')

# 125K records
CHUNK_SIZE = 131072

monthly_count = {}
todays_month = datetime.today().strftime('%b')
todays_year = datetime.today().strftime('%Y')

def get_show_info_command(*, series, segment_list, year, month_index, time_chunk):
    spec_series = None
    spec_todays_year = None
    spec_month_index = None

    if spec_series is None:
        spec_series = series

    if spec_todays_year is None:
        spec_todays_year = year

    if spec_month_index is None:
        spec_month_index = month_index

    start_day = None
    end_day = None
    time_delta = time_chunk.days
    days_in_month = calendar.monthrange(int(spec_todays_year), spec_month_index)[1]

    last_iteration = False
    while not last_iteration:
        start_day = 1 if start_day is None else start_day + time_delta
        end_day = start_day + time_delta if start_day + time_delta < days_in_month else days_in_month

        if end_day == days_in_month:
            end_month_index  = spec_month_index + 1
            end_day = 1
            last_iteration = True
        else:
            end_month_index  = spec_month_index

        cmd = [ f'/home/jsoc/cvs/Development/JSOC/bin/linux_avx/show_info', f'ds={spec_series}[][? t_obs >= $({spec_todays_year}.{str(spec_month_index)}.{str(start_day)}) AND t_obs < $({spec_todays_year}.{str(end_month_index)}.{str(end_day)}) ?]', f'-oPx', f'key=quality', f'seg={segment_list}' ]

        yield cmd

# iterate over series
with open('filtergram_series.json') as series_file_obj:
    jsoc_dict = json_load(series_file_obj)
    all_series_info = jsoc_dict['series']
    data_frame = None

    for series_info in all_series_info:
        series = series_info['name']
        segments = series_info['segments']
        start_month = series_info.get('start_month', None)
        end_month = series_info.get('end_month', None)

        segment_list = ','.join(segments)

        # how many records per day for this series
        cmd = [ f'/home/jsoc/cvs/Development/JSOC/bin/linux_avx/jsoc_info', f'op=rs_summary', f'ds={series.lower()}[][? t_obs >= $({reference_date.strftime("%Y.%m.%d")}) AND t_obs < $({(reference_date + timedelta(days=1)).strftime("%Y.%m.%d")}) ?]', f's=1' ]

        print(f'running {" ".join(cmd)}')

        try:
            response = run(cmd, stdout=PIPE)
            json_dict = json_loads(response.stdout.decode())
            records_per_day = json_dict.get('count', None)
            if records_per_day is None:
                raise ValueError
        except:
            # next series
            print(f'cannot obtain data for reference date for series {series}', file=sys_stderr)
            continue

        for month_index, month in enumerate(calendar.month_abbr):
            if len(month) > 0:
                # compare against start and end
                month_dt = datetime.strptime(month, '%b')
                if start_month is not None:
                    if month_dt < datetime.strptime(start_month, '%b'):
                        continue

                if end_month is not None:
                    if month_dt > datetime.strptime(end_month, '%b'):
                        break

                if month == todays_month:
                    break

                print(f'tallying {month}', flush=True)

                if month not in monthly_count:
                    monthly_count[month] = 0

                hours_per_chunk = CHUNK_SIZE * 24 / records_per_day
                get_chunk_generator = get_show_info_command(series=series, segment_list=segment_list, year=todays_year, month_index=month_index, time_chunk=timedelta(hours=hours_per_chunk))

                # run show_info (which can handle many more records than jsoc_info)
                for cmd in get_chunk_generator:
                    print(f'running {" ".join(cmd)}', flush=True)
                    with Popen(cmd, stdout=PIPE) as proc:
                        try:
                            data_frame = pd.read_csv(proc.stdout, sep='\t')
                        except pd.errors.EmptyDataError:
                            print(f'no data for this show_info command', flush=True)
                            continue

                        for segment in segments:
                            # count all online files
                            for image in data_frame[data_frame[segment].str.match('^/')][segment]:
                                if os.path.exists(image):
                                    monthly_count[month] += 1

                            # neither DRMS nor SUMS keeps track of offline segment files; the best we can
                            # do is assume a file exists on tape if:
                            #   - the online status is 'N'
                            #   - the archive status is 'Y' or 'Pending'
                            #   - the top bit of the quality integer is 0
                            for quality in data_frame[(data_frame.online.str.lower() == 'n') & ((data_frame.archive.str.lower() == 'y') | (data_frame.archive.str.lower() == 'pending'))].quality:
                                if int(quality, 0) & 2147483648 == 0:
                                    monthly_count[month] += 1

                    data_frame = None
print(f'done')
print(f'\n{str(monthly_count)}', flush=True)
