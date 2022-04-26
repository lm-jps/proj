from os.path import join as path_join
from subprocess import PIPE, Popen, run

NEW_DATA_FILE = 'new-data.txt'


DEFAULT_BIN_PATH = '/home/jsoc/cvs/Development/JSOC/bin/linux_avx'
SET_INFO_BIN = 'set_info'
SHOW_INFO_BIN = 'show_info'

with open(NEW_DATA_FILE) as new_data_f:
    for line in new_data_f:
        data_row = line.strip()
        if len(data_row) > 0:
            time_string, noaa_ar_id = data_row.split()
            command = [ path_join(DEFAULT_BIN_PATH, SHOW_INFO_BIN), '-ak', 'DRMS_DBUTF8CLIENTENCODING=1', f'su_rsb.noaa_activeregions[{time_string}][{noaa_ar_id}]' ]
            input_proc = Popen(command, stdout=PIPE)

            command = [ path_join(DEFAULT_BIN_PATH, SET_INFO_BIN), '-k', 'DRMS_DBUTF8CLIENTENCODING=1', 'ds=jsoc.noaa_active_regions', 'JSOC_DBHOST=hmidb2' ]
            ingest_proc = run(command, check=True, stdin=input_proc.stdout)
            # print(f'show_info stdout {input_proc.stdout.read().decode()}')
            input_proc.wait()


