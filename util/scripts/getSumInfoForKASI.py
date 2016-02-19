#!/usr/bin/env python

# This program takes as input a file containing columns of data. The last column must be a list of SUNUMs.
# For each of the SUNUMs presented, this program retrieves the path to the referenced SU (if it is exists and is online).
# The original, non-SUNUMs columns are printed without modification, followed by a column of SUMS paths.

from __future__ import print_function
import sys
import re
import os
import pwd
import psycopg2
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../../include'))
from drmsparams import DRMSParams


class SumsDrmsParams(DRMSParams):
    def __init__(self):
        super(SumsDrmsParams, self).__init__()

    def get(self, name):
        val = super(SumsDrmsParams, self).get(name)

        if val is None:
            raise Exception('drmsParams', 'Unknown DRMS parameter: ' + name + '.')
        return val
        
def getSumMainInfo(sunums, bytesMap, seriesMap, cursor):
    cmd = 'SELECT ds_index, bytes, owning_series FROM sum_main WHERE ds_index IN (' + ','.join(sunums) + ')'
    cursor.execute(cmd)
    rows = cursor.fetchall()
    for row in rows:
        sunum = row[0]
        sunumStr = str(sunum)
        bytesMap[sunumStr] = str(row[1])
        seriesMap[sunumStr] = row[2]

if __name__ == "__main__":
    sumsDrmsParams = SumsDrmsParams()
    if sumsDrmsParams is None:
        raise Exception('drmsParams', 'Unable to locate DRMS parameters file (drmsparams.py).')
            
    regexp = re.compile(r'^\s*(\d+)\s*$')
    file = sys.argv[1]
    with open(file, 'r') as fin:
        with psycopg2.connect(database=sumsDrmsParams.get('DBNAME') + '_sums', user='production', host=sumsDrmsParams.get('SUMS_DB_HOST'), port=sumsDrmsParams.get('SUMPGPORT')) as conn:
            with conn.cursor() as cursor:    
                sunums = []
                allSunums = []
                bytesMap = {}
                seriesMap = {}
                iloop = 1024
                                
                for line in fin:
                    line = line.rstrip()
                    match = regexp.match(line)
                    if match is not None:
                        sunumStr = match.group(1)
                        if int(sunumStr) >= 0:
                            sunums.append(sunumStr)
                            allSunums.append(sunumStr)

                            if iloop == 0:
                                getSumMainInfo(sunums, bytesMap, seriesMap, cursor)
                                sunums = []
                                iloop = 1024
                            else:
                                iloop -= 1
                
                if len(sunums) > 0:
                    getSumMainInfo(sunums, bytesMap, seriesMap, cursor)
                    sunums = []
                    
                # Match up the series information with the paths obtained.
                for sunumStr in allSunums:
                    if sunumStr in bytesMap and sunumStr in seriesMap:
                        bytesStr = bytesMap[sunumStr]
                        series = seriesMap[sunumStr]
                    else:
                        bytesStr = 'badSUNUM'
                        series = 'badSUNUM'
                        
                    print(sunumStr + '\t' + bytesStr + '\t' + series)
