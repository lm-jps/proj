#!/usr/bin/env python

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


def getPathsFromDb(sunums, sunumMap, cursor):
    cmd = 'SELECT ds_index, wd FROM sum_partn_alloc WHERE ds_index IN (' + ','.join(sunums) + ')'
    cursor.execute(cmd)
    rows = cursor.fetchall()
    for row in rows:
        sunum = row[0]                                    
        sunumStr = str(sunum)        
        sunumMap[sunumStr] = row[1]


if __name__ == "__main__":
    sumsDrmsParams = SumsDrmsParams()
    if sumsDrmsParams is None:
        raise Exception('drmsParams', 'Unable to locate DRMS parameters file (drmsparams.py).')
            
    regexp = re.compile(r'\s*(\S+)\s+(\S+)\s+(\S+)')
    file = sys.argv[1]
    with open(file, 'r') as fin:
        with psycopg2.connect(database=sumsDrmsParams.get('DBNAME') + '_sums', user='production', host=sumsDrmsParams.get('SUMS_DB_HOST'), port=sumsDrmsParams.get('SUMPGPORT')) as conn:
            with conn.cursor() as cursor:    
                sunums = []
                recInfo = []
                sunumMap = {}
                iloop = 32
                                
                for line in fin:
                    match = regexp.match(line)
                    if match is not None:
                        sunumStr = match.group(3)
                        if int(sunumStr) >= 0:
                            recInfo.append((match.group(1), match.group(2), sunumStr))
                            sunums.append(sunumStr)

                            if iloop == 0:
                                getPathsFromDb(sunums, sunumMap, cursor)
                                sunums = []
                                iloop = 32
                            else:
                                iloop -= 1
                
                if len(sunums) > 0:
                    getPathsFromDb(sunums, sunumMap, cursor)
                    sunums = []
                    
                # Match up the series information with the paths obtained.
                for rec in recInfo:
                    series = rec[0]
                    recnum = rec[1]
                    sunumStr = rec[2]
                    if sunumStr in sunumMap:
                        path = sunumMap[sunumStr]
                    else:
                        path = 'badSUNUM'
                        
                    print(series + '\t' + recnum + '\t' + sunumStr + '\t' + path)
