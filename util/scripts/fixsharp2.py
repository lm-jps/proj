#!/home/jsoc/bin/linux_x86_64/activepython

import sys
import os.path
import getopt
import pwd
import re
import psycopg2

# Return codes
RET_SUCCESS = 0
RET_INVALIDARG = 1
RET_DBCONNECT = 2
RET_SQL = 3

def GetArgs(args):
    istat = bool(0)
    optD = {}
    
    try:
        opts, remainder = getopt.getopt(args, "hs:t:k:u:n:h:p:d", ["source=", "targets=", "keys=", "unique=", "dbname=", "dbhost=", "dbport="])
    except getopt.GetoptError:
        print('Usage:\n  fixsharp2.py [-h] -s <source series> -t <target-series list> -k <keys> -u <prime key key-list> -n <db name> -h <db host> -p <db port> [-d]', file=sys.stderr)
        istat = bool(1)
    
    if istat == bool(0):
        for opt, arg in opts:
            if opt == '-h':
                print('Usage:\n  fixsharp2.py [-h] -s <source series> -t <target-series list> -k <keys> -u <prime key key-list> -n <db name> -h <db host> -p <db port> [-d]')
                sys.exit(0)
            elif opt in ("-s", "--source"):
                regexp = re.compile(r"\s*(\S+)\.(\S+)\s*")
                matchobj = regexp.match(arg)
                if matchobj is None:
                    istat = bool(1)
                else:
                    optD['ns'] = matchobj.group(1)
                    optD['table'] = matchobj.group(2)
            elif opt in ("-t", "--targets"):
                # Is the argument a file?
                if os.path.isfile(arg):
                    # If the argument is a file, parse it.
                    optD['targets'] = list()
                    
                    try:
                        with open(arg, 'r') as fin:
                            while True:
                                targetsRaw = fin.readlines(8192)
                                if not targetsRaw:
                                    break
                                targets = [target.strip(' \t\n,') for target in targetsRaw]
                                optD['targets'].extend(targets)
                    except IOError as exc:
                        type, value, traceback = sys.exc_info()
                        print(exc.strerror, file=sys.stderr)
                        print('Unable to open ' + "'" + value.filename + "'.", file=sys.stderr)
                        istat = bool(1)
                else:
                    # Otherwise, parse the argument itself.
                    optD['targets'] = arg.split(',') # a list
            elif opt in ("-k", "--keys"):
                # Is the argument a file?
                if os.path.isfile(arg):
                    # If the argument is a file, parse it.
                    optD['keys'] = list()
                    
                    try:
                        with open(arg, 'r') as fin:
                            while True:
                                keysRaw = fin.readlines(8192)
                                if not keysRaw:
                                    break
                                keys = [key.strip(' \t\n,') for key in keysRaw]
                                optD['keys'].extend(keys)
                    except IOError as exc:
                        type, value, traceback = sys.exc_info()
                        print(exc.strerror, file=sys.stderr)
                        print('Unable to open ' + "'" + value.filename + "'.", file=sys.stderr)
                        istat = bool(1)
                else:
                    # Otherwise, parse the argument itself.
                    optD['keys'] = arg.split(',') # a comma-separated list
            elif opt in ("-u", "--unique"):
                optD['unique'] = arg.split(',') # a comma-separated list
            elif opt in ("-n", "--dbname"):
                optD['dbname'] = arg
            elif opt in ("-h", "--dbhost"):
                optD['dbhost'] = arg
            elif opt in ("-p", "--dbport"):
                optD['dbport'] = arg
            elif opt == '-d':
                # DoIt!
                optD['doit'] = 1
            else:
                optD[opt] = arg
    
    if istat or not optD or not 'ns' in optD or not 'targets' in optD or not 'keys' in optD or not 'unique' in optD or not 'dbname' in optD or not 'dbhost' in optD or not 'dbport' in optD:
        print(optD)
        print('Missing required arguments.', file=sys.stderr)
        optD = list()
    return optD

def makeWhere(pkeys, source, target):
    res = ''
    for key in pkeys:
        if len(res) != 0:
            res += ' AND '
        res += source + '.' + key + ' = ' + target + '.' + key
    return res

rv = RET_SUCCESS

# Parse arguments
if __name__ == "__main__":
    optD = GetArgs(sys.argv[1:])
    if not optD:
        rv = RET_INVALIDARG
    else:
        source = optD['ns'] + '.' + optD['table']
        targets = optD['targets']
        keys = optD['keys']
        pkeys = optD['unique']
        dbuser = pwd.getpwuid(os.getuid())[0]
        dbname = optD['dbname']
        dbhost = optD['dbhost']
        dbport = optD['dbport']
        if 'doit' in optD:
            doit = 1
        else:
            doit = 0

if rv == RET_SUCCESS:
    # Connect to the database
    try:
        # The connection is NOT in autocommit mode. If changes need to be saved, then conn.commit() must be called.
        with psycopg2.connect(database=dbname, user=dbuser, host=dbhost, port=dbport) as conn:
            with conn.cursor() as cursor:
                where = makeWhere(pkeys, 'source', 'target')
                for target in targets:
                    for key in keys:
                        # UPDATE <target> AS target SET NOAA_ARS = source.NOAA_ARS FROM <source> AS source WHERE target.HARPNUM = source.HARPNUM AND target.T_REC_INDEX = source.T_REC_INDEX;
                        sql = 'UPDATE ' + target + ' AS target SET ' + key + ' = source.' + key + ' FROM ' + source + ' AS source WHERE ' + where
                        
                        if doit:
                            cursor.execute(sql)
                        else:
                            print('Executing SQL:\n  ==>' + sql)

    except psycopg2.Error as exc:
        # Closes the cursor and connection
        print(exc.diag.message_primary, file=sys.stderr)
        # No need to close cursor - leaving the with block does that.
        if not conn:
            rv = RET_DBCONNECT
        else:
            rv = RET_SQL

    # There is no need to call conn.commit() since connect() was called from within a with block. If an exception was not raised in the with block,
    # then a conn.commit() was implicitly called. If an exception was raised, then conn.rollback() was implicitly called.

sys.exit(rv)




    

