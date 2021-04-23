#!/usr/bin/env python3

# THIS IS A RENAMED VERSION OF manage-request.py. It turns out using an '-' for a python file name was a bad decision; a module name cannot have an hyphen in it, so it is not possible to import one that does; AND CVS does not permit file name changes, so, manage-request.py will be used for the legacy exportdata.html, and manage_request.py will be used for the new flask app

# The arguments to this script are parsed by cgi.FieldStorage(), which knows how to parse
# both HTTP GET and POST requests. A nice feature is that we can test the script as it runs in a CGI context
# by simply running on the command line with a single argument that is equivalent to an HTTP GET parameter string
# (e.g., address=gimli@mithril.com&addresstab=jsoc.export_addresses&domaintab=jsoc.export_addressdomains).

# Parameters:
#   address (required) - The email address to check or register.
#   addresstab (required) - The database table containing all registered (or registration-pending) email addresses.
#   domaintab (required) - The database table containing all email domains.
#   dbuser (optional) - The database account to be used when connecting to the database. The default is the value of the WEB_DBUSER parameter in DRMSParams.
#   checkonly (optional) - If set to 1, then no attept is made to register an unregistered email. In this case, if no error occurs then the possible return status codes are RV_REGISTEREDADDRESS, RV_REGISTRATIONPENDING, or RV_UNREGISTEREDADDRESS. The default is False (unknown addresses are registered).

# Testing:
#   to test providing CGI arguments, set `TEST_CGI` to True, and then provide a single QUERY_STRING-like argument, like:
#   manage-request.py 'address=person@domain&operation=check'

import sys
import os
from datetime import timedelta
import json
import psycopg2
from drmsparams import DRMSParams, DPMissingParameterError
from drmsCmdl import CmdlParser

# for use in export web app
from action import Action
class PendingRequestAction(Action):
    actions = [ 'check_pending_request', 'cancel_pending_request' ]
    def __init__(self, *, method, address):
        self._method = getattr(self, method)
        self._address = address

    def check_pending_request(self):
        # returns dict
        return check(self._address)

    def cancel_pending_request(self):
        # returns dict
        return cancel(self._address)

#test
TEST_CGI = False

# error codes
MR_ERROR_UNKNOWN = -1
MR_ERROR_PARAMETER = -2
MR_ERROR_ARGUMENT = -3
MR_ERROR_DB = -4
MR_ERROR_CHECK = -5
MR_ERROR_CANCEL = -6

# request statuses
MR_STATUS_UNKNOWN = 0
MR_STATUS_NOT_PENDING = 1
MR_STATUS_PENDING = 2
MR_STATUS_REQUEST_CANCELED = 3


# classes
class Arguments(object):
    def __init__(self, parser, args=None):
        # This could raise in a few places. Let the caller handle these exceptions.
        self.parser = parser

        # Parse the arguments.
        self.parse(args)

        # Set all args.
        self.set_all_args()

    def parse(self, args=None):
        try:
            self.parsed_args = self.parser.parse_args(args)
        except Exception as exc:
            if len(exc.args) == 2:
                type, msg = exc

                if type != 'CmdlParser-ArgUnrecognized' and type != 'CmdlParser-ArgBadformat':
                    raise # Re-raise

                raise ArgumentError(msg=msg)
            else:
                raise # Re-raise

    def __getattr__(self, name):
        # only called if object.__getattribute__(self, name) raises; and if that is true, then we want
        # to look in self.parsed_args for it, and set the instance attribute if it does exist in self.params
        value = None
        if name in vars(self.parsed_args):
            value = vars(self.parsed_args)[name]
            object.__setattr__(self, name, value)
            return value

        raise AttributeError('invalid argument ' + name)

    def set_all_args(self):
        # store in instance dict
        for name, value in vars(self.parsed_args).items():
            setattr(self, name, value)

    @classmethod
    def get_arguments(cls, program_args, drms_params):
        try:
            dbhost = drms_params.get_required('SERVER')
            dbport = drms_params.get_required('DRMSPGPORT')
            dbname = drms_params.get_required('DBNAME')
            dbuser = drms_params.get_required('WEB_DBUSER')
        except DPMissingParameterError as exc:
            raise ParameterError(str(exc))

        # check for arguments from cgi form
        args = None

        if program_args is not None and len(program_args) > 0:
            args = program_args

        parser = CmdlParser(usage='%(prog)s address=<registered email address> operation=<check, cancel> [ --dbhost=<db host> ] [ --dbport=<db port> ] [ --dbname=<db name> ] [ --dbuser=<db user>] ')
        # all arguments are considered optional in argparse (see `prefix_chars`); we can therefore do this:
        # parser.add_argument('-H', 'H', '--dbhost', ...)
        parser.add_argument('A', 'address', help='the export-registered email address', metavar='<email address>', dest='address', required=True)
        parser.add_argument('O', 'operation', help='the export-request operation to perform (check, cancel)', metavar='<operation>', dest='op', default='check', required=True)
        parser.add_argument('-H', 'H', '--dbhost', help='the host machine of the database that is used to manage pending export requests', metavar='<db host>', dest='dbhost', default=dbhost)
        parser.add_argument('-P', 'P', '--dbport', help='The port on the host machine that is accepting connections for the database', metavar='<db host port>', dest='dbport', default=dbport)
        parser.add_argument('-N', 'N', '--dbname', help='the name of the database used to manage pending export requests', metavar='<db name>', dest='dbname', default=dbname)
        parser.add_argument('-U', 'U', '--dbuser', help='the name of the database user account', metavar='<db user>', dest='dbuser', default=dbuser)

        arguments = Arguments(parser, args)

        return arguments


# exceptions
class ManageRequestError(Exception):
    _error_code = MR_ERROR_UNKNOWN

    def __init__(self, msg='generic error'):
        self._msg = msg
        self._header = None
        self._response = None

    def __str__(self):
        return self._msg

    def _generate_response(self):
        self._response = ErrorResponse(error_code=self._error_code, msg=self._header + ' ' + self._msg)

    @property
    def response(self):
        if self._response is None:
            self._generate_response()

        return self._response

class ParameterError(ManageRequestError):
    _error_code = MR_ERROR_PARAMETER

    def __init__(self, msg=None):
        super().__init__(msg=msg)

class ArgumentError(ManageRequestError):
    _error_code = MR_ERROR_ARGUMENT

    def __init__(self, msg=None):
        super().__init__(msg=msg)

class DbError(ManageRequestError):
    _error_code = MR_ERROR_DB

    def __init__(self, msg=None):
        super().__init__(msg=msg)

class CheckError(ManageRequestError):
    _error_code = MR_ERROR_CHECK

    def __init__(self, address='UNKNOWN_EXPORT_USER', msg=None):
        self._address = address

        super().__init__(msg='unable to check export request for export user ' + self._address)

        if msg is not None:
            self._msg += '(' + msg + ')'

class CancelError(ManageRequestError):
    _error_code = MR_ERROR_CANCEL

    def __init__(self, not_pending=False, address='UNKNOWN_EXPORT_USER', request_id='UNKNOWN_REQUEST_ID', start_time='UNKNOWN_START_TIME', msg=None):
        self._address = address
        self._request_id = request_id
        self._start_time = start_time

        if not_pending:
            # no error looking up address (the user has a pending request), but there was a problem deleting pending request
            super().__init__(msg='unable to cancel export request ' + self._address + ' [ request_id=' + request_id + ', start_time=' + start_time.strftime('%Y-%m-%d %T') + ' ]')
        else:
            # error looking up address (it is NOT the case that a look-up succeeded, but the address was not found)
            super().__init__(msg='unable to check export request for export user ' + self._address)

        if msg is not None:
            self._msg += '(' + msg + ')'

class OperationFactory(object):
    def __new__(cls, operation='UNKNOWN_OPERATION', address='UNKNOWN_EXPORT_USER', table='UNKNOWN_PENDING_REQUESTS_TABLE', timeout=timedelta(minutes=60)):
        if operation.lower() == CheckOperation._name:
            return CheckOperation(address, table, timeout)
        elif operation.lower() == CancelOperation._name:
            return CancelOperation(address, table, timeout)
        else:
            raise ArgumentError(msg='invalid operation type ' + operation)


# operations
class Operation(object):
    def __init__(self, address='UNKNOWN_EXPORT_USER', table='UNKNOWN_PENDING_REQUESTS_TABLE', timeout=timedelta(minutes=60)):
        self._address = address
        self._pending_requests_table = table
        self._timeout = timeout
        self._request_id = None
        self._start_time = None
        self._response = None

    def __str__(self):
        return self._name

    def process(self, cursor):
        cmd = 'SELECT request_id, start_time FROM ' + self._pending_requests_table + " WHERE address = '" + self._address + "' AND CURRENT_TIMESTAMP - start_time < interval '" + self._timeout + " minutes'"

        try:
            cursor.execute(cmd)
            rows = cursor.fetchall()
            if len(rows) > 1:
                raise DbError(msg='unexpected number of rows returned: ' + cmd)
        except psycopg2.Error as exc:
            # handle database-command errors
            import traceback
            raise DbError(msg=traceback.format_exc(8))

        if len(rows) != 0:
            self._request_id = rows[0][0]
            self._start_time = rows[0][1]

    @property
    def response(self):
        return self._response

class CheckOperation(Operation):
    _name = 'check'

    def __init__(self, address='UNKNOWN_EXPORT_USER', table='UNKNOWN_PENDING_REQUESTS_TABLE', timeout=timedelta(minutes=60)):
        super().__init__(address, table, timeout)

    def process(self, cursor):
        try:
            super().process(cursor)
        except DbError as exc:
            raise CheckError(address=self._address, msg=str(exc))

        if not self._request_id:
            self._response = NotPendingResponse(address=self._address)
        else:
            self._response = PendingResponse(address=self._address, request_id=self._request_id, start_time=self._start_time)

class CancelOperation(Operation):
    _name = 'cancel'

    def __init__(self, address='UNKNOWN_EXPORT_USER', table='UNKNOWN_PENDING_REQUESTS_TABLE', timeout=timedelta(minutes=60)):
        super().__init__(address, table, timeout)

    def process(self, cursor):
        # first run the Operation.process() code to obtain the request_id
        try:
            super().process(cursor)
        except DbError as exc:
            raise CancelError(not_pending=True, address=self._address, msg=str(exc))

        # then run the code to delete the pending request
        if not self._request_id:
            self._response = NotPendingResponse(address=self._address)
        else:
            cmd = 'DELETE FROM ' + self._pending_requests_table + " WHERE request_id = '" + self._request_id + "'"

            try:
                cursor.execute(cmd)
            except psycopg2.Error as exc:
                import traceback
                raise CancelError(not_pending=False, address=self._address, request_id=self._request_id, start_time=self._start_time, msg=traceback.format_exc(8))

            self._response = CancelResponse(address=self._address, request_id=self._request_id, start_time=self._start_time)


# responses
class ErrorResponse(Response):
    def __init__(self, error_code=None, msg=None):
        super().__init__(error_code=error_code, msg=msg)

class NotPendingResponse(Response):
    def __init__(self, address='UNKNOWN_EXPORT_USER'):
        super().__init__(error_code=MR_STATUS_NOT_PENDING, msg='no existing export request for export user ' + address)

class PendingResponse(Response):
    def __init__(self, address='UNKNOWN_EXPORT_USER', request_id='UNKNOWN_REQUEST_ID', start_time='UNKNOWN_START_TIME'):
        super().__init__(error_code=MR_STATUS_PENDING, msg='existing export request for export user ' + address + ' [ request_id=' + request_id + ', start_time=' + start_time.strftime('%Y-%m-%d %T') + ' ]')

class CancelResponse(Response):
    def __init__(self, address='UNKNOWN_EXPORT_USER', request_id='UNKNOWN_REQUEST_ID', start_time='UNKNOWN_START_TIME'):
        super().__init__(error_code=MR_STATUS_REQUEST_CANCELED, msg='existing export request for export user ' + address + ' [ request_id=' + request_id + ', start_time=' + start_time.strftime('%Y-%m-%d %T') + ' ] ' + 'was canceled')

# wrapper functions for use in a flask app
def check(address):
    resp = check_or_cancel(address=address, operation='check')
    return resp.generate_dict()

def cancel(address):
    resp = check_or_cancel(address=address, operation='cancel')
    return resp.generate_dict()

def check_or_cancel(**kwargs):
    args = None

    if len(kwargs) > 0:
        args = []
        for key, val in kwargs.items():
            args.append(key + '=' + val)

    try:
        drms_params = DRMSParams()

        if drms_params is None:
            raise ParameterError(msg='unable to locate DRMS parameters file (drmsparams.py)')

        arguments = Arguments.get_arguments(args, drms_params)

        try:
            operation = OperationFactory(operation=arguments.op, address=arguments.address, table=drms_params.EXPORT_PENDING_REQUESTS_TABLE, timeout=drms_params.EXPORT_PENDING_REQUESTS_TIME_OUT)
            resp = None

            try:
                with psycopg2.connect(host=arguments.dbhost, port=arguments.dbport, database=arguments.dbname, user=arguments.dbuser) as conn:
                    with conn.cursor() as cursor:
                        operation.process(cursor)
                        resp = operation.response
            except psycopg2.DatabaseError as exc:
                # Closes the cursor and connection
                import traceback
                raise DbError(msg='unable to connect to the database: ' + traceback.format_exc(8))
        except AttributeError as exc:
            raise ArgumentError(msg=str(exc))
    except ManageRequestError as exc:
        resp = exc.response

    return resp


if __name__ == "__main__":
    resp = check_or_cancel()

    json_response = resp.generate_json()

    # Do not print application/json here. This script may be called outside of a CGI context.
    print(json.dumps(json_response))

    # Always return 0. If there was an error, an error code (the 'status' property) and message (the 'statusMsg' property) goes in the returned HTML.
    sys.exit(0)
else:
    # stuff run when this module is loaded into another module; export things needed to call check() and cancel()
    # return json and let the wrapper layer convert that to whatever is needed by the API
    pass






# add py: base/export/scripts/checkexpdbserver, base/export/scripts/jsocextfetch, base/export/scripts/jsocextinfo, base/export/scripts/showextinfo, base/export/scripts/showextseries, base/export/scripts/seriesinfo, base/export/scripts/checkAddress
# add C programs (in drms/base): drms-export-as-fits, jsoc_info, show_info, drms_parserecset; AAAH! use securedrms.py

# dir structure - use links in code tree to create this
# proj/export/webapps/drmsexport
#   drmsexport
#     __init__.py
#     action.py
#     arguments.py
#     check_address.py
#     check_dbserver.py
#     drmsparams.py
#     error.py
#     export.py
#     export_wsgi.py
#     seriesinfo.py (like show_info/jsoc_info/showextinfo/jsocextinfo)
#     initiate_request.py (like jsoc_fetch/jsocextfetch)
#     manage_request.py
#     securedrms.py(copied from base/libs/py)
#     statuscode.py
#  setup.py
