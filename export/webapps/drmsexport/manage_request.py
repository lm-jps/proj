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
from drms_export import Error as ExportError, ErrorCode as ExportErrorCode, Response
from drms_parameters import DRMSParams, DPMissingParameterError
from drms_utils import Arguments as Args, CmdlParser, StatusCode as SC


#test
TEST_CGI = False

# error codes
MR_ERROR_UNKNOWN = -1
MR_ERROR_PARAMETER = -2
MR_ERROR_ARGUMENT = -3
MR_ERROR_DB = -4
MR_ERROR_CHECK = -5
MR_ERROR_CANCEL = -6


class StatusCode(SC):
    NOT_PENDING = 1, 'request {id} is not pending'
    PENDING = 2, 'request {id} is pending'
    REQUEST_CANCELED = 3, 'pending request {id} was canceled'

class ErrorCode(ExportErrorCode):
    PARAMETERS = 1, 'failure locating DRMS parameters'
    ARGUMENTS = 2, 'bad arguments'
    DB = 3, 'failure executing database command'
    DB_CONNECTION = 4, 'failure connecting to database'
    CHECK = 5, 'unable to check export request for export user {address}'
    CANCEL = 6, 'unable to cancel export request for export user {address}'


# exceptions
class ParametersError(ExportError):
    _error_code = ErrorCode(ErrorCode.PARAMETERS)
    # _header = f'if present, then `[cls.header]` will appear at the beginning of the error message'

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class ArgumentsError(ExportError):
    _error_code = ErrorCode(ErrorCode.ARGUMENTS)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class DBError(ExportError):
    _error_code = ErrorCode(ErrorCode.DB)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class DBConnectionError(ExportError):
    _error_code = ErrorCode(ErrorCode.DB_CONNECTION)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class CheckError(ExportError):
    _error_code = ErrorCode(ErrorCode.CHECK)

    def __init__(self, *, address, msg=None):
        error_msg = self._error_code.fullname(address=address)
        if msg is not None and len(msg) > 0:
            error_msg = f'{error_msg}: {msg}'
        super().__init__(msg=error_msg)

class CancelError(ExportError):
    _error_code = ErrorCode(ErrorCode.CANCEL)

    def __init__(self, *, address, request_id, start_time, msg=None):
        # not_pending == True ==> cannot find a pending request for `address`
        error_msg = self._error_code.fullname(address=address)
        if msg is not None and len(msg) > 0:
            error_msg = f'{error_msg}: {msg}'
        self._msg = error_msg

    def __init__(self, not_pending=False, address='UNKNOWN_EXPORT_USER', request_id='UNKNOWN_REQUEST_ID', start_time='UNKNOWN_START_TIME', msg=None):
        self._address = address
        self._request_id = request_id
        self._start_time = start_time

        if not_pending:
            # no error looking up address (the user has a pending request), but there was a problem deleting pending request
            super().__init__(msg=f'unable to cancel export request for user {self._address}')
        else:
            # error looking up address (it is NOT the case that a look-up succeeded, but the address was not found)
            super().__init__(msg=f'unable to check export request for user {self._address}')

        if msg is not None:
            self._msg += f'({msg})'


# classes
class Arguments(Args):
    _arguments = None

    @classmethod
    def get_arguments(cls, program_args, drms_params):
        if cls._arguments is None:
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

            cls._arguments = arguments

        return cls._arguments

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
                raise DbError(msg=f'unexpected number of rows returned: {cmd}')
        except psycopg2.Error as exc:
            # handle database-command errors
            raise DbError(msg=str(exc))

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
            # cannot locate the pending request in the pending-requests db table
            raise CancelError(address=self._address, msg=f'cannot locate pending request in database ({str(exc)})')

        # then run the code to delete the pending request
        if not self._request_id:
            self._response = NotPendingResponse(address=self._address)
        else:
            cmd = 'DELETE FROM ' + self._pending_requests_table + " WHERE request_id = '" + self._request_id + "'"

            try:
                cursor.execute(cmd)
            except psycopg2.Error as exc:
                # cannot delete the pending request from the pending-requests db table
                error_msg = f'cannot delete pending request with id={self._request_id} and start_time={self._start_time.strftime("%Y-%m-%d %T")} ({str(exc)})'
                raise CancelError(address=self._address, msg=error_msg)

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

# for use in export web app
from action import Action
class PendingRequestAction(Action):
    actions = [ 'check_pending_request', 'cancel_pending_request' ]
    def __init__(self, *, method, address, dbhost=None, dbport=None, dbname=None, dbuser=None):
        self._method = getattr(self, method)
        self._address = address
        self._options = {}
        self._options['dbhost'] = dbhost
        self._options['dbport'] = dbport
        self._options['dbname'] = dbname
        self._options['dbuser'] = dbuser

    def check_pending_request(self):
        # returns dict
        response = perform_action(operation='check', address=self._address, options=self._options)
        return response.generate_dict()

    def cancel_pending_request(self):
        # returns dict
        response = perform_action(operation='register', address=self._address, options=self._options)
        return response.generate_dict()

def perform_action(**kwargs):
    args = []

    for key, val in kwargs.items():
        if val is not None:
            if key == 'options':
                for option, option_val in val.items():
                    args.append(f'--{option}={option_val}')
            else:
                args.append(f'{key}={val}')

    try:
        drms_params = DRMSParams()

        if drms_params is None:
            raise ParameterError(msg='unable to locate DRMS parameters file (drmsparams.py)')

        arguments = Arguments.get_arguments(args, drms_params)

        try:
            operation = OperationFactory(operation=arguments.op, address=arguments.address, table=drms_params.EXPORT_PENDING_REQUESTS_TABLE, timeout=drms_params.EXPORT_PENDING_REQUESTS_TIME_OUT)
            response = None

            try:
                with psycopg2.connect(host=arguments.dbhost, port=arguments.dbport, database=arguments.dbname, user=arguments.dbuser) as conn:
                    with conn.cursor() as cursor:
                        operation.process(cursor)
                        response = operation.response
            except psycopg2.OperationalError as exc:
                # closes the cursor and connection
                raise DBConnectionError(f'unable to connect to the database: {str(exc)}')
        except AttributeError as exc:
            raise ArgumentError(msg=str(exc))
    except ExportError as exc:
        response = exc.response

    return response


if __name__ == "__main__":
    response = perform_action()
    json_response = response.generate_json()

    # do not print application/json here; this script may be called outside of a web abb context
    print(json_response)

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
