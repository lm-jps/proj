#!/usr/bin/env python3

# This is a version of checkAddress.py modified to run inside a Flask context.

# The arguments to this script are parsed by cgi.FieldStorage(), which knows how to parse
# both HTTP GET and POST requests. A nice feature is that we can test the script as it runs in a CGI context
# by simply running on the command line with a single argument that is equivalent to an HTTP GET parameter string
# (e.g., address=gimli@mithril.com&addresstab=jsoc.export_addresses&domaintab=jsoc.export_addressdomains).

# Parameters:
#   address (required) - The email address to check or register.
#   addresstab (required) - The database table containing all registered (or registration-pending) email addresses.
#   domaintab (required) - The database table containing all email domains.
#   dbuser (optional) - The database account to be used when connecting to the database. The default is the value of the WEB_DBUSER parameter in DRMSParams.
#   checkonly (optional) - If set to 1, then no attept is made to register an unregistered email. In this case, if no error occurs then the possible return status codes are StatusCode.REGISTERED_ADDRESS, StatusCode.REGISTRATION_PENDING, or StatusCode.UNREGISTERED_ADDRESS. The default is False (unknown addresses are registered).

import sys
import os
import pwd
import re
import uuid
from urllib.parse import unquote
import smtplib
import json
import psycopg2
from drmsCmdl import CmdlParser
from drmsparams import DRMSParams
from statuscode import StatusCode as SC
from arguments import Arguments as Args
from response import Response
from error import Error as ExportError


class StatusCode(SC):
    REGISTRATION_INITIATED = 1, f'registration of {address} initiated'
    REGISTRATION_PENDING = 2, f'registration of {address} is pending'
    REGISTERED_ADDRESS = 3, f'{address} is registered'
    UNREGISTERED_ADDRESS = 4, f'{address} is not registered'

class ErrorCode(SC):
    PARAMETERS = 2, 'failure locating DRMS parameters'
    ARGUMENTS = 3, 'bad arguments'
    MAIL = 4, 'failure sending mail'
    DB = 5, 'failure executing database command'
    DB_CONNECTION = 6, 'failure connecting to database'
    DUPLICATION = 7, 'address is already registered'

class ParametersError(ExportError):
    _status_code = StatusCode(ErrorCode.PARAMETERS)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class ArgumentsError(ExportError):
    _status_code = StatusCode(ErrorCode.ARGUMENTS)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class MailError(ExportError):
    _status_code = StatusCode(ErrorCode.MAIL)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class DBError(ExportError):
    _status_code = StatusCode(ErrorCode.DB)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class DBConnectionError(ExportError):
    _status_code = StatusCode(ErrorCode.DB_CONNECTION)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class DuplicationError(ExportError):
    _status_code = StatusCode(ErrorCode.DUPLICATION)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class RegistrationInitiatedResponse(Response):
    def __init__(self, *, address):
        super().__init__(status_code=StatusCode.REGISTRATION_INITIATED)

class RegistrationPendingResponse(Response):
    def __init__(self, *, address):
        super().__init__(status_code=StatusCode.REGISTRATION_PENDING)

class RegisteredResponse(Response):
    def __init__(self, *, address):
        super().__init__(status_code=StatusCode.REGISTERED_ADDRESS)

class UnregisteredResponse(Response):
    def __init__(self, *, address):
        super().__init__(status_code=StatusCode.UNREGISTERED_ADDRESS)

class UnquoteAction(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        unquoted = unquote(value)
        setattr(namespace, self.dest, unquoted)

class Arguments(Args):
    _arguments = None

    @classmethod
    def get_arguments(cls, *, program_args, drms_params):
        if cls._arguments is None:
            try:
                dbhost = drms_params.get_required('SERVER')
                dbport = int(drms_params.get_required('DRMSPGPORT'))
                dbname = drms_params.get_required('DBNAME')
                dbuser = drms_params.get_required('WEB_DBUSER')
                address_info_fn = drms_params.get_required('EXPORT_ADDRESS_INFO_FN')
                address_info_insert_fn = drms_params.get_required('EXPORT_ADDRESS_INFO_INSERT_FN')
                user_info_fn = drms_params.get_required('EXPORT_USER_INFO_FN')
                user_info_insert_fn = drms_params.get_required('EXPORT_USER_INFO_INSERT_FN')
                regemail_timeout = drms_params.get_required('REGEMAIL_TIMEOUT')
            except DPMissingParameterError as exc:
                raise ParametersError(msg=str(exc))

            args = None

            if program_args is not None and len(program_args) > 0:
                args = program_args

            parser = CmdlParser(usage='%(prog)s address=<email address to register/check> operation=<register/check> [ --name=<user\'s name> ] [ --snail=<user snail mail address> ] [ --dbhost=<db host> ] [ --dbport=<db port> ] [ --dbname=<db name> ] [ --dbuser=<db user>] ')
            # all arguments are considered optional in argparse (see `prefix_chars`); we can therefore do this:
            # parser.add_argument('-H', 'H', '--dbhost', ...)
            parser.add_argument('a', 'address', help='the email address to register or check', metavar='<email address>', dest='address', action=UnquoteAction, required=True)
            parser.add_argument('o', 'operation', help='the operation: register or check', metavar='<operation>', choices=[ 'register', 'check'], dest='operation', required=True)

            parser.add_argument('-n', '--name', help='the user name to register', metavar='<export user\'s name>', dest='name', action=UnquoteAction, default='NULL')
            parser.add_argument('-s', '--snail', help='the user snail-mail address to register', metavar='<export user\'s snail mail>', dest='snail', action=UnquoteAction, default='NULL')
            parser.add_argument('-H', '--dbhost', help='the host machine of the database that is used to manage pending export requests', metavar='<db host>', dest='dbhost', default=dbhost)
            parser.add_argument('-P', '--dbport', help='The port on the host machine that is accepting connections for the database', metavar='<db host port>', dest='dbport', type=int, default=int(dbport))
            parser.add_argument('-N', '--dbname', help='the name of the database used to manage pending export requests', metavar='<db name>', dest='dbname', default=dbname)
            parser.add_argument('-U', '--dbuser', help='the name of the database user account', metavar='<db user>', dest='dbuser', default=dbuser)

            arguments = Arguments(parser=parser, args=args)
            arguments.set_arg('address_info_fn', address_info_fn)
            arguments.set_arg('address_info_insert_fn', address_info_insert_fn)
            arguments.set_arg('user_info_fn', user_info_fn)
            arguments.set_arg('user_info_insert_fn', user_info_insert_fn)
            arguments.set_arg('regemail_timeout', regemail_timeout)

            # Do a quick validation on the email address.
            reg_exp = re.compile(r'\s*[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,6}')
            match_obj = reg_exp.match(arguments.address)
            if match_obj is None:
                raise ArgumentsError(f'{arguments.address} is not a valid email address')

            cls._arguments = arguments

        return cls._arguments

def SendMail(address, timeout, confirmation):
    subject = 'CONFIRM EXPORT ADDRESS'
    fromAddr = 'jsoc@solarpost.stanford.edu'
    toAddrs = [ address ]
    bccAddrs = [ 'art.amezcua@stanford.edu' ]
    msg = 'From: ' + fromAddr + '\nTo: ' + ','.join(toAddrs) + '\nSubject: ' + subject + '\nThis message was automatically generated by the JSOC export system at Stanford.\n\nYou have requested that data be exported from the JSOC. To do so, you must register your email address with the export system. To complete the registration process, please reply to this message within ' + timeout + ' minutes. Please do not modify the body of this message when replying. The server will extract the embedded confirmation code to verify that the provided email address is valid. You will receive another email message notifying you of the disposition of your registration.'
    msg += '\n[' + str(confirmation) + ']'

    try:
        server = smtplib.SMTP('solarpost.stanford.edu')
        server.sendmail(fromAddr, toAddrs + bccAddrs, msg)
        server.quit()
    except Exception as exc:
        # If any exception happened, then the email message was not received.
        raise MailError(f'unable to send email message to {','.join(toAddrs)} to confirm address')

def get_arguments(**kwargs):
    args = None

    if len(kwargs) > 0:
        args = []
        for key, val in kwargs.items():
            args.append(f'{key} = {val}'')

    drms_params = DRMSParams()

    if drms_params is None:
        raise ParametersError(msg=f'unable to locate DRMS parameters file')

    return Arguments.get_arguments(program_args=args, drms_params=drms_params)

def generate_registered_or_pending_response(arguments, confirmation):
    response = None

    # get requestor ID also
    cmd = f'SELECT id FROM {arguments.user_info_fn}(\'{arguments.address}\')'
    try:
        cursor.execute(cmd)
        rows = cursor.fetchall()
    except psycopg2.Error as exc:
        # handle database-command errors
        raise DBError(f'{str(exc)}')

    user_id = rows[0][0]

    if confirmation is None or len(confirmation) == 0:
        # if confirmation == None ==> registered
        response = Response.generate_response(status_code=StatusCode.REGISTERED_ADDRESS, user_id=user_id)
    else:
        # if confirmation != None ==> pending registration
        response = Response.generate_response(status_code=StatusCode.REGISTRATION_PENDING, user_id=user_id)

    return response

def check(**kwargs):
    response = None

    try:
        arguments = get_arguments(kwargs)

        try:
            with psycopg2.connect(database=arguments.dbname, user=arguments.dbuser, host=arguments.dbhost, port=str(arguments.dbport)) as conn:
                with conn.cursor() as cursor:
                    cmd = f'SELECT confirmation FROM {arguments.address_info_fn}(\'{arguments.address}\')'

                    try:
                        cursor.execute(cmd)
                        rows = cursor.fetchall()

                        if len(rows) == 0:
                            # the address is not in the db, and the user did not request registration ==> not registered
                            response = Response.generate_response(status_code=StatusCode.UNREGISTERED_ADDRESS, user_id=-1)
                        elif len(rows) == 1:
                            # the address is in the db
                            confirmation = rows[0][0]
                        else:
                            raise DBError(f'unexpected number of rows returned: {cmd}')
                    except psycopg2.Error as exc:
                        # handle database-command errors
                        raise DBError(f'{str(exc)}')

                    if response is None:
                        response = generate_registered_or_pending_response(arguments, confirmation)
        except psycopg2.OperationalError as exc:
            # closes the cursor and connection
            raise DBConnectionError(f'unable to connect to the database: {str(exc)}')
    except ExportError as exc:
        response = exc.response

    return response

def check_and_register(**kwargs):
    try:
        arguments = get_arguments(kwargs)

        try:
            with psycopg2.connect(database=arguments.dbname, user=arguments.dbuser, host=arguments.dbhost, port=str(arguments.dbport)) as conn:
                with conn.cursor() as cursor:
                    cmd = f'SELECT confirmation FROM {arguments.address_info_fn}(\'{arguments.address}\')'

                    try:
                        cursor.execute(cmd)
                        rows = cursor.fetchall()

                        if len(rows) == 0:
                            # the address is not in the db, and the user did request registration ==> registration initiated
                            confirmation = uuid.uuid4()

                            # ensure confirmation does not already exist in addresses_tab
                            cmd = f'SELECT confirmation FROM {arguments.address_info_fn}() WHERE confirmation = \'{str(confirmation)}\''

                            try:
                                cursor.execute(cmd)
                                rows = cursor.fetchall()
                                if len(rows) > 0:
                                    raise DuplicationError(f'cannot insert row into address table; confirmation {str(confirmation)} already exists')
                            except psycopg2.Error as exc:
                                # Handle database-command errors.
                                raise DBError(f'{str(exc)}')

                            # insert into the addresses table (and domains table if need be)
                            cmd = f'SELECT * FROM {arguments.address_info_insert_fn}(\'{arguments.address}\', \'{str(confirmation)}\')'
                            try:
                                cursor.execute(cmd)
                            except psycopg2.Error as exc:
                                # Handle database-command errors.
                                raise DBError(f'{str(exc)}')

                            # we have to also insert into the export user table since we have that information now, not after the user has replied to the registration email
                            # (which is processed by registerAddress.py); if a failure happens anywhere along the way, we need to delete the entry from the export user table
                            cmd = f'SELECT id FROM {arguments.user_info_insert_fn}(\'{arguments.address}\', \'{arguments.name}\', \'{arguments.snail}\')'

                            try:
                                cursor.execute(cmd)
                                rows = cursor.fetchall()
                            except psycopg2.Error as exc:
                                # Handle database-command errors.
                                raise DBError(f'{str(exc)}')

                            user_id = rows[0][0]

                            # send an email message out with a new confirmation code
                            SendMail(arguments.address, arguments.regemail_timeout, confirmation)

                            msg = f'Your email address is being registered for use with the export system. You will receive an email message from user jsoc. Please reply to this email message within {arguments.regemail_timeout} minutes without modifying the body.'

                            response = Response.generate_response(status_code=StatusCode.REGISTRATION_INITIATED, msg=msg, user_id=user_id)
                        elif len(rows) == 1:
                            # the address is in the db
                            confirmation = rows[0][0]
                        else:
                            raise DBError(f'unexpected number of rows returned: {cmd}')
                    except psycopg2.Error as exc:
                        # handle database-command errors
                        raise DBError(f'{str(exc)}')

                    if response is None:
                        response = generate_registered_or_pending_response(arguments, confirmation)
        except psycopg2.OperationalError as exc:
            # closes the cursor and connection
            raise DBConnectionError(f'unable to connect to the database: {str(exc)}')

    except ExportError as exc:
        response = exc.response

    return response

if __name__ == "__main__":
    try:
        arguments = get_arguments(kwargs)

        if arguments.operation == 'register':
            response = check_and_register()
        else:
            response = check()

    except ExportError as exc:
        response = exc.response

    print(response.generate_json())

    # Always return 0. If there was an error, an error code (the 'status' property) and message (the 'statusMsg' property) goes in the returned HTML.
    sys.exit(0)
