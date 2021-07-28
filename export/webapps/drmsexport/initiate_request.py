#!/usr/bin/env python3

# basically calls jsoc_fetch

from argparse import Action as ArgsAction
from drms_export import Error as ExportError, ErrorCode as ExportErrorCode, Response, SecureClientFactory
from drms_parameters import DRMSParams, DPMissingParameterError
from drms_utils import Arguments as Args, CmdlParser, MakeObject, StatusCode as SC


FILE_NAME_SIZE = 256

class StatusCode(SC):
    SUCCESS = 0, 'success'
    REQUEST_NOT_FOUND = 1, 'request not found'
    REQUEST_QUEUED = 2, 'request is in queue'
    REQUEST_PROCESSING = 3, 'request is being processed by manager'
    REQUEST_COMPLETE = 4, 'request has been serviced'

    @classmethod
    def get_status_code(cls, request_status):
        if request_status == 6:
            # this should go away - if in jsoc.export_new, but not in jsoc.export, then this is the same as REQUEST_QUEUED
            return StatusCode.REQUEST_NOT_FOUND
        if request_status == 2  or request_status == 12:
            return StatusCode.REQUEST_QUEUED
        if request_status == 1:
            return StatusCode.REQUEST_PROCESSING
        if request_status == 0:
            return StatusCode.REQUEST_COMPLETE

class ErrorCode(ExportErrorCode):
    PARAMETERS = 1, 'failure locating DRMS parameters'
    ARGUMENTS = 2, 'bad arguments'
    DRMS_CLIENT = 3, 'drms client error'
    STREAM_FORMAT = 4, 'format error in downloaded payload'
    REQUEST_REJECTED = 5, 'too many requests'
    REQUEST_NOT_ONLINE = 6, 'request no longer online'
    REQUEST_FATAL = 7, 'fatal error'
    REQUEST_TOO_LARGE = 8, 'request too large'

    @classmethod
    def get_error_code(cls, request_error):
        if request_status == 7:
            return ErrorCode.REQUEST_REJECTED
        if request_status == 5:
            return ErrorCode.REQUEST_NOT_ONLINE
        if request_status == 4:
            return ErrorCode.REQUEST_FATAL
        if request_status == 3:
            return ErrorCode.REQUEST_TOO_LARGE


class ParametersError(ExportError):
    _error_code = ErrorCode(ErrorCode.PARAMETERS)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class ArgumentsError(ExportError):
    _error_code = ErrorCode(ErrorCode.ARGUMENTS)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class DRMSClientError(ExportError):
    _error_code = ErrorCode(ErrorCode.DRMS_CLIENT)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class StreamFormatError(ExportError):
    _error_code = ErrorCode(ErrorCode.STREAM_FORMAT)

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class ExportTypeAction(ArgsAction):
    def __call__(self, parser, namespace, value, option_string=None):
        if value == 'premium':
            setattr(namespace, 'package', None)
        elif value == 'mini':
            setattr(namespace, 'processing', None)
            setattr(namespace, 'package', None)
            setattr(namespace, 'access', 'http')
            setattr(namespace, 'file_format', 'native')
        else:
            setattr(namespace, 'access', 'streamed')

        setattr(namespace, self.dest, value)

class ExportArgumentsAction(ArgsAction):
    def __call__(self, parser, namespace, value, option_string=None):
        # convert json arguments to dict
        arguments_dict = json.loads(value)
        if 'specification' not in arguments_dict:
            raise ArgumentsError(f'missing required export argument `specification`')
        setattr(namespace, self.dest, arguments_dict)

def create_webserver_action(drms_params):
    class WebserverAction(ArgsAction):
        def __call__(self, parser, namespace, value, option_string=None):
            webserver_dict = {}
            webserver_dict['host'] = value
            webserver_dict['public'] = True if webserver['public'].lower() != drms_params.get_required('WEB_DOMAIN_PRIVATE') else False
            webserver_obj = MakeObject('webserver', webserver_dict)
            setattr(namespace, self.dest, arguments_dict)

class Arguments(Args):
    _arguments = None

    @classmethod
    def get_arguments(cls, *, program_args, drms_params):
        if cls._arguments is None:
            try:
                db_host = drms_params.get_required('SERVER')
                db_port = int(drms_params.get_required('DRMSPGPORT'))
                db_name = drms_params.get_required('DBNAME')
                db_user = drms_params.get_required('WEB_DBUSER')
            except DPMissingParameterError as exc:
                raise ParametersError(msg=str(exc))

            args = None

            if program_args is not None and len(program_args) > 0:
                args = program_args

            parser = CmdlParser(usage='%(prog)s [ --dbhost=<db host> ] [ --dbport=<db port> ] [ --dbname=<db name> ] [ --dbuser=<db user>] [ --skiptar=<boolean value> ] [ --compression=<boolean value ] [ --file_name_format=<> ] webserver=<> spec=<>')

            # optional
            parser.add_argument('-r', '--requestor', help='the name of the export user', metavar='<requestor>', dest='requestor', default=None)
            parser.add_argument('-H', '--dbhost', help='the machine hosting the database that contains export requests from this site', metavar='<db host>', dest='db_host', default=dbhost)
            parser.add_argument('-P', '--dbport', help='the port on the host machine that is accepting connections for the database', metavar='<db host port>', dest='db_port', type=int, default=int(dbport))
            parser.add_argument('-N', '--dbname', help='the name of the database that contains export requests', metavar='<db name>', dest='db_name', default=dbname)
            parser.add_argument('-U', '--dbuser', help='the name of the database user account', metavar='<db user>', dest='dbuser', default=db_user)

            # required
            parser.add_argument('export_type', help='the export type: premium, mini, or streamed', metavar='<webserver>', choices=Choices(['premium', 'mini', 'streamed']), action=ExportTypeAction, dest='export_type', required=True)
            parser.add_argument('webserver', help='the webserver invoking this script', metavar='<webserver>', action=create_webserver_action(drms_params), dest='webserver', required=True)
            parser.add_argument('address', help='the email addressed registered for export', metavar='<email address>', dest='address', required=True)
            parser.add_argument('arguments', help='export arguments', action=ExportArgumentsAction, dest='export_arguments', reqired=True)

            arguments = Arguments(parser=parser, args=args)
            cls._arguments = arguments

        return cls._arguments

def jsonify(data_frame):
    # input:
    #         record                    filename
    #   0   hmi.M_720s[2017.1.23]       hmi__m_720s_20170123.fits
    #   1   hmi.M_720s[2017.1.24]       hmi__m_720s_20170124.fits
    #   2   hmi.M_720s[2017.1.25]       hmi__m_720s_20170125.fits
    #
    # output:
    # [
    #    {
    #      "record" : "hmi.M_720s[2017.1.23]",
    #      "filename" : "hmi__m_720s_20170123.fits"
    #    },
    #    {
    #      "record" : "hmi.M_720s[2017.1.24]",
    #      "filename" : "hmi__m_720s_20170124.fits"
    #    },
    #    {
    #      "record" : "hmi.M_720s[2017.1.25]",
    #      "filename" : "hmi__m_720s_20170125.fits"
    #    }
    # ]
    return data_frame.to_json(orient='records')

def export_premium(*, drms_client, address, requestor=None, export_arguments):
    '''
    export_arguments:
    {
      specification : <DRMS record-set specification>,
      processing : <dict of processing argument>,
      package : <type and file name>,
      access : <http, ftp>,
      file_format <exported file type>
      file_name_format : <format string for name of exported file>,
      size_ratio : <multiplier to obtain image size from original image size>
      number_records : <maximum number of records exported>
    }
    '''

    try:
        if kwargs['package'].type == 'tar':
            if kwargs['access'] == 'ftp':
                method = 'ftp-tar'
            else:
                method = 'url-tar'
        else:
            method = 'url'

        request = drms_client.export(ds=specification, email=address, requestor=requestor, method=method, protocol=kwargs['file_format'], protocol_args=kwargs['protocol_args'], filename_fmt=kwargs['file_name_format'], n=kwargs['number_records'], synchronous_export=False)
    except SecureDRMSError as exc:
        raise DRMSClientError(str(exc))

    if (status_code == StatusCode.REQUEST_COMPLETE):
        data = jsonify(request.data)
        response = Response.generate_response(status_code=StatusCode.SUCCESS, requestid=request.id, status=request.status, method=request.method, dir=request.dir, data=data, count=len(data), size=request.size)
    else:
        response = ErrorResponse.generate_response(status_code=status_code, requestid=request.id, status=request.status, method=request.method, dir=request.dir, data=None, count=len(data), size=request.size, error=request.error_msg, contact=request.contact)

    return response

def export_mini(*, drms_client, address, specification, export_arguments, requestor=None, **kwargs):
    '''
    export_arguments:
    {
      specification : <DRMS record-set specification>,
      file_name_format : <format string for name of exported file>,
      size_ratio : <multiplier to obtain image size from original image size>
      number_records : <maximum number of records exported>
    }
    '''

    try:
        # export() runs the export job asynchronously, returning a SecureExportRequest instance; the following attributes are returned and used by the export web page to display a table of links to the exported data files:
        #   ExportRequest.data - a pandas.DataFrame with two columns: record (a record spec), and filename (the name of the file in the export SU)
        #   ExportRequest.method - the original export type (url, url_direct, url_quick, ftp, url-tar, ftp-tar); totally not needed since the only allowable `method` is 'url_quick'
        #   ExportRequest.dir - the export SU (the SUMS SU that contains the export files)

        request = drms_client.export(ds=specification, email=address, requestor=requestor, method='url_quick', protocol='as-is', protocol_args=kwargs['protocol_args'], filename_fmt=swargs['file_name_format'], n=['number_records'], synchronous_export=False)

        # XXX since the original drms package does not provide attributes other than request_id and status, override ExportRequest in SecureExportRequest to provide the other needed bits of information (size, error, contact); these three attributes do not necessarily exist; SecureExportRequest._d has that information
    except SecureDRMSError as exc:
        raise DRMSClientError(str(exc))

    status_code = get_status_code(request.status)

    if (status_code == StatusCode.REQUEST_COMPLETE):
        data = jsonify(request.data)
        response = Response.generate_response(status_code=StatusCode.SUCCESS, requestid=request.id, status=request.status, method=request.method, dir=request.dir, data=data, count=len(data), size=request.size)
    else:
        # not all statuses are possible (a premium request will not be initiated)
        response = ErrorResponse.generate_response(status_code=status_code, requestid=request.id, status=request.status, method=request.method, dir=request.dir, data=None, count=len(data), size=request.size, error=request.error_msg, contact=request.contact)

    return response

def export_streamed(*, drms_client, address, specification, export_arguments, requestor=None, **kwargs):
    '''
    export_arguments:
    {
      specification : <DRMS record-set specification>,
      file_name_format : <format string for name of exported file>,
      size_ratio : <multiplier to obtain image size from original image size>
    }
    '''
    try:
        download_stream = io.BytesIO()

        # spawn a thread to download and write to `download_stream`; will return while streaming is occurring;
        # drms-export-to-stdout will verify email address
        drms_client.export_package(email=address, spec=specificaton, filename_fmt=kwargs['file_name_fmt'], synchronous=download_stream)

        file_name = ''

        # a FITS file prepended with filename; read the filename first
        while len(file_name) < FILE_NAME_SIZE:
            pipe_bytes = download_stream.read(FILE_NAME_SIZE - len(file_name))

            if len(pipe_bytes) == 0:
                raise StreamFormatError('improper formatting for FITS-file-name header')
            file_name = file_name + pipe_bytes.decode('UTF8')

        # truncate padding (0 bytes)
        file_name = file_name.rstrip('\x00')

        if len(file_name) == 0:
            # FITS file, which may have a filename, or it may not
            file_name = DUMMY_FITS_FILE_NAME

        # force alphanumeric (plus '_'), preserving the file extension
        matches = re.search(r'(^.+)[.]([^.]+$)', file_name)
        if matches is not None:
            base = matches.group(1)
            extension = matches.group(2)
            file_name = re.sub(r'[^a-zA-Z0-9_]', '_', base) + '.' + extension

        # we are streaming either a single FITS file; write HTTP header
        sys.stdout.buffer.write(b'Content-type: application/octet-stream\n')
        sys.stdout.buffer.write(b'Content-Disposition: attachment; filename="' + file_name.encode() + b'"\n')
        sys.stdout.buffer.write(b'Content-transfer-encoding: binary\n\n')

        # dump the remainder of the FITS file
        sys.stdout.buffer.write(download_stream.read())
        sys.stdout.buffer.flush()
    except SecureDRMSError as exc:
        raise DRMSClientError(str(exc))

    response = Response.generate_response(status_code=StatusCode.SUCCESS)
    return response

def get_arguments(**kwargs):
    args = []
    for key, val in kwargs.items():
        args.append(f'{key} = {val}')

    drms_params = DRMSParams()

    if drms_params is None:
        raise ParametersError(msg=f'unable to locate DRMS parameters file')

    return Arguments.get_arguments(program_args=args, drms_params=drms_params)

def perform_action(**kwargs):
    response = None

    try:
        arguments = get_arguments(kwargs)

        try:
            factory = SecureClientFactory(debug=False)

            connection_info = { 'dbhost' : arguments.db_host, 'dbport' : arguments.db_port, 'dbname' : arguments.db_name, 'dbuser' : arguments.db_user }
            sshclient_internal = factory.create_ssh_client(use_internal=True)
            sshclient_external = factory.create_ssh_client(use_internal=False)

            # this operates on internal DB - external DB is a subset of internal, so the former is appropriate
            spec_info = sshclient_internal.parse_spec(arguments.specification)

            # if this script is invoked from the external web server, determine the db host that can handle all series in the request; some series might be 'pass-through' series
            has_passthrough_series = False
            if arguments.webserver.public:
                # `arguments.dbhost` is the DB server for the public site (`arguments.webserver.host` is the httpd host, e.g., jsoc.stanford.edu)
                series = []
                for subset in spec_info.subsets:
                    series.append(subset.seriesname)

                # call check_dbserver.py code
                action = Action.action(action_type='determine_db_server', args={ "--dbport" : str(dbport), "--dbname" : dbname, "--dbuser" : dbuser, "public_dbhost" : arguments.dbhost, "series" : json.dumps(series) })
                response = action()
                server_info = response.generate_json_obj()

                # map server_info.server to internal/external DB server and use the correct server with jsoc_info; we need to create both an internal drms_client and an external one, and we use the external one if server_info.server maps to the external server, and use the internal one otherwise; when configuring securedrms.py, get rid of jsocextfetch.py - just use jsoc_fetch always, but use the correct client to set the JSOC_DBHOST argument to jsoc_info correctly

                if server_info.public:
                    drms_client = sshclient_external
                else:
                    drms_client = sshclient_internal

            # there are four potential attributes an export request can have:
            #   processing (processed ==> original image is modified); NA for keyword-only exports
            #   package (tar ==> exported files are packaged into a tar file)
            #   access (ftp ==> exported produce is accessible from ftp server; http ==> exported produce is accessible from http server; stream ==> data are streamed back to user)
            #   file format (NATIVE, FITS, JPEG, MPEG, MP4); NA for keyword-only exports
            # depending on the operation selected, some of these attributes may be determined and immutable; also, some file-format conversions are not possible (for example, if an export contains metadata, but no image data, then none of the image file formats are relevant)

            # possible export arguments (from webapp):
            #   specification [RecordSet] (record-set specification)
            #   address [Notify] (registered email address)
            #   requestor [Requestor] export user's name, either entered in text edit, or from cookie
            #   export_type [Method - url, url_direct, url_quick, ftp, url-tar, ftp-tar] these map to...
            #     url ==> export_premium (no tar)
            #     url_direct ==> export_streamed
            #     url_quick ==> export_mini
            #     ftp ==> export_premium (ftp-access)
            #     url-tar ==> export_premium (tar)
            #     ftp-tar ==> export_premium (tar, ftp access)
            #   processing_info [Processing] (processing-program-specific arguments)
            #   tar [from Method (name contains '-tar' string)]
            #   access [from Method]
            #   file_format [Protocol - as-is, FITS, JPEG, MPEG, MP4] (exported file type)
            #   file_name_format [Filename Format]
            #   size_ratio [not shown] determined by code in processing.js plus keyword metadata

            export_arguments = arguments.arguments # dict

            if arguments.export_type == 'premium':
                # supports all export four export attributes; use securedrms.SecureClient.export(method='url')
                response = export_premium(drms_client=drms_client, address=arguments.address, requestor=arguments.requestor, **export_arguments)
            elif arguments.export_type == 'mini_export':
                # supports no processing, no tar, http access, native file format attributes; use securedrms.SecureClient.export(method='url_quick')
                response = export_mini(drms_client=sshclient, specification=arguments.specification, address=arguments.address, requestor=arguments.requestor, **export_arguments)
            elif arguments.export_type == 'streamed_export':
                # supports no-processing, no-tar, stream-access, native file format attributes; exports a single file only, all SUs must be online; use securedrms.SecureClient.export_package()
                response = export_streamed(drms_client=sshclient, specification=arguments.specification, address=arguments.address, requestor=arguments.requestor, **export_arguments)

        except SecureDRMSError as exc:
            raise DRMSClientError(str(exc))
    except ExportError as exc:
        response = exc.response

    return response


# for use in export web app
from action import Action
class InitiateRequestAction(Action):
    actions = [ 'start_premium_export', 'start_mini_export', 'start_streamed_export' ]
    def __init__(self, *, method, webserver, address, export_arguments, requestor=None, dbhost=None, dbport=None, dbname=None, dbuser=None):
        self._method = getattr(self, method)
        self._webserver = webserver # dict
        self._address = adderss
        self._export_arguments = export_arguments
        self._options = {}
        self._options['requestor'] = requestor
        self._options['dbhost'] = dbhost
        self._options['dbport'] = dbport
        self._options['dbname'] = dbname
        self._options['dbuser'] = dbuser

    def start_premium_export(self):
        '''
        export_arguments:
        {
          specification : <DRMS record-set specification>,
          processing : <dict of processing argument>,
          package : <type and file name>,
          access : <http, ftp>,
          file_format <exported file type>,
          file_name_format : <format string for name of exported file>,
          size_ratio : <multiplier to obtain image size from original image size>,
          number_records : <maximum number of records exported>
        }
        '''
        response = perform_action(export_type='premium', webserver=self._webserver, address=self._address, export_arguments=self._export_arguments, options=self._options)
        return response

    def start_mini_export(self):
        '''
        export_arguments:
        {
          specification : <DRMS record-set specification>,
          file_name_format : <format string for name of exported file>,
          size_ratio : <multiplier to obtain image size from original image size>
          number_records : <maximum number of records exported>
        }
        '''
        response = perform_action(export_type='mini', webserver=self._webserver, address=self._address, export_arguments=self._export_arguments, options=self._options)
        return response

    def start_streamed_export(self):
        # returns dict
        response = perform_action(export_type='streamed', webserver=self._webserver, address=self._address, export_arguments=self._export_arguments, options=self._options)
        return response





if __name__ == "__main__":
    try:
        response = perform_action()
    except ExportError as exc:
        response = exc.response

    print(response.generate_json())

    # Always return 0. If there was an error, an error code (the 'status' property) and message (the 'statusMsg' property) goes in the returned HTML.
    sys.exit(0)
else:
    pass
