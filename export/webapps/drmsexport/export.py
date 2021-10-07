#!/usr/bin/env python3

from flask import Flask, request
from flask_restful import Api, Resource
from logging import DEBUG as LOGGING_LEVEL_DEBUG, INFO as LOGGING_LEVEL_INFO
from urllib.parse import urlparse
from webargs import fields, validate, flaskparser, ValidationError
from webargs.flaskparser import use_args, use_kwargs, parser, abort

from action import Action
from drms_export import ErrorCode as ExportErrorCode, ErrorResponse as ExportErrorResponse
from drms_parameters import DRMSParams


# the drmsexport web application comprises five components:
#   + a set of HTML web pages that contain forms to collect information from the export-system user; the web pages use JavaScript,
#     HTML elements, and AJAX tools to create HTTP AJAX requests
#   + a browsers/network tools that send the form data contained within the AJAX requests to an HTTP server; the
#     browsers/network tools receive HTTP responses, updating the web pages displayed
#   + an HTTP server that acts as a reverse proxy: it receives HTTP requests from browsers/network tools on the internet, forwards
#     the contained data as uwsgi requests to an upstream WSGI server, receives uwsgi responses from the WSGI server, and finally sends
#     HTTP responses back to the originating browsers/tools
#   + a WSGI server that receives uwsgi requests from the reverse-proxy server, sends corresponding WSGI requests to the
#     Flask web-application entry point, and receives WSGI responses from the Flask web-application entry point, and sends uwsgi responses
#     back to the originating reverse-proxy server
#   + a Flask app that services WSGI requests it receives from the WSGI server, and sends WSGI responses back to the WSGI server

#########################
#    HTML web pages     #
#########################
# the drmsexport web application HTML and JavaScript files are in proj/export/webapps
#   + the static web pages are exportdata.html and export_request_form.html; they contain in-line JavaScript, as
#     well as references to JavaScript contained in separate files
#     - exportdata.html: this file contains JavaScript only; it includes a JavaScript script that contains a single string variable
#       that contains a text representation of the export_request_form.html
#     - export_request_form.html: this file contains the definitions of the HTML elements that compose the export web page
#   + the export JavaScript files are:
#     - export_request_form.htmlesc: this is a version of export_request_form.html that has been converted into a single JavaScript string
#       variable, where whitespace has been removed and characters have been percent-escaped if necessary; exportdata.html
#       'includes' this file as a JavaScript script
#     - export_email_verify.js: this file contains code that makes HTTP requests that access the email-registration system
#     - no_processing.js: this file contains a single array JavaScript variable that lists all DRMS data series for which
#       export-processing is prohibited
#     - processing.js: this file contains code that makes HTTP requests that cause export processing to occur during the export
#       process
#     - protocols.js: this file contains code that makes HTTP requests that gets image-protocol (export file type) parameters
#     - prototype.js: this is a JavaScript framework that provides additional functionality, such as the ability to make AJAX requests;
#       we probably are using version 1.7
#     - prototip.js: this is a JavaScript microframework that provides "tooltip" functionality - when the user clicks on a tooltip element
#       a text bubble appears
#     - cookies.js: this is a JavaScript microframework that provides cookie support; it requires prototype js
#     - prototip.css: css sheet for prototip microframework


#########################
# browser/network tool  #
#########################
# the main export-system web page is http://solarweb2.stanford.edu:8080/export; there are several endpoints:
#   + http://solarweb2.stanford.edu:8080/export/address-registration: this endpoint provides access to
#     services that check the registration status of an email address, and register a new email address; arguments:
#     - address: the email address to check on/register
#     - [ db-host ]: the INTERNAL/PRIVATE database server that hosts the registered export-system user address information
#     - [ db-name ]: the name of the database that contains email address and user informaton
#     - [ db-port ]: the port on the database host machine accepting connections
#     - [ db-user ]: the name of the database user account to use
#     - [ user-name ]: the full name of the export-system user
#     - [ user-snail ]: the physical address of the export-system user
#   + http://solarweb2.stanford.edu:8080/export/series-server: this endpoint provides access to
#     services that
#     - public_db_host: the EXTERNAL/PUBLIC database server that hosts the DRMS data-series data
#     - series_set: the set of DRMS data series for which information is to be obtained
#     - webserver: the webserver of this endpoint
#     - [ client-type ]: the securedrms client type (ssh, http)
#     - [ db-name ]: the name of the database that contains DRMS data-series information
#     - [ db-port ]: the port on the database host machine accepting connections
#     - [ db-user ]: the name of the database user account to use
#   + http://solarweb2.stanford.edu:8080/export/record-set: this endpoint provides access to
#     services that provide keyword, segment, and link informaton about DRMS record sets; arguments:
#     - specification: the DRMS record-set specification identifying the records for which information is to be obtained
#     - db-host: the database server that hosts the DRMS record-set data
#     - webserver: the webserver of this endpoint
#     - [ parse-only ]: if True, then parse record-set string only
#     - [ client-type ]: the securedrms client type (ssh, http)
#     - [ keywords ]: the list of keywords for which information is to be obtained
#     - [ segments ]: the list of segments for which information is to be obtained
#     - [ links ]: the list of links for which information is to be obtained
#     - [ db-name ]: the name of the database that contains DRMS record-set information
#     - [ db-port ]: the port on the database host machine accepting connections
#     - [ db-user ]: the name of the database user account to use
#   + http://solarweb2.stanford.edu:8080/export/series: this endpoint provides access to
#     services that provide informaton about DRMS data series; arguments:
#     - series: the DRMS series for which information is to be obtained
#     - db-host: the database server that hosts the DRMS data-series information
#     - webserver: the webserver of this endpoint
#     - [ client-type ]: the securedrms client type (ssh, http)
#     - [ db-name ]: the name of the database that contains DRMS data-series information
#     - [ db-port ]: the port on the database host machine accepting connections
#     - [ db-user ]: the name of the database user account to use
#   + http://solarweb2.stanford.edu:8080/export/new-premium-request: this endpoint provides access to
#     services that export DRMS data-series data; the full suite of export options is available; arguments:
#     - address: the email address registered for export
#     - db-host: the database server that hosts the DRMS data series
#     - webserver: the webserver of this endpoint
#     - arguments: the export-request arguments
#     - [ client-type ]: the securedrms client type (ssh, http)
#     - [ db-name ]: the name of the database that contains DRMS data-series information
#     - [ db-port ]: the port on the database host machine accepting connections
#     - [ db-user ]: the name of the database user account to use
#     - [ requestor ]: the full name of the export-system user
#   + http://solarweb2.stanford.edu:8080/export/new-mini-request: this endpoint provides access to
#     services that export DRMS data-series data; a reduced suite of export options is available
#     to allow for quicker payload delivery; arguments:
#     - address: the email address registered for export
#     - db-host: the database server that hosts the DRMS data series
#     - webserver: the webserver of this endpoint
#     - arguments: the export-request arguments
#     - [ client-type ]: the securedrms client type (ssh, http)
#     - [ db-name ]: the name of the database that contains DRMS data-series information
#     - [ db-port ]: the port on the database host machine accepting connections
#     - [ db-user ]: the name of the database user account to use
#     - [ requestor ]: the full name of the export-system user
#   + http://solarweb2.stanford.edu:8080/export/new-streamed-request: this endpoint provides access to
#     services that stream export DRMS data-series data; a reduced suite of export options is available
#     to allow for quicker payload delivery; arguments:
#     - address: the email address registered for export
#     - db-host: the database server that hosts the DRMS data series
#     - webserver: the webserver of this endpoint
#     - arguments: the export-request arguments
#     - [ client-type ]: the securedrms client type (ssh, http)
#     - [ db-name ]: the name of the database that contains DRMS data-series information
#     - [ db-port ]: the port on the database host machine accepting connections
#     - [ db-user ]: the name of the database user account to use
#     - [ requestor ]: the full name of the export-system user
#   + http://solarweb2.stanford.edu:8080/export/pending-request: this endpoint provides access to
#     services that check for the presence of pending requests; arguments:
#     - address: the email address registered for export
#     - db-host: the database server that hosts the DRMS data series
#     - webserver: the webserver of this endpoint
#     - [ db-name ]: the name of the database that contains pending-request information
#     - [ db-port ]: the port on the database host machine accepting connections
#     - [ db-user ]: the name of the database user account to use
#     - [ pending_requests_table ]: the database table of pending requests
#     - [ timeout ]: after this number of minutes have elapsed, requests are no longer considered pending
#   + http://solarweb2.stanford.edu:8080/export/pending-request-status: this endpoint provides access to
#     services that return the export status of a pending request; arguments:
#     - address: the email address registered for export
#     - db-host: the database server that hosts export-request information
#     - webserver: the webserver of this endpoint
#     - request-id: the export system request ID
#     - [ client-type ]: the securedrms client type (ssh, http)
#     - [ db-name ]: the name of the database that contains export-request information
#     - [ db-port ]: the port on the database host machine accepting connections
#     - [ db-user ]: the name of the database user account to use
#     - [ pending_requests_table ]: the database table of pending requests
#     - [ timeout ] after this number of minutes have elapsed, requests are no longer considered pending

#########################
# reverse-proxy (nginx) #
#########################
#   + to configure nginx to act as a proxy-server that sends uwsgi requests to the WSGI server, edit the default nginx
#     configuration file, /etc/nginx/nginx.conf:
#       - comment-out the existing `server` block; make sure the following `include` directive exists in the `html` block:
#         include /etc/nginx/conf.d/*.conf;
#       - create the directory /etc/nginx/conf.d if it does not already exist
#       - create /etc/nginx/conf.d/export.conf; this will contain a new server block for the drm_export app; the include
#         statement (above) will ensure that the drms_export server is defined; the content of export.conf should be
#         as follows:
#
#         server {
#             listen 8080;
#             charset utf-8;
#
#             location / {
#                 root /usr/share/nginx/html;
#                 try_files $uri $uri/index.html $uri.html;
#             }
#
#             location /export/ {
#                 # provided by nginx installation - relative to dir containing nginx.conf
#                 include uwsgi_params;
#
#                 # must match file specified in export.ini `socket` parameter; unix socket for speed
#                 # of communication between nginx and uWSGI running on the same host
#                 uwsgi_pass unix:///tmp/export.sock;
#
#                 # Redefine the header fields that NGINX sends to the upstream server
#                 proxy_set_header Host $host;
#                 proxy_set_header X-Real-IP $remote_addr;
#                 proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
#
#                 # Define the maximum file size on file uploads
#                 client_max_body_size 5M;
#             }
#
#             error_page 404 /404.html;
#             location = /404.html {
#             }
#
#             error_page 500 502 503 504 /50x.html;
#             location = /50x.html {
#             }
#         }
#   + the nginx `uwsgi_pass` directive ensures communication protocol between nginx and uWSGI is WSGI; the HTTP protocol
#     is also an option, with `proxy_pass`, but it is slower is not preferred

#########################
#  uWSGI (WSGI server)  #
#########################
# the WSGI server used by drmsexport is uWSGI
#   + upon receiving a uwsgi request from the reverse-proxy server, the uWSGI calls the Flask app's entry point, a WSGI
#     callable, which is exported from export_wsgi.py
#   + when calling the entry point, the WSGI server passes a dict of environment variables and a callable `start_response`,
#     encapsulated in a flask.Request instance; the Flask.Request contains the request arguments, such as an email
#     address or export-system request ID; the callable is used by the Flask app to return the application's response to
#     the WSGI server
#   + the WSGI server creates a uwsgi response from the arguments provided to the `start_response` callable, and sends the
#     response to the reverse-proxy server
#
# to configure the WSGI server, uWSGI, create a n .ini configuration file:
#   + each parameter of the configuration file has the format <parameter> = <value>:
#     [uwsgi]
#     # the location where uWSGI will search for the entry point (ensure the project directory exists)
#     chdir = /home/netdrms_production/export
#     module = export_wsgi:export
#
#     # run in master mode, spawning processes to handle requests
#     master = true
#     processes = 5
#
#     # listen for requests on a unix socket; ideal for configurations where the reverse-proxy server and uWSGI run on the
#     # same host
#     socket = /tmp/export.sock
#     chmod-socket = 664
#     vacuum = true

#########################
#   WSGI entry point    #
#########################
# the WSGI entry point is a Flask app callable that accepts two arguments: `environ` and `start_response`
#   + the Flask app is itself a callable (it has a __call__() method), so the entry point is simply a module,
#     export_wsgi.py, that imports the Flask app `export` from this script, export.py
#   + the WSGI server sends a request to the Flask app by calling the entry point:
#     export_wsgi.export(environ, start_response)
#     - export_wsgi.export(environ, start_response) calls export_wsgi.export.wsgi_app(environ, start_response)
#     - `environ` is a dict that contains the WSGI environment (e.g., type of HTTP request)
#     - `start_response` is a callable (a callback) that accepts a status code, a list of headers (each element
#        is a 2-tuple (<HTTP header name>, <HTTP header value>)), and an optional exception context; the
#        Flask app calls `start_response` to return a response to the WSGI server
#     - Flask encapsulates the receiving of requests from the WSGI server; it will convert the `environ` dict into a Flask.Request
#       and route the URL endpoint to the appropriate view function (the function that generates a response)
#     - Flask encapsulates the sending of responses to the WSGI server; the view function can return data in one of
#       many formats (iterable, dict, str, tuple, json text, html, flask.Response) and the results will be converted to a
#       Flask.Response; if a tuple is returned, then the view function can return (response, status), (response, headers),
#       or (response, headers, status) - status overrides the response status code, and headers is a list or dict

#########################
# Flask app (`export`)  #
#########################
# the Flask application is an instance of flask.Flask
#   + this script, export.py, creates the flask.Flask instance `export`; it provides access to "view functions"
#     (methods that generate responses to app requests) and "URL routes" (mappings from each URL endpoint)
#     to the view functions
#   + it creates the view functions by creating instances of flask_restful.Resource, one for each endpoint
#     - flask_restful.Resource is a subclass of flask.views.MethodView, which is what is called a "view class";
#       view classes offer one way to declare app view functions; flask_restful.Resource.dispatch_request() method
#       intercepts the request and then calls the appropriate method to handle the request (e.g., it calls the
#       `get` method to handle an HTTP GET request)
#     - the script uses flask_restful to create a RESTful web service; each endpoint names a web resource (object, noun)
#       that can be accessed by one or more of the HTTP protocol's request types (GET, POST, DELETE, PUT, PATCH)
#     - this RESTful interface is encapsulated in flask_restful.Resource subclasses, one for each web resource:
#
#       class PendingRequestResource(Resource):
#       '''
#       interface to an existing export request that is currently being processed by the export system
#       '''
#         ...
#        # view function to handle HTTP GET requests to the pending-request resource
#        def get(self):
#           ...
#   + it establishes the URL routing, registering the view class with the app, by creating a flask_restful.Api
#     instance and calling its add_resource() method
#     - the flask_restful.Api wraps the app
#     - the script calls flask_restful.Api.add_resource(), once for each flask_restful.Resource; for example, the
#       following registers the view class `AddressRegistrationResource` to be called when the WSGI server
#       sends a request to the Flask app:
#
#       export_api.add_resource(AddressRegistrationResource, '/export/address-registration')
#
#       flask_restful.Api.add_resource() registers the view class by calling Flask.flask.add_url_rule(); the actual
#       request-type-specific method to be called in the view class is not determined at this time; instead, the
#       view class's dispatch_request() will be called, and it is up to this method to call the request-type-specific
#       method (e.g., the `get` method)
#   + each view function returns a suitable value (iterable, dict, str, tuple, json text, html, flask.Response),
#     which the Flask app then converts to a flask.Flask.Response; in particular flask.jsonify() makes a flask.Response
#     from a dict and sets the flask.Reponse mimetype to `application/json`, which is necessary for the proper interpretation
#     by the original browser/network tool
#   + the flask.Flask.Response is then returned to the app's callable code, which then creates an iterable from
#     this Response and passes it to the `start_response` callable passed to the callable code
#   + the URL routing that leads to the execution of a view function works as follows:
#     - the app receives a WSGI request from the WSGI server, creating a global flask.Request instance
#     - the flask.Request instance contains a `method` property whose value is the type of HTTP request, like `GET`
#     - in response to receiving this request, the app consults the list of registered URL routes and then calls
#       the appropriate flask_restful.Resource.dispatch_request() method
#     - flask_restful.Resource.dispatch_request() examines flask.Request.method and calls the appropriate
#       view function - if the method is 'GET', then flask_restful.Resource.dispatch_request() calls the
#       flask_restful.Resource.get() view function

#############################
# launching the web service #
#############################
# launching the Flask app entails creating a project directory owned by the production user, copying
# the uWSGI configuration file and Flask-app files into place, and starting the export web service:
#   + create the project directory:
#     $ mkdir /home/netdrms_production/export
#   + copy the uWSGI configuration file, the python-module file that contains the Flask app definition, and
#     the file that contains app entry point to this project directory:
#     $ cd /opt/jsoc_drms/proj/export/webapps/drmsexport
#     $ cp export.ini /home/netdrms_production/export
#     $ cp *.py /home/netdrms_production/export
#   + start the web service

############################
# starting the web service #
############################
# starting the export web service entails starting the WSGI service (uWSGI), and starting the reverse-proxy
# service (nginx)
#   + to start uWSGI, run the uwsgi executable, providing export.ini as an argument:
#     uwsgi --ini=/home/netdrms_production/export/export.ini
#   + to start the nginx service:
#     sudo service nginx start

############################
# stopping the web service #
############################
# stopping the export web service entails stopping the WSGI service (uWSGI), and stopping the reverse-proxy
#   + to stop the nginx service, run the service executable:
#     $ sudo service nginx stop
#   + to stop the uWSGI service, run the kill command:
#     $ kill -9 <uwsgi pid>

# `export` is a Flask app;
export = Flask(__name__)
export_api = Api(export)
APP_DEBUG = True
export.logger.setLevel(LOGGING_LEVEL_DEBUG)

class ErrorCode(ExportErrorCode):
    # fetch error codes
    WEBARGS_PARSER = (1, 'failure parsing web arguments') # can only happen with jsoc_fetch op=exp_request call

@parser.error_handler
def handle_request_parser_error(error, req, schema, *, error_status_code, error_headers):
    print(f'webarg error message: {error.messages}')
    response_dict = ExportErrorResponse.generate_response(error_code=ErrorCode.WEBARGS_PARSER, error_message=error.messages).generate_serializable_dict()
    abort(428, **response_dict)

@export.errorhandler(428)
def generate_error_response(error):
    print(f'here XXX')

@export.before_request
def log_request_info():
    export.logger.debug(f'Headers: {request.headers}')
    export.logger.debug(f'Body: {request.get_data()}')

from check_address import CheckAddressAction
class AddressRegistrationResource(Resource):
    _arguments = { 'address' : fields.Str(required=True, validate=lambda a: a.find('@') >= 0), 'db_host' : fields.Str(required=False, data_key='db-host'), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'db_user' : fields.Str(required=False, data_key='db-user'), 'user_name' : fields.Str(required=False, data_key='user-name'), 'user_snail' : fields.Str(required=False, data_key='user-snail') }

    @use_kwargs(_arguments, location='querystring')
    def get(self, address, db_host=None, db_name=None, db_port=None, db_user=None, user_name=None, user_snail=None):
        # HTTP GET - get address registration information (if the address is registered)
        arguments = { 'address' : address, 'db_host' : db_host, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user, 'user_name' : user_name, 'user_snail' : user_snail }
        if APP_DEBUG:
            arguments['logging_level'] = 'debug'

        action = Action.action(action_type='check_address', args=arguments)
        return action().generate_serializable_dict()

    @use_kwargs(_arguments)
    def post(self, address, db_host=None, db_name=None, db_port=None, db_user=None, user_name=None, user_snail=None):
        # HTTP POST - register address (if it is not alrelady registered)
        arguments = { 'address' : address, 'db_host' : db_host, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user, 'user_name' : user_name, 'user_snail' : user_snail }
        if APP_DEBUG:
            arguments['logging_level'] = 'debug'

        action = Action.action(action_type='register_address', args=arguments)
        return action().generate_serializable_dict()

# called from public website only
from check_dbserver import DetermineDbServerAction
class ServerResource(Resource):
    _arguments = { 'public_db_host' : fields.Str(required=True, data_key='public-db-host'), 'series' : fields.List(fields.Str, required=True, validate=lambda a: DetermineDbServerAction.is_valid_series_set(a, None, urlparse(request.base_url).hostname, 'debug' if APP_DEBUG else None)), 'client_type' : fields.Str(required=False, data_key='client-type'), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'db_user' : fields.Str(required=False, data_key='db-user') }

    @use_kwargs(_arguments)
    def get(self, public_db_host, series, client_type=None, db_name=None, db_port=None, db_user=None):
        arguments = { 'public_db_host' : public_db_host, 'series' : series, 'drms_client_type' : client_type, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user }
        if APP_DEBUG:
            # no logging_level in determine_db_server action
            pass

        action = Action.action(action_type='determine_db_server', args=arguments)
        return action().generate_serializable_dict()

from get_record_info import GetRecordInfoAction
class RecordSetResource(Resource):
    _arguments = { 'specification' : fields.Str(required=True, validate=lambda a: GetRecordInfoAction.is_valid_specification(a, None, urlparse(request.base_url).hostname, 'debug' if APP_DEBUG else None)), 'db_host' : fields.Str(required=True, data_key='db-host'), 'parse_only' : fields.Bool(required=False, data_key='parse-only'), 'client_type' : fields.Str(required=False, data_key='client-type'), 'keywords' : fields.List(fields.Str, required=False), 'segments' : fields.List(fields.Str, required=False), 'links' : fields.List(fields.Str, required=False), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'db_user' : fields.Str(required=False, data_key='db-user') }

    _parse_response = None

    @use_kwargs(_arguments)
    def get(self, specification, db_host, parse_only=False, client_type=None, keywords=None, segments=None, links=None, db_name=None, db_port=None, db_user=None):
        if parse_only:
            if self._parse_response is None:
                arguments = { 'specification' : specification, 'db_host' : db_host, 'drms_client_type' : client_type, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user }
                if APP_DEBUG:
                    arguments['logging_level'] = 'debug'

                action = Action.action(action_type='parse_specification', args=arguments)
                self._parse_response = action()

            return self._parse_response.generate_serializable_dict()
        else:
            arguments = { 'specification' : specification, 'db_host' : db_host, 'webserver' : urlparse(request.base_url).hostname, 'drms_client_type' : client_type, 'keywords' : keywords, 'segments' : segments, 'links' : links, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user }
            if APP_DEBUG:
                arguments['logging_level'] = 'debug'

            action = Action.action(action_type='get_record_set_info', args=arguments)
            return action().generate_serializable_dict()

from get_series_info import GetSeriesInfoAction
class SeriesResource(Resource):
    _arguments = { 'series' : fields.Str(required=True, validate=lambda a: GetSeriesInfoAction.is_valid_series_set([ a ], None, urlparse(request.base_url).hostname, 'debug' if APP_DEBUG else None)), 'db_host' : fields.Str(required=True, data_key='db-host'), 'client_type' : fields.Str(required=False, data_key='client-type'), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'db_user' : fields.Str(required=False, data_key='db-user') }

    @use_kwargs(_arguments)
    def get(self, series, db_host, client_type=None, db_name=None, db_port=None, db_user=None):
        arguments = { 'series' : series, 'db_host' : db_host, 'webserver' : urlparse(request.base_url).hostname, 'drms_client_type' : client_type, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user }
        if APP_DEBUG:
            arguments['logging_level'] = 'debug'

        action = Action.action(action_type='get_series_info', args=arguments)
        return action().generate_serializable_dict()

from initiate_request import InitiateRequestAction, ErrorCode as IRErrorCode
class PremiumExportRequestResource(Resource):
    _arguments = { 'address' : fields.Str(required=True, validate=lambda a: a.find('@') >= 0), 'db_host' : fields.Str(required=True, data_key='db-host'), 'export_arguments' : fields.Str(required=True, validate=lambda a: InitiateRequestAction.is_valid_arguments(a, 'debug' if APP_DEBUG else None), data_key='export-arguments'), 'client_type' : fields.Str(required=False, data_key='client-type'), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'requestor' : fields.Str(required=False), 'db_user' : fields.Str(required=False, data_key='db-user') }

    @use_kwargs(_arguments)
    def post(self, address, db_host, export_arguments, client_type=None, db_name=None, db_port=None, requestor=None, db_user=None):
        arguments = { 'address' : address, 'db_host' : db_host, 'webserver' : urlparse(request.base_url).hostname, 'export_arguments' : export_arguments, 'drms_client_type' : client_type, 'db_name' : db_name, 'db_port' : db_port, 'requestor' : requestor, 'db_user' : db_user }
        if APP_DEBUG:
            arguments['logging_level'] = 'debug'

        action = Action.action(action_type='start_premium_export', args=arguments)
        return action().generate_serializable_dict()

class MiniExportRequestResource(Resource):
    _arguments = { 'address' : fields.Str(required=True, validate=lambda a: a.find('@') >= 0), 'db_host' : fields.Str(required=True, data_key='db-host'), 'export_arguments' : fields.Str(required=True, validate=lambda a: InitiateRequestAction.is_valid_arguments(a, 'debug' if APP_DEBUG else None), data_key='export-arguments'), 'client_type' : fields.Str(required=False, data_key='client-type'), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'requestor' : fields.Str(required=False), 'db_user' : fields.Str(required=False, data_key='db-user') }

    @use_kwargs(_arguments)
    def post(self, address, db_host, export_arguments, client_type=None, db_name=None, db_port=None, requestor=None, db_user=None):
        arguments = { 'address' : address, 'db_host' : db_host, 'webserver' : urlparse(request.base_url).hostname, 'export_arguments' : export_arguments, 'drms_client_type' : client_type, 'db_name' : db_name, 'db_port' : db_port, 'requestor' : requestor, 'db_user' : db_user }
        if APP_DEBUG:
            arguments['logging_level'] = 'debug'

        action = Action.action(action_type='start_mini_export', args=arguments)
        return action().generate_serializable_dict()

class StreamedExportRequestResource(Resource):
    _arguments = { 'address' : fields.Str(required=True, validate=lambda a: a.find('@') >= 0), 'db_host' : fields.Str(required=True, data_key='db-host'), 'export_arguments' : fields.Str(required=True, validate=lambda a: InitiateRequestAction.is_valid_arguments(a, 'debug' if APP_DEBUG else None), data_key='export-arguments'), 'client_type' : fields.Str(required=False, data_key='client-type'), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'requestor' : fields.Str(required=False), 'db_user' : fields.Str(required=False, data_key='db-user') }

    @use_kwargs(_arguments)
    def post(self, address, db_host, export_arguments, client_type=None, db_name=None, db_port=None, requestor=None, db_user=None):
        arguments = { 'address' : address, 'db_host' : db_host, 'webserver' : urlparse(request.base_url).hostname, 'export_arguments' : export_arguments, 'drms_client_type' : client_type, 'db_name' : db_name, 'db_port' : db_port, 'requestor' : requestor, 'db_user' : db_user }
        if APP_DEBUG:
            arguments['logging_level'] = 'debug'

        # creates a SecureExportRequest object, but does not start the child process;
        action = Action.action(action_type='start_streamed_export', args=arguments)

        # starts the child process and creates a generator to return data from the output content stream
        # response, when successful, will not be sent back to browser, but it can be used here to check for errors
        action_response = action()

        if isinstance(action_response, IRErrorCode):
            # return an error response using the `action_response` description
            return make_response(action_response.attributes.drms_export_status_description, 500)

        # call the generator's first iteration; this reads the download file name from the
        # child process' output stream and stores it in the destination object
        headers = action.generate_headers()

        # `action.generator` is the generator object (whose first iteration has already been invoked - to
        # get the download file name from the child process' output stream); the remaining iterations
        # will return all the download file data, one block at a time
        return export.response_class(action.generator, content_type='application/octet-stream', headers=headers)

from manage_request import PendingRequestAction
class PendingRequestResource(Resource):
    _arguments = { 'address' : fields.Str(required=True, validate=lambda a: a.find('@') >= 0), 'db_host' : fields.Str(required=True, data_key='db-host'), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'db_user' : fields.Str(required=False, data_key='db-user'), 'pending_requests_table' : fields.Str(required=False, data_key='pending-requests-table'), 'timeout' : fields.Int(required=False) }

    @use_kwargs(_arguments)
    def get(self, address, db_host, db_name=None, db_port=None, db_user=None, pending_requests_table=None, timeout=None):
        # HTTP GET - check for the existence of a pending request
        arguments = { 'address' : address, 'db_host' : db_host, 'webserver' : urlparse(request.base_url).hostname, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user, 'pending_requests_table' : pending_requests_table, 'timeout' : timeout }
        if APP_DEBUG:
            arguments['logging_level'] = 'debug'

        action = Action.action(action_type='check_pending_request', args=arguments)
        return action().generate_serializable_dict()

    @use_kwargs(_arguments)
    def post(self, address, db_host, db_name=None, db_port=None, db_user=None, pending_requests_table=None, timeout=None):
        # HTTP POST - cancel an export request
        arguments = { 'address' : address, 'db_host' : db_host, 'webserver' : urlparse(request.base_url).hostname, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user, 'pending_requests_table' : pending_requests_table, 'timeout' : timeout }
        if APP_DEBUG:
            arguments['logging_level'] = 'debug'

        action = Action.action(action_type='cancel_pending_request', args=arguments)
        return action().generate_serializable_dict()

class PendingRequestStatusResource(Resource):
    _arguments = { 'address' : fields.Str(required=True, validate=lambda a: a.find('@') >= 0), 'db_host' : fields.Str(required=True, data_key='db-host'), 'request_id' : fields.Str(required=True, data_key='request-id', validate=lambda a: PendingRequestAction.is_valid_request_id(a)), 'client_type' : fields.Str(required=False, data_key='client-type'), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'db_user' : fields.Str(required=False, data_key='db-user'), 'pending_requests_table' : fields.Str(required=False, data_key='pending-requests-table'), 'timeout' : fields.Int(required=False) }

    @use_kwargs(_arguments)
    def get(self, address, db_host, request_id, client_type=None, db_name=None, db_port=None, db_user=None, pending_requests_table=None, timeout=None):
        # HTTP GET - get the export status of a request
        arguments = { 'address' : address, 'db_host' : db_host, 'webserver' : urlparse(request.base_url).hostname, 'request_id' : request_id, 'drms_client_type' : client_type, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user, 'pending_requests_table' : pending_requests_table, 'timeout' : timeout }
        if APP_DEBUG:
            arguments['logging_level'] = 'debug'

        action = Action.action(action_type='get_export_status', args=arguments)
        return action().generate_serializable_dict()

class RootResource(Resource):
    def get(self):
        return "<h1 style='color:blue'>Hello There!</h1>"

export_api.add_resource(AddressRegistrationResource, '/export/address-registration')
export_api.add_resource(ServerResource, '/export/series-server')
export_api.add_resource(RecordSetResource, '/export/record-set')
export_api.add_resource(SeriesResource, '/export/series')

export_api.add_resource(PremiumExportRequestResource, '/export/new-premium-request')
export_api.add_resource(MiniExportRequestResource, '/export/new-mini-request')
export_api.add_resource(StreamedExportRequestResource, '/export/new-streamed-request')

export_api.add_resource(PendingRequestResource, '/export/pending-request')
export_api.add_resource(PendingRequestStatusResource, '/export/pending-request-status')
export_api.add_resource(RootResource, '/export')

if __name__ == '__main__':
    drms_params = DRMSParams()
    app.run(host=drms_params.EXPORT_WEB_SERVER, port=drms_parms.EXPORT_WEB_SERVER_PORT, debug=True)
else:
    # export was imported by export_wsgi
    pass
