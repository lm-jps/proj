#!/usr/bin/env python3

from sys import path as sp
from os import path as op

from flask import Flask, request, jsonify, make_response
from flask_restful import Api, Resource
from webargs import fields, validate, flaskparser, ValidationError
from webargs.flaskparser import use_args, use_kwargs, parser, abort

# __file__ does not exist when run in an interpreter
# sp.append('/opt/netdrms/include')
# sp.append('/opt/netdrms/base/export/scripts')

from action import Action
from drms_export import ErrorResponse

# `export` is a flask app;

# the drms_export web application comprises X components:
#   + a set of HTML web pages that contain forms to collect information from the export-system user; the web pages use JavaScript,
#     HTML elements, and AJAX tools to create HTTP requests which are then sent by the browser to an HTTP server
#   + a browser/network tool that sends an HTTP server HTTP requests, and receives HTTP responses from the server
#   + an HTTP server that acts as a reverse proxy: it receives HTTP requests from browsers/network tools on the internet, forwards
#     the requests as uwsgi requests to an upstream WSGI server, receives uwsgi responses from the WSGI server, and finally sends
#     HTML responses back to the originating browsers/tools
#   + a WSGI server that receives uwsgi requests from the reverse-proxy server, sends corresponding WSGI requests to the
#     Flask web application, receives WSGI responses from the Flask web application, and sends uwsgi responses
#     back to the originating reverse-proxy server
#   + a Flask app that services WSGI requests it receives from the WSGI server, and sends WSGI responses back to the WSGI server

#########################
# reverse-proxy (NGINX) #
#########################


#########################
#  uWSGI (WSGI server)  #
#########################
# the WSGI server used by drms_export is uWSGI
#   + upon receiving a uwsgi request from the reverse-proxy server, the uWSGI calls the Flask app's entry point, a WSGI
#     callable, which is exported from export_wsgi.py
#   + when calling the entry point, the WSGI server passes a dict of environment variables and a callable `start_response`,
#     encapsulated in a flask.Request instance; the Flask.Request contains the request arguments, such as an email
#     address or export-system request ID; the callable is used by the Flask app to return the application's response to
#     the WSGI server
#   + the WSGI server creates a uwsgi response from the arguments provided to the `start_response` callable, and sends the
#     response to the reverse-proxy server

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


for the URL http://solarweb2.stanford.edu:8080/export/request, then
#     get() method of the resource's flask_restful.Resource object is invoked



export = Flask(__name__)
export_api = Api(export)

class AddressRegistrationResource(Resource):
    _arguments = { 'address' : fields.Str(required=True, validate=lambda a: a.find('@') >= 0), 'db_host' : fields.Str(required=False, data_key='db-host'), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'db_user' : fields.Str(required=False, data_key='db-user'), 'user_name' : fields.Str(required=False, data_key='user-name'), 'user_snail' : fields.Str(required=False, data_key='user-snail') }

    @use_kwargs(_arguments)
    def get(self, address, db_host, db_name, db_port, db_user, user_name, user_snail):
        # HTTP GET - get address registration information (if the address is registered)
        action = Action.action(action_type='check_address', args={ 'address' : address, 'db_host' : db_host, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user, 'user_name' : user_name, 'user_snail' : user_snail })
        response_dict = action().generate_dict()
        return jsonify(**response_dict)

    @use_kwargs(_arguments)
    def post(self, address,  db_host, db_name, db_port, db_user, user_name, user_snail):
        # HTTP POST - register address (if it is not alrelady registered)
        action = Action.action(action_type='register_address', args={ 'address' : address, 'db_host' : db_host, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user, 'user_name' : user_name, 'user_snail' : user_snail })
        response_dict = action().generate_dict()
        return jsonify(**response_dict)

# called from public website only
class ServerResource(Resource):
    _arguments = { 'public_db_host' : fields.Str(required=True, data_key='public-db-host'), 'series_set' : fields.Str(required=True, validate=lambda a: ServerResource._is_valid_series_set(a), data_key='series-set'), 'webserver' : fields.Str(required=True), 'client_type' : fields.Str(required=False, data_key='client-type'), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'db_user' : fields.Str(required=False, data_key='db-user') }

    @use_kwargs(_arguments)
    def get(self, public_db_host, series_set, webserver, client_type,  db_name, db_port, db_user):
        action = Action.action(action_type='determine_db_server', args={ 'public_db_host' : public_db_host, 'series' : series_json, 'drms_client_type' : client_type, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user })
        response_dict = action().generate_dict()
        return jsonify(**response_dict)

    @classmethod
    def _is_valid_series_set(cls, series_set):
        is_valid = None
        try:
            # the only check we can realistically perform is to check the json structure
            is_valid = True if type(json_loads(series_set)['series']) == list else False
        except:
            is_valid = False

        return is_valid

class RecordSetResource(Resource):
    _arguments = { 'parse_only' : fields.Bool(required=False, data_key='parse-only'), 'specification' : fields.Str(required=True, validate=lambda a: RecordSetResource._is_valid_specification(a)), 'db_host' : fields.Str(required=True, data_key='db-host'), 'webserver' : fields.Str(required=True), 'client_type' : fields.Str(required=False, data_key='client-type'), 'keywords' : fields.Str(required=False), 'segments' : fields.Str(required=False), 'links' : fields.Str(required=False), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'db_user' : fields.Str(required=False, data_key='db-user') }

    _parse_response = None

    @use_kwargs(_arguments)
    def get(self, parse_only, specification, db_host, webserver, client_type, keywords, segments, links, db_name, db_port, db_user):
        if parse_only:
            if self._parse_response is None:
                action = Action.action(action_type='parse_specification', args={ 'specification' : specification, 'db_host' : db_host, 'drms_client_type' : client_type, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user })
                self._parse_response = action()

            response_dict = self._parse_response.generate_dict()
            return jsonify(**response_dict)
        else:
            action = Action.action(action_type='get_record_set_info', args={ 'specification' : specification, 'db_host' : db_host, 'webserver' : webserver, 'drms_client_type' : client_type, 'keywords' : keywords, 'segments' : segments, 'links' : links, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user })
            response_dict = action().generate_dict()
            return jsonify(**response_dict)

    @classmethod
    def _is_valid_specification(cls, specifcation):
        is_valid = None
        try:
            # parse specification
            cls._parse_response = None
            action = Action.action(action_type='parse_specification', args={ 'specification' : specification })
            response = action()
            cls._parse_response = response
            is_valid = False if isinstance(response, ErrorResponse) else True
        except:
            is_valid = False

        return is_valid

class SeriesResource(Resource):
    _arguments = { 'series' : fields.Str(required=True, validate=lambda a: SeriesResource._is_valid_series(a)), 'db_host' : fields.Str(required=True, data_key='db-host'), 'webserver' : fields.Str(required=True), 'client_type' : fields.Str(required=False, data_key='client-type'), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'db_user' : fields.Str(required=False, data_key='db-user') }

    @use_kwargs(_arguments)
    def get(self, series, db_host, webserver, client_type, db_name, db_port, db_user):
        action = Action.action(action_type='get_series_info', args={ 'series' : series, 'db_host' : db_host, 'webserver' : webserver, 'drms_client_type' : client_type, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user })
        response_dict = action().generate_dict()
        return jsonify(**response_dict)

    @classmethod
    def _is_valid_series(cls, series):
        is_valid = None
        try:
            # parse specification
            action = Action.action(action_type='parse_specification', args={ 'specification' : series })
            response = action()
            if isinstance(response, ErrorResponse):
                is_valid = False
            else:
                if len(response.subsets) != 1:
                    # expected a single series
                    is_valid = False
                elif response.subsets[0].filter != None or response.subsets[0].segments != None:
                    # expect a series, not a record set
                    is_valid = False
        except:
            is_valid = False

        return is_valid

class ExportRequestResource(Resource):
    @classmethod
    def _is_valid_arguments(cls, arguments_json):
        is_valid = None
        try:
            json_loads(arguments_json)
            is_valid = True
        except:
            is_valid = False

        return is_valid

class PremiumExportRequestResource(ExportRequestResource):
    _arguments = { 'address' : fields.Str(required=True, validate=lambda a: a.find('@') >= 0), 'db_host' : fields.Str(required=True, data_key='db-host'), 'webserver' : fields.Str(required=True), 'arguments' : fields.Str(required=True, validate=lambda a: InitiateRequestResourse._is_valid_arguments(a)), 'client_type' : fields.Str(required=False, data_key='client-type'), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'requestor' : fields.Str(required=False), 'db_user' : fields.Str(required=False, data_key='db-user') }

    @use_kwargs(_arguments)
    def post(self, address, db_host, webserver, arguments, client_type, db_name, db_port, requestor, db_user):
        action = Action.action(action_type='start_premium_export', args={ 'address' : address, 'db_host' : db_host, 'webserver' : webserver, 'arguments' : arguments, 'drms_client_type' : client_type, 'db_name' : db_name, 'db_port' : db_port, 'requestor' : requestor, 'db_user' : db_user })
        response_dict = action().generate_dict()
        return jsonify(**response_dict)

class MiniExportRequestResource(ExportRequestResource):
    _arguments = { 'address' : fields.Str(required=True, validate=lambda a: a.find('@') >= 0), 'db_host' : fields.Str(required=True, data_key='db-host'), 'webserver' : fields.Str(required=True), 'arguments' : fields.Str(required=True, validate=lambda a: InitiateRequestResourse._is_valid_arguments(a)), 'client_type' : fields.Str(required=False, data_key='client-type'), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'requestor' : fields.Str(required=False), 'db_user' : fields.Str(required=False, data_key='db-user') }

    @use_kwargs(_arguments)
    def post(self, address, db_host, webserver, arguments, client_type, db_name, db_port, requestor, db_user):
        action = Action.action(action_type='start_mini_export', args={ 'address' : address, 'db_host' : db_host, 'webserver' : webserver, 'arguments' : arguments, 'drms_client_type' : client_type, 'db_name' : db_name, 'db_port' : db_port, 'requestor' : requestor, 'db_user' : db_user })
        response_dict = action().generate_dict()
        return jsonify(**response_dict)

class StreamedExportRequestResource(ExportRequestResource):
    _arguments = { 'address' : fields.Str(required=True, validate=lambda a: a.find('@') >= 0), 'db_host' : fields.Str(required=True, data_key='db-host'), 'webserver' : fields.Str(required=True), 'arguments' : fields.Str(required=True, validate=lambda a: InitiateRequestResourse._is_valid_arguments(a)), 'client_type' : fields.Str(required=False, data_key='client-type'), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'requestor' : fields.Str(required=False), 'db_user' : fields.Str(required=False, data_key='db-user') }

    @use_kwargs(_arguments)
    def post(self, address, db_host, webserver, arguments, client_type, db_name, db_port, requestor, db_user):
        action = Action.action(action_type='start_streamed_export', args={ 'address' : address, 'db_host' : db_host, 'webserver' : webserver, 'arguments' : arguments, 'drms_client_type' : client_type, 'db_name' : db_name, 'db_port' : db_port, 'requestor' : requestor, 'db_user' : db_user })
        # the action will dump exported-file data to stdout
        response_dict = action().generate_dict()
        return jsonify(**response_dict)

class PendingRequestResource(Resource):
    _arguments = { 'address' : fields.Str(required=True, validate=lambda a: a.find('@') >= 0), 'db_host' : fields.Str(required=True, data_key='db-host'), 'webserver' : fields.Str(required=True), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'db_user' : fields.Str(required=False, data_key='db-user'), 'pending_requests_table' : fields.Str(required=False, data_key='pending-requests-table'), 'timeout' : fields.Int(required=False) }

    @use_kwargs(_arguments)
    def get(self, address, db_host, webserver, db_name, db_port, db_user, pending_requests_table, timeout):
        # HTTP GET - check for the existence of a pending request
        action = Action.action(action_type='check_pending_request', args={ 'address' : address, 'db_host' : db_host, 'webserver' : webserver, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user, 'pending_requests_table' : pending_requests_table, 'timeout' : timeout })
        response_dict = action().generate_dict()
        return jsonify(**response_dict)

    @use_kwargs(_arguments)
    def post(self, address, db_host, webserver, db_name, db_port, db_user, pending_requests_table, timeout):
        # HTTP POST - cancel an export request
        action = Action.action(action_type='cancel_pending_request', args={ 'address' : address, 'db_host' : db_host, 'webserver' : webserver, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user, 'pending_requests_table' : pending_requests_table, 'timeout' : timeout })
        response_dict = action().generate_dict()
        return jsonify(**response_dict)

class PendingRequestStatusResource(Resource):
    _arguments = { 'address' : fields.Str(required=True, validate=lambda a: a.find('@') >= 0), 'db_host' : fields.Str(required=True, data_key='db-host'), 'webserver' : fields.Str(required=True), 'request_id' : fields.Str(required=True, data_key='request-id', validate=lambda a: PendingRequestStatusResource.is_valid_request_id(a)), 'client_type' : fields.Str(required=False, data_key='client-type'), 'db_name' : fields.Str(required=False, data_key='db-name'), 'db_port' : fields.Int(required=False, data_key='db-port'), 'db_user' : fields.Str(required=False, data_key='db-user'), 'pending_requests_table' : fields.Str(required=False, data_key='pending-requests-table'), 'timeout' : fields.Int(required=False) }

    @use_kwargs(_arguments)
    def get(self, address, db_host, webserver, request_id, client_type, db_name, db_port, db_user, pending_requests_table, timeout):
        # HTTP GET - get the export status of a request
        action = Action.action(action_type='get_export_status', args={ 'address' : address, 'db_host' : db_host, 'webserver' : webserver, 'request_id' : request_id, 'drms_client_type' : client_type, 'db_name' : db_name, 'db_port' : db_port, 'db_user' : db_user, 'pending_requests_table' : pending_requests_table, 'timeout' : timeout  })
        response_dict = action().generate_dict()
        return jsonify(**response_dict)

    @classmethod
    def is_valid_request_id(self, address):
        reg_ex = PendingRequestAction.get_reg_ex()
        return reg_ex.match(address) is not None

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
    sp.append(op.join(op.dirname(op.realpath(__file__)), '../../../include'))
    from drmsparams import DRMSParams

    drms_params = DRMSParams()
    app.run(host=drms_params.EXPORT_WEB_SERVER, port=drms_parms.EXPORT_WEB_SERVER_PORT, debug=True)
else:
    # export was imported by export_wsgi
    pass
