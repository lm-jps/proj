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

# `export` is a flask app; upon receiving an HTTP or WSGI request, the uWSGI web server calls export.run() by calling the 'callable' (the object with run() defined) in the WSGI entry point (export_wsgi.export); the WSGI server is configured by specifying the flask app `export` (defined here), imported into the entry point - the entry point contains from export import export
# to export.run(), the uWSGI web server passes a dict of environment variables, which contain the request arguments (such as an email address), and a callback method;
# export.run() then calls the appropriate flask_restful.Resource method to handle the request to the specified web resource; for example, if an HTTP GET request is received for the URL http://solarweb2.stanford.edu:8080/export/request, then get() method of the resource's flask_restful.Resource object is invoked;
# a flask_restful.Api object, `export_api`, is defined to establish this mapping from a URL resource to a flask_restful.Resource object
# each flask_restful.Resource method returns a response that export.run() then converts into a flask.Response object;
# export.run() then calls the callback method uWSGI server provided, providing the request status, and a list of tuples where each tuple is an (HTTP header, value)
# export.run() returns an iterable formed from the flask.Response that it received from the flask_restful.Resource method
# the uWSGI server then creates a response (in the agreed format, such as HTTP) from the iterable and sends it to the caller (curl, or an HTTP web server, for example)
export = Flask(__name__)

# `export_api` is a flask_restful.Api object that establishes API access to the wrapped flask app `export`; each web resource (i.e., a URL, like http://solarweb2.stanford.edu:8080/export/request) has an associated flask_restful.Resource object that handles all requests that resource; the calls to export_api.add_resource() define these associations, and in doing so, define the complete API;
# the various methods of each flask_restful.Resource object handle each type of request to the associated web resource; for example, if an HTTP GET request is received for the URL http://solarweb2.stanford.edu:8080/export/request, then the get() method of the associated flask_restful.Resource is invoked
# each flask_restful.Resource method returns to export.run() a response that can take one of several forms; a dict is acceptable, in which case export.run() will create a flask.Response object that contains JSON; a flask.Response object is also an acceptable response - the flask.jsonify() method is one way to create such an object; in particular, flask.jsonify() sets the flask.Response mimetype to `application/json`; it appears that Flask.run() can also accept list, or even HTML responses
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
