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

class RequestResource(Resource):
    user_args = { 'address' : fields.Str(required=True, validate=lambda a: a.find('@') >= 0 ) }

    @use_kwargs(user_args)
    def post(self, address):
        # HTTP POST - cancel an export request
        action = Action.action(action_type='cancel_pending_request', args={ "address" : address })
        return jsonify(action())

    @use_kwargs(user_args)
    def get(self, address):
        # HTTP GET - check for the existence of a pending request
        action = Action.action(action_type='check_pending_request', args={ "address" : address })
        return jsonify(action())

class RootResource(Resource):
    def get(self):
        return "<h1 style='color:blue'>Hello There!</h1>"

export_api.add_resource(RequestResource, '/export/request')
export_api.add_resource(RootResource, '/export')

if __name__ == '__main__':
    sp.append(op.join(op.dirname(op.realpath(__file__)), '../../../include'))
    from drmsparams import DRMSParams

    drms_params = DRMSParams()
    app.run(host=drms_params.EXPORT_WEB_SERVER, port=drms_parms.EXPORT_WEB_SERVER_PORT, debug=True)
else:
    # export was imported by export_wsgi
    pass
