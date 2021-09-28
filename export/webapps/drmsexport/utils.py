#!/usr/bin/env python3

from action import Action
from drms_export import securedrms

__all__ = [ 'extract_program_and_module_args', 'create_drms_client' ]

def extract_program_and_module_args(*, is_program, **kwargs):
    program_args = None
    module_args = None

    if is_program:
        program_args = []
        for key, val in kwargs.items():
            if val is not None:
                if key == 'options':
                    for option, option_val in val.items():
                        if option_val is not None:
                            program_args.append(f'--{option}={option_val}')
                else:
                    program_args.append(f'{key}={val}')
    else:
        # a loaded module
        module_args = {}
        for key, val in kwargs.items():
            if val is not None:
                if key == 'options':
                    for option, option_val in val.items():
                        if option_val is not None:
                            module_args[option] = option_val
                else:
                    module_args[key] = val

    return (program_args, module_args)


# determines which type of drms client - a private one, or a public one - is needed to serve drms client requests, and returns the
# needed client; the determination is based upon the provided webserver (which has a `public` property), drms_client, series, and
# specification arguments
# `webserver` can be None; if so, then db_host must be public
# `series` is a comma-separated list of DRMS data series
def create_drms_client(*, webserver, address=None, series, specification, drms_client_type, drms_client=None, public_drms_client_server, private_drms_client_server, private_db_host, db_host, db_port, db_name, db_user, debug=False, log):
    from check_dbserver import StatusCode as CdbStatusCode

    # parse specification to obtain series so we can check for pass-through series (relevant only if the user is on a public webserver);
    # parsing checks syntax, it does not check for the existence of series in the specification, so either a public or private drms client can be used
    drms_client = None
    factory = None
    public_drms_client = None
    private_drms_client = None
    series_resolved = None

    if webserver is None or webserver.public:
        # db_host must be a public db server; use external drms client
        if webserver is not None:
            log.write_debug([ f'[ create_drms_client ] public webserver {webserver.host} using DRMS client' ])

        if drms_client is None:
            log.write_debug([ f'[ create_drms_client ] no securedrms client provided; creating public one' ])
            factory = securedrms.SecureClientFactory(debug=debug, email=address)
            use_ssh = True if drms_client_type == 'ssh' else False
            connection_info = { 'dbhost' : db_host, 'dbport' : db_port, 'dbname' : db_name, 'dbuser' : db_user }
            public_drms_client = factory.create_client(server=public_drms_client_server, use_ssh=use_ssh, use_internal=False, connection_info=connection_info)
        else:
            # public client (since the webserver is public)
            log.write_debug([ f'[ create_drms_client ] public securedrms client provided' ])
            public_drms_client = drms_client

        # need to determine if pass-through series have been specified; if so, use securedrms client that uses private db
        if series is not None:
            series_resolved = [ series ]
        elif specification is not None:
            log.write_debug([ f'[ create_drms_client ] parsing record-set specification `{specification}`' ])
            response_dict = public_drms_client.parse_spec(specification)

            if response_dict['errMsg'] is None:
                series_resolved = []

                # `subsets` exists if status == PsStatusCode.SUCCESS
                subsets = response_dict['subsets']

                for subset in subsets:
                    series_resolved.append(subset['seriesname'])

        if series_resolved is not None:
            log.write_debug([ f'[ create_drms_client ] determining DB server suitable for requested data from series `{", ".join(series_resolved)}`' ])

            action_type = 'determine_db_server'
            action_args = { 'public_db_host' : db_host, 'series' : series_resolved, 'drms_client' : public_drms_client }
            action = Action.action(action_type=action_type, args=action_args)
            response = action()

            if response.attributes.drms_export_status_code != CdbStatusCode.SUCCESS:
                log.write_error([ f'[ create_drms_client ] failure calling `{action_type}` action; status: `{response.attributes.drms_export_status_code.description()}`' ])
            elif response.attributes.server is None:
                log.write_error([ f'[ create_drms_client cannot service any series in `{", ".join(series_resolved)}`' ])
            else:
                db_host = response.attributes.server
    else:
        log.write_debug([ f'[ create_drms_client ] private webserver `{webserver.host}` initiating series-information request' ])
        private_drms_client = drms_client # could be None

    if db_host is not None:
        # now we know whether we should be using a public or private drms client
        use_public_db_host = True if (webserver is None or webserver.public) and db_host != private_db_host else False
        if use_public_db_host:
            log.write_debug([ f'[ create_drms_client ] public db host will be used to service request' ])
            if public_drms_client is not None:
                drms_client = public_drms_client
            else:
                # we must have had a public client - we used it to determine the db server
                log.write_error([ f'[ create_drms_client ] missing public drms client'])
        else:
            log.write_debug([ f'[ create_drms_client ] private db host must be used to service request' ])
            # private drms client needed; if a public one was passed in, then private_drms_client is None
            if private_drms_client is not None:
                log.write_debug([ f'[ create_drms_client ] no securedrms client provided; creating private one' ])
                drms_client = private_drms_client
            else:
                try:
                    if factory is None:
                        factory = securedrms.SecureClientFactory(debug=debug, email=arguments.address)
                    use_ssh = True if drms_client_type == 'ssh' else False
                    connection_info = { 'dbhost' : private_db_host, 'dbport' : db_port, 'dbname' : db_name, 'dbuser' : db_user }
                    log.write_debug([ f'[ create_drms_client ] creating private securedrms client' ])
                    private_drms_client = factory.create_client(server=private_drms_client_server, use_ssh=use_ssh, use_internal=True, connection_info=connection_info)
                    drms_client = private_drms_client
                except Exception as exc:
                    log.write_error([ f'[ create_drms_client ] {str(exc)}'])

    return drms_client
