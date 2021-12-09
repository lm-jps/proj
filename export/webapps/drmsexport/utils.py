#!/usr/bin/env python3

from action import Action
from drms_export import ErrorResponse

__all__ = [ 'extract_program_and_module_args', 'get_db_host' ]

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


# determines which db host - the private one, or a public one - is needed to serve export-client requests; the
# determination is based upon the provided webserver (which has a `public` property), series, and
# specification arguments
# `webserver` can be None; if so, then `webserver` and `db_host` must be public
# `series` is a comma-separated list of DRMS data series
# `db_host` the database host selected to support requests from `webserver` (public webserver --> public db host)
def get_db_host(*, webserver, series, private_db_host, db_host, db_port, db_name, db_user, exc, log):
    from check_dbserver import DetermineDbServerAction

    # parse specification to obtain series so we can check for pass-through series (relevant only if the user is on a public webserver);
    # parsing checks syntax, it does not check for the existence of series in the specification, so either a public or private drms client can be used
    resolved_db_host = None

    if webserver is None or webserver.public:
        # db_host must be a public db server
        if webserver is not None:
            log.write_debug([ f'[ get_db_host ] public webserver `{webserver.host}` initiating series-information request' ])

        # need to determine if pass-through series have been specified; if so, use securedrms client that uses private db
        if series is not None:
            log.write_debug([ f'[ get_db_host ] determining DB server suitable for requested data from series `{", ".join(series)}`' ])

            action_type = 'determine_db_server'
            action_args = { 'log' : log, 'public_db_host' : db_host, 'series' : series }
            action = Action.action(action_type=action_type, args=action_args)
            response = action()

            if isinstance(response, ErrorResponse):
                log.write_error([ f'[ get_db_host ] {response.attributes.error_message}'])
                raise exc(error_message=f'failure calling `{action_type}` action')

            # could be either public or private db host
            resolved_db_host = response.attributes.server
    else:
        log.write_debug([ f'[ get_db_host ] private webserver `{webserver.host}` initiating series-information request' ])
        resolved_db_host = private_db_host

    return resolved_db_host
