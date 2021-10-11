function email_initVars(element_address, element_requestor)
{
    // 'valid' indicates if the string currently in the field is valid or not; this gets set only when the field value
    //     changes
    // 'addresses' is an object of validated email address (the properties are email addresses, and the values are booleans that
    //     indicate a registered email address or not, or 'pending' when checkAddress is in the process of checking
    //     the address registration status):
    //       'addresses' : { 'fred@bedrock.com' : true, 'barney@bedrock.com' : false, 'mrslate@bedrock.com' : 'pending' }
    element_address.store({ 'valid' : null, 'addresses' : {} });
    element_requestor.store({ 'valid' : null });
}

// Global constants
var JSOC_CHECK_EMAIL = "checkAddress.sh";
var EXPORT_PATH = '/export';
var REGISTRATION_RESOURCE = 'address-registration';
var REGISTRATION_PATH = [ EXPORT_PATH, REGISTRATION_RESOURCE ].join('/')
var MAX_NOTIFY_TIMER = 10000; // ms, 10 seconds
var MAX_REGISTRATION_TIME = 180; // seconds

function email_getargs(element_address, element_requestor)
{
    var requestor = null;
    var address = null;
    var addresses = null;

    // these cookies return undefined if they are not set
    requestor = Cookie.getData("user");
    address = Cookie.getData("notify");

    if (requestor === null || requestor === undefined || requestor.trim() == "undefined")
    {
        element_requestor.value = '';
    }
    else
    {
        element_requestor.value = requestor;
    }

    if (address === null || address === undefined || address.trim() == "undefined")
    {
        element_address.value = '';
    }
    else
    {
        element_address.value = address;
    }
}

// [both exportdata.html and register_email.html]
function startEmailCheck(element_address, element_requestor, element_info_msg, callback_update_ui_fn, callback_fn)
{
    // This is called when the check params button is pressed.
    // if not ready, complain and do nothing.
    var address = null;
    var requestor = null;
    var snail_address = null; // no UI yet
    var isRequestorValid = null;
    var isAddressValid = null;

    requestor = element_requestor.value.trim();
    isRequestorValid = ValidateExportRequestor(element_requestor);
    address = element_address.value;
    isAddressValid = ValidateNotificationAddress(element_address);

    if (isRequestorValid && isAddressValid)
    {
        // will avoid making an AJAX call if possible (checked email addresses are cached in element_address.retrieve('addresses')
        // false --> do registration
        register_address(address.slice(0), requestor, snail_address, element_address, element_info_msg, callback_update_ui_fn, callback_fn);
    }
}

function ValidateNotificationAddress(element_address)
{
    var retrievedValid = null;
    var address = element_address.value.trim();
    var valid = true;

    retrievedValid = element_address.retrieve('valid', null);

    if (retrievedValid === null)
    {
        if (address.length == 0)
        {
            valid = false;
        }
        else if (address.indexOf("@") == -1 || address.indexOf("@") == 0 || address.indexOf("@") != address.lastIndexOf("@"))
        {
            element_address.store('error_msg', 'Invalid notification address: missing or invalid domain name');
            valid = false;
        }

        element_address.store('valid', valid);
    }
    else
    {
        valid = retrievedValid;
    }

    return valid;
}

function ValidateExportRequestor(element_requestor)
{
    var requestor = element_requestor.value.trim();
    var retrievedValid = null;
    var valid = true;

    retrievedValid = element_requestor.retrieve('valid', null);

    if (retrievedValid === null)
    {
        if (requestor.indexOf("@") >= 0)
        {
            element_requestor.store('error_msg', 'An email address belongs in the Notify field, not in the Requestor field');
            valid = false;
        }
        else
        {
            // valid!
        }

        element_requestor.store('valid', valid);
    }
    else
    {
        valid = retrievedValid;
    }

    return valid;
}

function call_address_ajax(address, requestor, element_address, callback_update_ui_fn, check_only, registration_callback)
{
    var export_app_arguments = null;
    var method = null;
    var options = null;
    var addresses = null;

    export_app_arguments = { "address" : address, "user-name" : requestor };
    method = check_only ? 'get' : 'post';
    options =
    {
        method: method,
        onSuccess: function(response)
        {
            var ca_status_obj = response.responseJSON
            var status = ca_status_obj.drms_export_status;
            var error_msg = null;
            var address_internal = null;
            var addresses = null;
            var check_msg = null;
            var check_address_timer = null;
            var registration_timer = null;
            var stop_registration_timer = false;
            var callback_list = null;

            address_internal = address.trim();
            addresses = element_address.retrieve('addresses', {});

            if (status.search(/errorcode/) == 0)
            {
                addresses[address_internal].registration_status = 'error';
                error_msg = ca_status_obj.error_message;

                // stop registration timer
                stop_registration_timer = true;
                element_address.store('error_msg', error_msg + '; enter new email address');
            }
            else if (status == 'StatusCode.REGISTRATION_INITIATED')
            {
                // operation was `register` - registration check initiated (email not found in db, and checkOnly == false)
                addresses[address_internal].registration_status = 'registering';
                alert('adding ' + address_internal + ' to addresses')

                element_address.store('seconds_remaining', MAX_REGISTRATION_TIME);
                timer = setInterval(registration_callback, 1000);
                element_address.store('registration_timer', timer);
            }
            else if (status == 'StatusCode.UNREGISTERED_ADDRESS')
            {
                // operation was `check`
                addresses[address_internal].registration_status = false;
            }
            else if (status == 'StatusCode.REGISTERED_ADDRESS')
            {
                // operation was either `check` or `register` - address is registered
                addresses[address_internal].registration_status = true;

                Cookie.setData("user", name);
                Cookie.setData("notify", address);

                // stop registration timer
                stop_registration_timer = true;
            }
            else if (status == 'StatusCode.REGISTRATION_PENDING')
            {
                // operation was `register` - registration pending
                addresses[address_internal].registration_status = 'pending';
            }
            else
            {
                addresses[address_internal].registration_status = 'error';
                error_msg = ca_status_obj.error_message;
                check_msg = 'unexpected error: ' + error_msg

                element_address.store('error_msg', chk_msg);

                // stop registration timer
                stop_registration_timer = true;
            }

            if (callback_update_ui_fn)
            {
                callback_update_ui_fn();
            }

            // clear time-out that cancels check_address endpoint call
            check_address_timer = element_address.retrieve('check_address_timer', null);
            if (check_address_timer !== null)
            {
                clearInterval(check_address_timer);
                element_address.store('check_address_timer', null);
            }

            if (stop_registration_timer)
            {
                registration_timer = element_address.retrieve('registration_timer', null);
                if (registration_timer !== null)
                {
                    clearInterval(registration_timer);
                    element_address.store('registration_timer', null);
                }

                // since registration is complete, successful or not, call the callbacks
                callback_list = element_address.retrieve('callback_list', null);
                if (callback_list)
                {
                    call_callback_list(callback_list);
                    element_address.store('callback_list', null);
                }
            }
        },
        onFailure: function(response)
        {
            var ca_status_obj = response.responseJSON;
            var error_msg = ca_status_obj.error_message;
            var check_address_timer = null;
            var registration_timer = null;
            var callback_list = null;

            element_address.store('error_msg', 'Internal failure checking for email registration: ' + error_msg);

            // clear time-out that cancels check_address endpoint call
            check_address_timer = element_address.retrieve('check_address_timer', null);
            if (check_address_timer !== null)
            {
                clearInterval(check_address_timer);
                element_address.store('check_address_timer', null);
            }

            // stop registration timer
            registration_timer = element_address.retrieve('registration_timer', null);
            if (registration_timer !== null)
            {
                clearInterval(registration_timer);
                element_address.store('registration_timer', null);
            }

            // since registration is complete, successful or not, call the callbacks
            callback_list = element_address.retrieve('callback_list', null);
            if (callback_list)
            {
                call_callback_list(callback_list);
                element_address.store('callback_list', null);
            }
        },
        on428: function(response)
        {
            // arguments to check_address endpoint could not be parsed or were invalid
            // clear time-out that cancels check_address endpoint call
            var ca_status_obj = response.responseJSON
            var error_msg = ca_status_obj.error_message;
            var check_address_timer = null;
            var registration_timer = null;
            var callback_list = null;

            element_address.store('error_msg', 'Invalid arguments to check_address endpoint: ' + error_msg);

            check_address_timer = element_address.retrieve('check_address_timer', null);
            if (check_address_timer !== null)
            {
                clearInterval(check_address_timer);
                element_address.store('check_address_timer', null);
            }

            // stop registration timer
            registration_timer = element_address.retrieve('registration_timer', null);
            if (registration_timer !== null)
            {
                clearInterval(registration_timer);
                element_address.store('registration_timer', null);
            }

            // since registration is complete, successful or not, call the callbacks
            callback_list = element_address.retrieve('callback_list', null);
            if (callback_list)
            {
                call_callback_list(callback_list);
                element_address.store('callback_list', null);
            }
        }
    };

    if (method == 'post')
    {
        options['postBody'] = JSON.stringify(export_app_arguments);
        options['contentType'] = 'application/json';
    }
    else
    {
        options['parameters'] = export_app_arguments;
    }

    addresses = element_address.retrieve('addresses', {});
    addresses[address.trim()].registration_status = 'checking';

    new Ajax.Request(window.location.origin + REGISTRATION_PATH, options);
}

function call_callback_list(callback_list)
{
    var cb = -1;

    if (callback_list !== null)
    {
        for (cb = 0; cb < callback_list.length; cb++)
        {
            if (callback_list[cb])
            {
                (callback_list[cb])();
            }
        }
    }
}

// callback_fn is called once registration status is not pending
function check_or_register_address(address, requestor, snail_address, element_address, element_info_msg, callback_update_ui_fn, callback_fn, check_only)
{
    var statuses = null;
    var registration_status = null;
    var callback_list = null;
    var check_address_callback = null;
    var registration_callback = null;
    var timer = null;
    var pending = false;

    addresses = element_address.retrieve('addresses', {});
    callback_list = element_address.retrieve('callback_list', null);

    // if the email registration status has already been cached, return disposition
    if (Object.prototype.hasOwnProperty.call(addresses, address))
    {
        statuses = addresses[address];

        if (Object.prototype.hasOwnProperty.call(statuses, 'registration_status'))
        {
            registration_status = statuses.registration_status;
        }
    }
    else
    {
        addresses[address] = {};
        statuses = addresses[address];
    }

    if (registration_status !== null && typeof(registration_status) == 'boolean')
    {
        // we have the disposition; shortcut - do not call ajax
        if (callback_fn)
        {
            callback_fn();
        }
    }
    else if (registration_status !== null && typeof(registration_status) == 'string' && (registration_status == 'registering' || registration_status == 'pending'))
    {
        // do nothing, there is a pending registration - do not accept a new callback
    }
    else
    {
        // append the new callback function (regardless if a registration is pending for this address)
        if (callback_list === null)
        {
            callback_list = [];
        }

        if (callback_fn !== null)
        {
            callback_list[callback_list.length] = callback_fn;
            element_address.store('callback_list', callback_list);
        }

        if (true)
        {
            function create_check_address_callback(address)
            {
                check_address_callback = function()
                {
                    var addresses = null;
                    var address = null;
                    var statuses = null;
                    var registration_status = null;
                    var timer = null;
                    var registration_timer = null;
                    var callback_list = null;

                    addresses = element_address.retrieve('addresses', {});
                    timer = element_address.retrieve('check_address_timer', null);
                    registration_timer = element_address.retrieve('registration_timer', null);
                    callback_list = element_address.retrieve('callback_list', null);

                    if (Object.prototype.hasOwnProperty.call(addresses, address))
                    {
                        statuses = addresses[address];
                        registration_status = statuses.registration_status

                        if (typeof(registration_status) === 'string' && (registration_status == 'registering' || registration_status == 'pending' ))
                        {
                            statuses.registration_status = 'timed_out_server';

                            // timeout error message
                            element_address.store('error_msg', 'Timeout waiting for response from registration server');
                        }
                        else
                        {
                            // error
                            element_address.store('error_msg', 'Invalid email-registration state for ' + address);
                        }
                    }

                    // clear interval function
                    if (timer !== null)
                    {
                        clearInterval(timer);
                        element_address.store('check_address_timer', null);
                    }

                    // a registration must be pending, so stop the registration timer and error it out
                    if (registration_timer !== null)
                    {
                        clearInterval(registration_timer);
                        element_address.store('registration_timer', null);
                    }

                    // call callback_list
                    if (callback_list)
                    {
                        call_callback_list(callback_list);
                        element_address.store('callback_list', null);
                    }
                };
            }

            check_address_callback = create_check_address_callback(address);

            if (!check_only)
            {
                // set the registration time-out timer function; if the email is already registered,
                // then checkAddress.py will return status == registered and it will not start the
                // registration process; the success() function will then set seconds_remaining
                // to 0 and this callback will clear the registration function interval
                function create_registration_callback(address)
                {
                    alert('create_registration_callback ' + address);
                    var clean_address = address.trim()
                    registration_callback = function()
                    {
                        var addresses = null;
                        var seconds_remaining = null;
                        var registration_status = null;
                        var check_address_timer = null;
                        var registration_timer = null;
                        var callback_list = null;
                        var error = false;

                        addresses = element_address.retrieve('addresses', {});
                        seconds_remaining = element_address.retrieve('seconds_remaining', null);
                        check_address_timer = element_address.retrieve('check_address_timer', null)
                        registration_timer = element_address.retrieve('registration_timer', null);
                        callback_list = element_address.retrieve('callback_list', null);

                        if (seconds_remaining === null)
                        {
                            // error - seconds_remaining initialized when registration initiated
                            error = true;
                        }
                        else if (!Object.prototype.hasOwnProperty.call(addresses, clean_address))
                        {
                            // error - address must have been added to addresses during registration initiation
                            error = true;
                        }
                        else if (!Object.prototype.hasOwnProperty.call(addresses[clean_address], 'registration_status'))
                        {
                            // error
                            error = true;
                        }
                        else if (typeof(addresses[clean_address].registration_status) !== 'string')
                        {
                            // error
                            error = true;
                        }
                        else if (addresses[clean_address].registration_status != 'registering' && addresses[clean_address].registration_status != 'pending')
                        {
                            // error
                            error = true;
                        }

                        if (!error)
                        {
                            if (seconds_remaining <= 0)
                            {
                                // timeout waiting for registration to complete
                                // set to timed_out so that the callback check_registered_callback_fn knows how to set UI
                                addresses[clean_address].registration_status = 'timed_out_client';

                                // timeout error message
                                element_address.store('error_msg', 'Timeout waiting for response to email message');

                                // gotta call callbacks (registration timed out)
                                if (callback_list)
                                {
                                    call_callback_list(callback_list);
                                    element_address.store('callback_list', null);
                                }

                                if (callback_update_ui_fn)
                                {
                                    callback_update_ui_fn();
                                }

                                if (registration_timer !== null)
                                {
                                    clearInterval(registration_timer);
                                    element_address.store('registration_timer', null);
                                }

                                // reset seconds_remaining (in case the user runs the registration CGI again)
                                element_address.store('seconds_remaining', null);
                            }
                            else
                            {
                                // the registration success function will set seconds_remaining to 0 so the interval function will be cleared
                                element_info_msg.store('cp_message', (seconds_remaining - 1).toString() + ' seconds remaining, waiting for your email reply from ' + clean_address);
                                element_address.store('seconds_remaining', seconds_remaining - 1);

                                if (callback_update_ui_fn)
                                {
                                    callback_update_ui_fn();
                                }

                                // the registration is pending, set up a new check address timer and call check address ajax()
                                // if an existing check_address_timer is set, clear it first
                                if (check_address_timer !== null)
                                {
                                    clearInterval(check_address_timer);
                                    element_address.store('check_address_timer', null);
                                }

                                timer = setInterval(check_address_callback, MAX_NOTIFY_TIMER);
                                element_address.store('check_address_timer', timer);
                                call_address_ajax(clean_address, requestor, element_address, callback_update_ui_fn, check_only, null);
                            }
                        }
                    };

                    return registration_callback;
                }

                registration_callback = create_registration_callback(address);
            }

            // callbacks all set up, now do the actual AJAX if it has not been already initiated
            // if an existing check_address_timer is set, clear it first
            check_address_timer = element_address.retrieve('check_address_timer', null);
            if (check_address_timer !== null)
            {
                clearInterval(check_address_timer);
                element_address.store('check_address_timer', null);
            }

            timer = setInterval(check_address_callback, MAX_NOTIFY_TIMER);
            element_address.store('check_address_timer', timer);
            call_address_ajax(address, requestor, element_address, callback_update_ui_fn, check_only, registration_callback);
        }
    }
}

// register an address if it is not already registered
function register_address(address, requestor, snail_address, element_address, element_info_msg, callback_update_ui_fn, callback_fn)
{
    return check_or_register_address(address, requestor, snail_address, element_address, element_info_msg, callback_update_ui_fn, callback_fn, false);
}

// only check for the registration status of an address
function check_registration(address, requestor, snail_address, element_address, element_info_msg, callback_update_ui_fn, callback_fn)
{
    return check_or_register_address(address, requestor, snail_address, element_address, element_info_msg, callback_update_ui_fn, callback_fn, true);
}

// checkOnly is true if an email registration check is pending (CheckAddressRegistration() was called with checkOnly == false,
// and the checkAddress.sh CGI has not completed the registration process); otherwise, the user has clicked on the "check
// params for export" button after CheckAddressRegistration(checkOnly == false) has been called
//

function NotificationAddressRegistered(element_address)
{
    var addresses;
    var address;

    address = element_address.value.trim();

    if (element_address.retrieve('valid', false))
    {
        // email address is valid, but is it registered?
        addresses = element_address.retrieve('addresses', {});
        if (addresses.hasOwnProperty(address))
        {
            return addresses[address];
        }
    }

    return false;
}
