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

function call_address_ajax(address, requestor, element_address, callback_update_ui_fn, check_only)
{
    var cgiArgs = null;

    // do AJAX
    cgiArgs = { "address" : address, "name" : requestor, "checkonly" : check_only ? 1 : 0 };

    new Ajax.Request('http://' + Host + '/cgi-bin/ajax/' + JSOC_CHECK_EMAIL,
    {
        method: 'post',
        parameters: cgiArgs,
        onSuccess: function(response)
        {
            var ca_status_obj = response.responseJSON
            var status = ca_status_obj.status;
            var error_msg = null;
            var address_internal = null;
            var addresses = null;
            var check_msg = null;
            var check_address_timer = null;

            address_internal = address.trim();
            addresses = element_address.retrieve('addresses', {});

            if (status < 0)
            {
                addresses[address_internal].registration_status = 'error';
                error_msg = ca_status_obj.msg;

                // stop registration timer
                element_address.store('seconds_remaining', 0);

                element_address.store('error_msg', error_msg + '; enter new email address');
            }
            else if (status == 1)
            {
                // registration check initiated (email not found in db, and checkOnly == false)
                addresses[address_internal].registration_status = 'registering';
            }
            else if (status == 2)
            {
                // address is valid
                addresses[address_internal].registration_status = true;

                // ART - temp so I don't have to delete the cookie a milion times
                // Cookie.setData("user", name);
                // Cookie.setData("notify", address);

                // stop registration timer
                element_address.store('seconds_remaining', 0);
            }
            else if (status == 3)
            {
                // registration pending (email not found in db, and checkOnly == false)
                addresses[address_internal].registration_status = 'pending';
            }
            else
            {
                // status == 4 ==> invalid address
                addresses[address_internal].registration_status = false;
                check_msg = 'Notify address provided, "' + address_internal + '" is not a valid address, correct and retry.';

                if (status == 4)
                {
                    check_msg += '<br>A prior attempt timed-out before email Reply';
                }

                element_address.store('error_msg', chk_msg);

                // stop registration timer
                element_address.store('seconds_remaining', 0);
            }

            if (callback_update_ui_fn)
            {
                callback_update_ui_fn();
            }

            // clear time-out timer (since success was called)
            check_address_timer = element_address.retrieve('check_address_timer', null);
            if (check_address_timer !== null)
            {
                clearInterval(check_address_timer);
                element_address.store('check_address_timer', null);
            }
        },
        onFailure: function()
        {
            alert('oops, our code is broken');
            element_address.store('error_msg', 'Internal failure checking for email registration');
        }
    });
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

function check_or_register_address(address, requestor, snail_address, element_address, element_info_msg, callback_update_ui_fn, callback_fn, check_only)
{
    var statuses = null;
    var registration_status = null;
    var registration_pending = null;
    var callback_list = null;
    var check_address_callback = null;
    var registration_callback = null;
    var timer = null;
    var pending = false;

    addresses = element_address.retrieve('addresses', {});
    callback_list = element_address.retrieve('callback_list', null);
    registration_pending = (callback_list !== null);

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

        // if pending, then we only add the callback to the list of callbacks
        if (!registration_pending)
        {
            registration_status = 'registering';
            statuses.registration_status = registration_status;

            check_address_callback = function()
            {
                var addresses = null;
                var address = null;
                var statuses = null;
                var registration_status = null;
                var timer = null;
                var callback_list = null;

                addresses = element_address.retrieve('addresses', {});
                timer = element_address.retrieve('check_address_timer', null);
                callback_list = element_address.retrieve('callback_list', null);

                // email has already been checked for registration already
                for (address in addresses)
                {
                    // address is a reference
                    if (Object.prototype.hasOwnProperty.call(addresses, address))
                    {
                        statuses = addresses[address];
                        registration_status = statuses.registration_status;

                        if (typeof(registration_status) === 'string')
                        {
                            statuses.registration_status = 'timed_out_server';

                            // timeout error message
                            element_address.store('error_msg', 'Timeout waiting for response from registration server');
                        }
                        else if (typeof(registration_status) === 'boolean' && (registration_status == true || registration_status == false))
                        {
                            // registration completed; the success function will call the callbacks if it gets called before the
                            // timeout (for check_only == true only)
                        }
                        else
                        {
                            // error
                            element_address.store('error_msg', 'Invalid email-registration state for ' + address);
                        }
                    }
                }

                // clear interval function
                if (timer !== null)
                {
                    clearInterval(timer);
                    element_address.store('check_address_timer', null);
                }

                // must stop the registration process too
                if (!check_only)
                {
                    // stop registration timer
                    element_address.store('seconds_remaining', 0);

                    // call callback_list
                    if (callback_list)
                    {
                        call_callback_list(callback_list);
                        element_address.store('callback_list', null);
                    }
                }
            };

            // set time-out timer function, which will be called only if a time-out happens; otherwise,
            // the success() function will run instead; if the success() runs first, it will
            // clear the time-out function

            // if called, this will orphan the check/register call
            timer = setInterval(check_address_callback, MAX_NOTIFY_TIMER);
            element_address.store('check_address_timer', timer);

            if (!check_only)
            {
                // set the registration time-out timer function; if the email is already registered,
                // then checkAddress.py will return status == registered and it will not start the
                // registration process; the success() function will then set seconds_remaining
                // to 0 and this callback will clear the registration function interval
                registration_callback = function()
                {
                    var addresses = null;
                    var seconds_remaining = null;
                    var statuses = null;
                    var registration_status = null;
                    var check_address_timer = null;
                    var registration_timer = null;
                    var callback_list = null;
                    var clear_interval = false;

                    //addresses = element_address.retrieve('addresses', {});
                    addresses = element_address.retrieve('addresses', {});
                    seconds_remaining = element_address.retrieve('seconds_remaining', MAX_REGISTRATION_TIME);
                    check_address_timer = element_address.retrieve('check_address_timer', null)
                    registration_timer = element_address.retrieve('registration_timer', null);
                    callback_list = element_address.retrieve('callback_list', null);

                    if (seconds_remaining === null)
                    {
                        seconds_remaining = MAX_REGISTRATION_TIME;
                        element_address.store('seconds_remaining', seconds_remaining);
                    }

                    if (seconds_remaining <= 0)
                    {
                        // timeout OR complete (reg status is true/false)
                        for (address in addresses)
                        {
                            // address is a reference
                            if (Object.prototype.hasOwnProperty.call(addresses, address))
                            {
                                statuses = addresses[address];
                                registration_status = statuses.registration_status;

                                if (typeof(registration_status) === 'string')
                                {
                                    if (registration_status == 'registering' || registration_status == 'pending')
                                    {
                                        // set to timed_out so that the callback check_registered_callback_fn knows how
                                        // to set UI
                                        statuses.registration_status = 'timed_out_client';

                                        // timeout error message
                                        element_address.store('error_msg', 'Timeout waiting for response to email message');
                                    }
                                    else
                                    {
                                        // error during registration process (registration_status should be error - keep it error)
                                    }
                                }
                                else if (typeof(registration_status) === 'boolean' && (registration_status == true || registration_status == false))
                                {
                                    // registration completed already, do nothing; onsuccess()
                                }
                                else
                                {
                                    // error
                                }
                            }
                        }

                        // gotta call callbacks (reg is complete, success/failure/time_out)
                        if (callback_list)
                        {
                            call_callback_list(callback_list);
                            element_address.store('callback_list', null);
                        }

                        clear_interval = true;

                        // reset seconds_remaining (in case the user runs the registration CGI again)
                        element_address.store('seconds_remaining', null);
                    }
                    else
                    {
                        // the registration success function will set seconds_remaining to 0 so the interval function will be cleared
                        element_info_msg.store('cp_message', (seconds_remaining - 1).toString() + ' seconds remaining, waiting for your email reply from ' + address);
                        element_address.store('seconds_remaining', seconds_remaining - 1);

                        if (callback_update_ui_fn)
                        {
                            callback_update_ui_fn();
                        }

                        if (check_address_timer === null)
                        {
                            if (callback_list === null)
                            {
                                // done with registration
                                clear_interval = true;
                            }
                            else
                            {
                                // if check_address_timer is null, then the success() function ran;
                                // if callback_list is null, then we are completely done with a registration,
                                // otherwise registration is pending;
                                // if registration is pending, set up a new check address timer and call check address ajax()
                                timer = setInterval(check_address_callback, MAX_NOTIFY_TIMER);
                                element_address.store('check_address_timer', timer);
                                call_address_ajax(address, requestor, element_address, callback_update_ui_fn, check_only);
                            }
                        }
                    }

                    if (clear_interval)
                    {
                        if (registration_timer !== null)
                        {
                            clearInterval(registration_timer);
                            element_address.store('registration_timer', null);
                        }
                    }
                };

                // if called, this will update
                element_address.store('seconds_remaining', MAX_REGISTRATION_TIME);
                timer = setInterval(registration_callback, 1000);
                element_address.store('registration_timer', timer);
            }

            // callbacks all set up, now do the actual AJAX if it has not been already initiated
            call_address_ajax(address, requestor, element_address, callback_update_ui_fn, check_only);
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
