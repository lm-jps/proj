function email_initVars()
{
    // 'valid' indicates if the string currently in the field is valid or not; this gets set only when the field value
    //     changes
    // 'addresses' is an object of validated email address (the properties are email addresses, and the values are booleans that 
    //     indicate a registered email address or not, or 'pending' when checkAddress is in the process of checking
    //     the address registration status):
    //       'addresses' : { 'fred@bedrock.com' : true, 'barney@bedrock.com' : false, 'mrslate@bedrock.com' : 'pending' }
    $("ExportNotify").store({ 'valid' : null, 'addresses' : {} });
    $("ExportRequestor").store({ 'valid' : null });

    ExportNotifyTimeLeft = 0;
    ExportNotifyTimeDelta = 2000;
    
    // we want to hold of on auto-checking notification email address until the user has touched both of these the first time;
    // after that, we can check when they modify one or the other
    UserEnteredNotify = false;
    UserEnteredRequestor = false;
}

// Global vars
var JSOC_CHECK_EMAIL = "checkAddress.sh";
var MAX_NOTIFY_TIMER = 16; // seconds
var ExportNotifyTimer = null;
var ExportNotifyTimeLeft;  // Max seconds before giveup on any AJAX checkAddress call
var ExportNotifyTimeDelta; // milliseconds between checking
var MAX_REGISTRATION_TIME = 180; // seconds
var RegistrationTimer = null;
var RegistrationTimeLeft; // Max seconds before giveup on waiting for pending registration to complete
var UserEnteredNotify;
var UserEnteredRequestor;


function email_getargs()
{
    var cookieState = Cookie.getData("emailOK");
    var addresses = null;
      
    if (cookieState == "undefined")
    {
        cookieState = 0;
    }
 
    // these cookies return undefined if they are not set
    $("ExportRequestor").value = Cookie.getData("user");
    $("ExportNotify").value = Cookie.getData("notify");

    // defaults to use, in case the notification address and/or requestor have not been set in the cookie    
    if ($("ExportNotify").value == "undefined")
    {
        $("ExportNotify").value = "solarmail";
    }
    
    if ($("ExportRequestor").value == "undefined")
    {
        $("ExportRequestor").value = "";
    }
    
    if (cookieState == 2)
    {
        // valid and registered email address (we store only a valid and registered email address and requestor in the cookie)
        addresses = $("ExportNotify").retrieve('addresses', {});
        addresses[$("ExportNotify").value.trim()] = true;
        $("ExportNotify").store({ 'valid' : true, 'addresses' : addresses });
        $("ExportRequestor").store('valid', true);
        $("ExportNotifyMsg").style.color = colorDarkGreen;
        $("ExportNotifyMsg").innerHTML = "OK";
        $("ExportCheckMsg").innerHTML = "Email address " + $("ExportNotify").value.trim() + " is registered.";
    }
    else
    {
        // unregistered email in cookie (if any address in cookie) - make sure the state is 1, and not 0 or undefined
        Cookie.setData("emailOK", 1);
    }
}

function startEmailCheck(callbackFxn)
{   // This is called when the check params button is pressed.
    // if not ready, complain and do nothing.
    var requestor = null;
    var isRequestorValid = null;
    var address = null;
    var isAddressValid = null;
    
    
    requestor = $("ExportRequestor").value.trim();
    isRequestorValid = ValidateExportRequestor();
    address = $("ExportNotify").value;    
    isAddressValid = ValidateNotificationAddress();
    
    if (isRequestorValid && isAddressValid)
    {
        $("ExportCheckMsg").innerHTML = 'Submitted email address is: "' + address + '".';
        
        // will avoid making an AJAX call if possible (checked email addresses are cached in $("ExportNotify").retrieve('addresses')
        CheckAddressRegistration(address.slice(0), callbackFxn, true);
    }
}

function ValidateNotificationAddress()
{
    var retrievedValid = null;
    var address = $("ExportNotify").value.trim();
    var isSolarMailAddress = false;
    var valid = true;

    
    retrievedValid = $("ExportNotify").retrieve('valid', null);
    
    if (retrievedValid === null)
    {
        if (address.indexOf("@") == -1 || address.indexOf("@") == 0 || address.indexOf("@") != address.lastIndexOf("@"))
        {
            $("ExportCheckMsg").innerHTML = "Invalid notification address: missing or invalid domain name.";
            valid = false;
        }
        
        $("ExportNotify").store('valid', valid);
    }
    else
    {
        valid = retrievedValid;
    }

    return valid;
}

function ValidateExportRequestor()
{
    var address = $("ExportNotify").value.trim();
    var requestor = $("ExportRequestor").value.trim();
    var retrievedValid = null;
    var isSolarMailAddress = false;
    var valid = true;


    retrievedValid = $("ExportRequestor").retrieve('valid', null);
    
    if (retrievedValid === null)
    {
        if (address.length == 0 || address.toLowerCase() === 'solarmail')
        {
            isSolarMailAddress = true;
        }
    
        if (requestor.toLowerCase() === 'none' || requestor.toLowerCase() === "no" || requestor.length == 0)
        {
            if (isSolarMailAddress)
            {
                $("ExportCheckMsg").innerHTML = "Invalid notification address: Requestor cannot be 'none'.";
                valid = false;
            }
        }
        else if (requestor.indexOf("@") >= 0)
        {
            $("RequestorMessage").innerHTML = "An email address belongs in the Notify field, not in the Requestor field.";
            valid = false;
        }
        else 
        {
            // valid!
        }

        $("ExportRequestor").store('valid', valid);
    }
    else
    {
        valid = retrievedValid;
    }

    return valid;
}

// This function gets called when the Notify field gets modified.
function SetExportNotify(clickedByUser)
{
    var address = null;
    var requestor = null;
    var doValidation = true;
    var valid = false;

    
    if (clickedByUser)
    {
        UserEnteredNotify = true;
    }
    
    // invalidate address
    $("ExportNotify").store('valid', null);
    
    address = $("ExportNotify").value.trim();
    requestor = $("ExportRequestor").value.trim();
    if (address.toLowerCase() === 'solarmail')
    {
        $("RequestorMessage").innerHTML = "Provide your SolarMail name.";
        
        // if the user already entered a Requestor, modify the Notify, then validate and register
        if (UserEnteredRequestor && requestor.length > 0)
        {
            address = $("ExportRequestor").value.trim() + "@spd.aas.org";
            $("ExportNotify").value = address;
        }
        else
        {
            doValidation = false;
        }
    }
    else
    {
        $("RequestorMessage").innerHTML = "Provide an optional identifier."
    }
    
    if (doValidation)
    {
        valid = ValidateNotificationAddress();
        if (valid)
        {
            $("ExportNotify").store('valid', true);
            CheckAddressRegistration(address.slice(0), null, false);
        }
    }

}

// SetExportUser is called from change in either the Requester or Notify input text boxes.
// There is no call made to check the validity of the notification email address provided
// however. Instead, the user must click on the "Check Params for Export" button.
function SetExportUser(clickedByUser)
{ 
    var requestor = null;
    var address = null;
    var newAddress = null;
    var valid = false;
    var isSolarMailAddress = false;


    if (clickedByUser)
    {
        UserEnteredRequestor = true;
    }
    
    // invalidate requestor
    $("ExportRequestor").store('valid', null);

    requestor = $("ExportRequestor").value.trim();
    address = $("ExportNotify").value.trim();
    
    valid = ValidateExportRequestor();
    if (valid)
    {
        // we may be able to form a SolarMail address
        if (address.length == 0 || address.toLowerCase() === 'solarmail')
        {
            isSolarMailAddress = true;
        }

        if (isSolarMailAddress)
        {
            if (requestor.length == 0)
            {
                if (UserEnteredNotify)
                {
                    $("ExportCheckMsg").innerHTML = "Invalid solarmail address: Requestor must be specified.";
                }
            }
            else
            {
                newAddress = $("ExportRequestor").value.trim() + "@spd.aas.org";
                
                if (address.toLowerCase() !== newAddress.toLowerCase())
                {
                    $("ExportNotify").value = newAddress;

                    // call ExportNotify onChange()
                    $("ExportNotify").onchange();
                }
            }
        }
    }
}

function CallCheckAddressRegistration(address)
{
    ExportNotifyTimeLeft -= ExportNotifyTimeDelta/1000;
    
    // always do a simple check for the registration or check process to be complete; do not start a new registration
    // process
    CheckAddressRegistration(address, null, true);
}

function RegistrationTimeout(address)
{
    RegistrationTimeLeft -= 1;
    
    if (RegistrationTimeLeft <= MAX_REGISTRATION_TIME)
    {
        $("ExportCheckMsg").innerHTML = RegistrationTimeLeft + " seconds remaining, still waiting for your email reply from " + address;
    }

    if (RegistrationTimeLeft <= 0)
    {
        // timer expired, address check failed.
        addresses = $("ExportNotify").retrieve('addresses', {});
        addresses[address] = 'timeout';
        $("ExportNotify").store( 'addresses', addresses );
        $("ExportNotifyMsg").style.color = colorDarkRed;
        $("ExportNotifyMsg").innerHTML = "Timeout";
        $("ExportCheckMsg").innerHTML = "Timeout - the registration process did not complete in the allowed time. Try again with a correct Notify address.";
        
        clearInterval(RegistrationTimer);
        RegistrationTimer = null;
    }
}

// checkOnly is true if an email registration check is pending (CheckAddressRegistration() was called with checkOnly == false, 
// and the checkAddress.sh CGI has not completed the registration process); otherwise, the user has clicked on the "check
// params for export" button after CheckAddressRegistration(checkOnly == false) has been called
//
// it is possible that the user changed the email address while a registration process or check was pending; make sure
// to examine ExportNotifyTimer - if it not null and the user is attempting to register a different address
function CheckAddressRegistration(address, callbackFxn, checkOnly)
{
    var cgiArgs = null;
    var pending = false;
    
    
    // if this is a new registration (checkOnly == false), then end the existing ExportNotifyTimer (yes, orphan it) and let it
    // play out; then go ahead and pretend that there is not an existing registration pending
    if (!checkOnly)
    {
        if (ExportNotifyTimer !== null)
        {
            clearInterval(ExportNotifyTimer);
            ExportNotifyTimer = null;
        }
        
        if (RegistrationTimer !== null)
        {
            // stop the countdown-till-timeout message
            clearInterval(RegistrationTimer);
            RegistrationTimer = null;
            $("ExportCheckMsg").innerHTML = "";
        }
    }

    // this is called when a valid format Notify or Requestor and Notify==solarmail is present.
    // address = $("ExportNotify").value.trim();
  
    // make sure the email address is valid
    if (!$("ExportNotify").retrieve('valid', false))
    {
        if (ExportNotifyTimer !== null)
        {
            clearInterval(ExportNotifyTimer);
            ExportNotifyTimer = null;
        }

        return;
    }
  
    // check to see if the email address already got registration-checked (the results are cached)
    addresses = $("ExportNotify").retrieve('addresses', {});
    if (addresses.hasOwnProperty(address))
    {
        // the AJAX call has completed 
        if (ExportNotifyTimer !== null)
        {
            clearInterval(ExportNotifyTimer);
            ExportNotifyTimer = null;
        }

        // email has already been checked for registration
        if (typeof(addresses[address]) === 'string')
        {
            if (addresses[address] === 'pending')
            {
                // we can remove timer since the asynchronous call has completed - this is done below
                pending = true;
            }
            else if (addresses[address] === 'timeout')
            {
                // registration process timeout (not timeout due to no response from AJAX call)
                return;
            }
        }
        else if (typeof(addresses[address]) === 'boolean')
        {
            if (addresses[address])
            {
                // registered 
                $("ExportNotifyMsg").style.color = colorDarkGreen;
                $("ExportNotifyMsg").innerHTML = "OK";
                $("ExportCheckMsg").innerHTML = "Email address " + address + " is registered.";
            }
            else
            {
                // not registered
                $("ExportNotifyMsg").style.color = colorDarkRed;
                $("ExportNotifyMsg").innerHTML = "Failed";
                $("ExportCheckMsg").innerHTML = address + ' is not a registered address - please provide a registered address and retry.';
            }
            
            // we have our answer
            return;
        }
        else
        {
            // error
            $("ExportCheckMsg").innerHTML = address + ' invalid email-registration state.';
            return;
        }        
    }
    else if (ExportNotifyTimer !== null)
    {
        // AJAX call has not completed; see if it has been too long
        if (ExportNotifyTimeLeft <= 0)
        {
            // time-out waiting for any response from checkAddress server
            addresses = $("ExportNotify").retrieve('addresses', {});
            addresses[address] = 'timeout';
            $("ExportNotify").store('addresses', addresses);
            $("ExportNotifyMsg").style.color = colorDarkRed;
            $("ExportNotifyMsg").innerHTML = "Timeout";
            $("ExportCheckMsg").innerHTML = "Timeout - A reply from the address-registration server was not received in the allowed time.";
            
            if (ExportNotifyTimer !== null)
            {
                clearInterval(ExportNotifyTimer);
                ExportNotifyTimer = null;
            }
            
            // do not continue checking for registration completion
            return;
        }
    }

    if (ExportNotifyTimer === null)
    {
        // we need to do either a registration check (checkOnly == true), or a check followed by a registration (checkOnly == false)
        ExportNotifyTimeLeft = MAX_NOTIFY_TIMER;
        ExportNotifyTimer = setInterval(function () { CallCheckAddressRegistration(address.slice(0)); }, ExportNotifyTimeDelta);
        // timer will run until valid email found or time limit expired or user changes email address
        // *** JS will NOT call CallCheckAddressRegistration() until after CheckAddressRegistration() because JS does not 
        // interrupt function calls to switch to a different function call; so we do not need to worry about entering 
        // CheckAddressRegistration() while CheckAddressRegistration() is already executing
    }
    else
    {
        // there is a either registration check or registration process in progress - do not issue a new one        
        return;
    }
    
    if (!pending)
    {
        // if we are doing a fresh check or check + register, then set to Checking..., but if we are checking on a pending 
        // check or check + register, then $("ExportNotifyMsg").value will already be appropriate (e.g., Registering...)
        $("ExportNotifyMsg").style.color = colorDarkBlue;
        $("ExportNotifyMsg").innerHTML = "Checking...";
    }

    cgiArgs = { "address" : address, "checkonly" : checkOnly ? 1 : 0 };
    
    new Ajax.Request('http://' + Host + '/cgi-bin/ajax/' + JSOC_CHECK_EMAIL,
    {
        method: 'post',
        parameters: cgiArgs,
        onSuccess: function(transport, json)
        {
            var response = transport.responseText || null;
            var parseInfo = response.evalJSON();
            var status = parseInfo.status;
            var addressInternal = null;


            addressInternal = address.slice(0);
            if (status == 1)
            {
                alert('starting a reg');
                // registration check initiated (email not found in db, and checkOnly == false)
                RegistrationTimeLeft = MAX_REGISTRATION_TIME + 5; // 5 seconds so the user can read the message from checkAddress
                RegistrationTimer = setInterval(function () { RegistrationTimeout(addressInternal) }, 1000);
                
                addresses = $("ExportNotify").retrieve('addresses', {});
                addresses[addressInternal] = 'pending';
                $("ExportNotify").store( 'addresses', addresses );
                $("ExportCheckMsg").innerHTML = parseInfo.msg;
                $("ExportNotifyMsg").style.color = colorDarkBlue;
                $("ExportNotifyMsg").innerHTML = "Registering...";
                
            }
            else if (status == 2)
            {
                // address is valid
                addresses = $("ExportNotify").retrieve('addresses', {});
                addresses[addressInternal] = true;
                $("ExportNotify").store( 'addresses', addresses );
                $("ExportNotifyMsg").style.color = colorDarkGreen;
                $("ExportNotifyMsg").innerHTML = "OK";
                // $("ExportCheckMsg") is set during callback

                Cookie.setData("emailOK", 2);
                Cookie.setData("user", $("ExportRequestor").value);
                Cookie.setData("notify", $("ExportNotify").value);
            }
            else if (status == 3)
            {
                // registration pending            
                $("ExportNotifyMsg").innerHTML = "";

                // keep waiting, check each ExportNotifyTimeDelta seconds
                addresses = $("ExportNotify").retrieve('addresses', {});
                addresses[addressInternal] = 'pending';
                $("ExportNotify").store( 'addresses', addresses );
                $("ExportNotifyMsg").style.color = colorDarkBlue;
                $("ExportNotifyMsg").innerHTML = "Pending...";
            }
            else 
            {
                // status == 4 ==> invalid address                
                addresses = $("ExportNotify").retrieve('addresses', {});
                addresses[addressInternal] = false;
                $("ExportNotify").store( 'addresses', addresses );
                $("ExportNotifyMsg").style.color = colorDarkRed;
                $("ExportNotifyMsg").innerHTML = "Failed";
                $("ExportCheckMsg").innerHTML = 'Notify address provided, "' + addressInternal + '" is not a valid address, correct and retry.';
                if (status == 4)
                {
                    $("ExportCheckMsg").innerHTML += '<br>A prior attempt timed-out before email Reply';
                }
            }

            if (callbackFxn)
            {
                callbackFxn();
            }
        },
        onFailure: function() 
        { 
            alert('oops, our code is broken'); 
            
            $("ExportCheckMsg").innerHTML = 'Internal failure checking for email registration.';
            
            if (ExportNotifyTimer !== null)
            {
                clearInterval(ExportNotifyTimer);
                ExportNotifyTimer = null;
            }
        }
    });

  retval = 4;
  return(retval);
  }

function NotificationAddressRegistered()
{
    var addresses;
    var address;
    
    address = $("ExportNotify").value.trim();
    
    if ($("ExportNotify").retrieve('valid', false))
    {
        // email address is valid, but is it registered?
        addresses = $("ExportNotify").retrieve('addresses', {});
        if (addresses.hasOwnProperty(address))
        {
            return addresses[address];
        }
    }
    
    return false;
}
