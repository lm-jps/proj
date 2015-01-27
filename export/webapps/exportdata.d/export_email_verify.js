function email_initVars()
  {
  ExportNotifyOK = 0;
  ExportNotifyTimeLeft = 0;
  ExportNotifyTimeDelta = 2000;
  ExportNotifyValid = -9;
  cookieState = 0;
  }

// Global vars
var JSOC_CHECK_EMAIL = "checkAddress.sh";
var MAX_NOTIFY_TIMER = 180;
var ExportEmail = "not set"; // current email to test
var EXPORTUSER;
var EXPORTNOTIFY;
var ExportNotifyValid;
var ExportNotifyTimer;
var ExportNotifyTimeLeft;  // Max seconds before giveup
var ExportNotifyTimeDelta; // milliseconds between checking
var cookieState; 

function email_getargs()
  {
  cookieState = Cookie.getData("emailOK");
  if (cookieState == "undefined")
    cookieState=0;
  $("ExportRequestor").value = Cookie.getData("user");
  $("ExportNotify").value = Cookie.getData("notify");
  if ($("ExportRequestor").value == "undefined")
    $("ExportRequestor").value = "";
  if ($("ExportRequestor").value == "undefined")
    $("ExportRequestor").value = "solarmail";
  if (cookieState == 2)
    {
    ExportNotifyOK = 1;
    $("ExportNotifyMsg").innerHTML = "OK";
    }
  else
    {
    cookieState=1;
    Cookie.setData("emailOK",cookieState);
    Cookie.setData("user", $("ExportRequestor").value);
    Cookie.setData("notify", $("ExportNotify").value);
    }
  }

function startEmailCheck()
  { // This is called when the check params button is pressed.
    // if not ready, complain and do nothing.
  if ( $("ExportNotify").value.toUpperCase() ===  "SOLARMAIL" )
    ExportEmail = $("ExportRequestor").value + "@spd.aas.org";
  else
    ExportEmail = $("ExportNotify").value;

  $("ExportCheckMsg").innerHTML = "Ready to verify email, press Check Params.";

  if (ExportNotifyOK == 0)
    {
    SetExportUser();
    if (ExportNotifyOK == 0)
      {
      $("ExportCheckMsg").innerHTML = "Notify field must be set before submit.";
      return;
      }
    }
  if (ExportNotifyOK == 1)
    {
    $("ExportCheckMsg").innerHTML = 'The current address, "' + ExportEmail + '", has been verified OK.  You are done.';
    return;
    }
  else if (ExportNotifyOK == 4)
    {
    $("ExportCheckMsg").innerHTML = "Still waiting for reply to our email, wait or change Notify or Requester fields.";
    return;
    } 
  else if (ExportNotifyOK == 5)
    {
    $("ExportCheckMsg").innerHTML = 'Submitted address, "' + ExportEmail + '", is not a valid email address format, fix and try again.';
    return;
    } 
    $("ExportCheckMsg").innerHTML = 'Submitted address is: "' + ExportEmail + '".';
  ExportNotifyOK = CheckNotifyValidity();
  }

// SetExportUser is called from change in either the Requester or Notify input text boxes.
function SetExportUser()
  {  
  if (ExportNotifyOK == 4)
    {
    ExportNotifyOK = 0;
    // Previous verify submit still in progress.
    }
  ExportUserOK = 0;
  EXPORTUSER = $("ExportRequestor").value.toUpperCase();
  EXPORTNOTIFY = $("ExportNotify").value.toUpperCase();
  if (EXPORTUSER.length != 0 && EXPORTUSER !== "NONE")
    {
    if (EXPORTUSER.indexOf("@") >= 0)
      {
      $("RequestorMessage").innerHTML = "Name only here, E-mail address goes in the Notify field";
      ExportUserOK = 0;
      }
    else
      {
      ExportUserOK = 1;
      $("RequestorMessage").innerHTML = "Provide an identifier for you, e.g. your SolarMail name.";
      }
    }
  else
    {
    $("RequestorMessage").innerHTML = "Provide an identifier for you, e.g. your SolarMail name.";
    $("ExportRequestor").value = "none";
    EXPORTUSER = "NONE";
    ExportUserOK = 0;
    }
  SetExportNotify();
  }

function SetExportNotify()
  {   
  // ExportNotifyOK values:
  //  0 - not OK to proceed
  //  1 - OK to proceed, email is validated 
  //  2 - OK, email is SolarMail and Requestor not empty, address not yet validated
  //  3 - OK, quick or direct and no RequestID expected.
  //  4 - Notify address pending validation.
  //  5 - Notify address failed verification within allowed time, MAX_NOTIFY_TIMER seconds.
  //  At exit from SetExportNotify, only 0 or 2 will be set.

  ExportNotifyOK = 0;
  ExportNotifyValid = 0;  // Any change in text boxes invalidates knowledge of email
  var tmp1 = $("ExportNotify").value.toUpperCase();
  if (EXPORTNOTIFY.length == 0)
    {
    ExportNotifyOK = 0;
    $("ExportNotify").value = "solarmail";
    EXPORTNOTIFY = "SOLARMAIL";
    if (EXPORTUSER.length > 0 && EXPORTUSER != "NONE")
      ExportUserOK = 1;
    else
      ExportUserOK = 0;
    }
  if (EXPORTNOTIFY === "SOLARMAIL")
    {
    if (ExportUserOK)
      ExportNotifyOK = 2;  // is SolarMail and Requestor may be valid.
    else
      ExportNotifyOK = 0;
    }
  else if (EXPORTNOTIFY == "NONE" || tmp1 == "NO")
    {
    ExportNotifyOK = 0;
    }
  else if (EXPORTNOTIFY.indexOf("@") > 0 && 
           EXPORTNOTIFY.indexOf("@") == EXPORTNOTIFY.lastIndexOf("@"))
    {
    $("ExportNotifyMsg").innerHTML = "Checking";
    ExportNotifyOK = 2;
    }
  else
    {
    ExportNotifyOK = 0; // this case eliminates two '@' in Notify 
    }
  if ( EXPORTNOTIFY ===  "SOLARMAIL" )
    ExportEmail = $("ExportRequestor").value + "@spd.aas.org";
  else
    ExportEmail = $("ExportNotify").value;
  }

function NotifyTimer()
  {
  ExportNotifyTimeLeft -= ExportNotifyTimeDelta/1000;
$("ExportStatus").innerHTML="ExportNotifyOK="+ExportNotifyOK+" timer="+ExportNotifyTimeLeft;
  CheckNotifyValidity();
  }


function CheckNotifyValidity()
  {
  // this is called when a valid format Notify or Requestor and Notify==solarmail is present.
  // On entry ExportNotifyOK is 2 or 4.
  // If 2 then initiate address check, if 4 then sleep a bit and test for completion of verification.
  // If validity check completed OK then return OK

  var checkOnly = 1; // default
  $("ExportNotifyMsg").innerHTML = "";

  cookieState=1;
  Cookie.setData("emailOK",cookieState);

  if (ExportNotifyOK != 4)
    {
    clearInterval(ExportNotifyTimer);
    }

  if (ExportNotifyOK == 0)
    {
    $("ExportCheckMsg").innerHTML = ExportEmail + " not verified within, fix or try again.";
    $("ExportNotifyMsg").innerHTML = "";
    return;
    }

  if (ExportNotifyValid == 1) // If already verified, done.
    {
    $("ExportNotifyMsg").innerHTML = "OK";
    ExportNotifyOK = 1;
    return;
    }

  if (ExportNotifyOK == 5)
    {
    $("ExportCheckMsg").innerHTML = "Notify address " + ExportEmail + " not verified within " + MAX_NOTIFY_TIMER + " seconds, try again or correct email address then try again.";
    $("ExportNotifyMsg").innerHTML = "";
    ExportNotifyOK = 0;
    return;
    }

  if (ExportNotifyOK == 2)  // New check of validity required
    {
    checkOnly = 0;
    $("ExportNotifyMsg").innerHTML = "";
    ExportNotifyTimeLeft = MAX_NOTIFY_TIMER;
    ExportNotifyTimer = setInterval(function () {NotifyTimer()}, ExportNotifyTimeDelta);
    // timer will run until valid email found or time limit expired or user changes email address
    }
    
  // Always set to pending before starting the verify process
  ExportNotifyOK = 4;

  var paramObj = {"address" : ExportEmail, "checkonly" : checkOnly};
  exportparameters = new Hash(paramObj);
  
  new Ajax.Request('http://' + Host + '/cgi-bin/ajax/' + JSOC_CHECK_EMAIL,
    {
    method: 'post',
    parameters: exportparameters.toObject(),
    onSuccess: function(transport, json)
      {
      // if Notify or Requester has been changed then ExportNotifyOK will
      // have been changed so ignore this now obsolete response.
      if (ExportNotifyOK == 4)
        {
        var response = transport.responseText || null;
        var parseInfo = response.evalJSON();
        var status = parseInfo.status;
        if (status == 1) // check initiated OK
          {
          ExportNotifyValid = 0;
          $("ExportCheckMsg").innerHTML = parseInfo.msg;
          $("ExportNotifyMsg").innerHTML = "";
          }
        else if (status == 2) // address is valid
          {
          clearInterval(ExportNotifyTimer);
          ExportNotifyValid = 1;
          ExportNotifyOK = 1;
          $("ExportCheckMsg").innerHTML = parseInfo.msg;
          $("ExportNotifyMsg").innerHTML = "OK";
          cookieState=2;
          Cookie.setData("emailOK",cookieState);
          Cookie.setData("user", $("ExportRequestor").value);
          Cookie.setData("notify", $("ExportNotify").value);
          }
        else if (status == 3)
          {
          $("ExportNotifyMsg").innerHTML = "";
$("ExportLocation").innerHTML = "status==3, notify="+ExportNotifyOK;
          if (ExportNotifyTimeLeft <= 0)
            { // timer expired, address check failed.
            clearInterval(ExportNotifyTimer);
            ExportNotifyValid = -1;
            ExportNotifyOK = 5;
            $("ExportCheckMsg").innerHTML = "Timeout - A reply was not received in the allowed time.  Try again with a correct Notify address.";
            }
          else
            {
            if (ExportNotifyOK != 4)
              { // User abort of checking due to change of box contents.
              clearInterval(ExportNotifyTimer);
              $("ExportCheckMsg").innerHTML = "Current attempt cancelled by typing, try again when ready.";
              }
            else
              { // keep waiting, check each ExportNotifyTimeDelta seconds
$("ExportLocation").innerHTML = "status==3, notify="+ExportNotifyOK+" setting missing message";
              ExportNotifyValid = 0;
              $("ExportCheckMsg").innerHTML = ExportNotifyTimeLeft + " seconds remaining, still waiting for your email reply from "+EXPORTMAIL;
              }
            }
          }
        else
          {
          ExportNotifyValid = -2; // immediate failure
          ExportNotifyOK = 0;
          clearInterval(ExportNotifyTimer);
          $("ExportNotifyMsg").innerHTML = "";
          $("ExportCheckMsg").innerHTML = 'Notify address provided, "' + ExportEmail + '" is not a valid address, correct and retry.';
          if (status == 4)
            {
            $("ExportCheckMsg").innerHTML += '<br>A prior attempt timed-out before email Reply';
            }
          }
        }
      else
        {
        clearInterval(ExportNotifyTimer);
        ExportNotifyOK = 0;
        $("ExportNotifyMsg").innerHTML = "";
        $("ExportCheckMsg").innerHTML = "Current attempt cancelled by typing, try again when ready.";
        }
      },
      onFailure: function() { alert('oops, our code is broken'); }
    });

  retval = 4;
  return(retval);
  }

