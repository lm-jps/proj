// Functions for export options
//
// CheckXXXX creates the contents of ExportProcessingArgs from parameters for
// this export type
//
// Each processing type needs a common name, e.g. ImPatch, used for the first part of
// function and variable names associated with that processing type.
// For processing 'Xxxx' the following functions should be defined in this file.
//
//  XxxxInit(isActive) - if isActive == 0 reset all visible elements to default state, reset associated global variables.
//                       if isActive == 1 reset only quantities needed between user actions, if any.
//                       XxxxInit(0) should be called from ProcessingInit() at the end of this file to initialize
//                       the html for this processing type.
//  XxxxCheck() - Examine all user set parameters for correctness and consistency, mark incorrect entries with color.
//                Build string to use in call to jsoc_fetch in the "process" parameter for this processing type.
//                this string will start with the process name, e.g. "hg_patch" and contains all needed command line
//                params with comma delimiters.  Flags should be e.g. "c=1" rather than "-c".
//                Returns the string if all is OK, or "" if user action is still required.
//  XxxxSet(xx) - Optional, may be used in the onChange or onClick events to accept user input,
//                or the XxxxCheckX function may be used for this purpose.
//  Other functions may be defined as needed and may be invoked in the XxxxSet() function where this processing type
//                is tested for completion, but only the XxxxInit and CheckXxxx funtions will be used elsewhere.
//
//  Finally, make entries in ProcessingInit() at the end of this file.
//                

//
// Processing details for AIA normalized Scaling
//

function AiaScaleSet()
  {
  }

function AiaScaleInit(isActive)
  {
  if (!isActive)
    {
    $("ProcessAiaScale").style.display="none";
    }
  }

function AiaScaleCheck()
  {
  var isok = 1;
  var args = "aia_scale";
  if (SeriesName === 'aia.lev1')
    {
    $("ExportFilenameFmt").value = $("ExportFilenameFmt").value.replace('T_REC','T_OBS');
    }
  else
    {
    isok = 0;
    alert("Error - aia_scale can only be used for series aia.lev1");
    }
  ExportProcessingOK = isok;
  CheckRediness();
  return (isok ? args : "");
  }

// End AiaScale

//
// Processing details for rebin using jsoc_rebin
//

function RebinInit(isActive)
  {
  if (!isActive)
    {
    // Clear out any old settings.
    $("RebinMethod").selectedIndex=0;
    $("RebinCrop").checked = false;
    $("RebinRotate").checked = false;
    $("RebinScale").value = "1.0";
    $("RebinFWHM").value = "-1.0";
    $("RebinNvector").value = "-1.0";
    $("RebinFWHM").style.backgroundColor = colorRed;
    $("RebinNvector").style.backgroundColor = colorRed;
    $("RebinFWHMRow").style.display = "none";
    $("RebinNvectorRow").style.display = "none";
    $("ProcessRebin").style.display = "none";
    
    // Display the rebin options table.
    }
  }

function RebinSet(control)
  {
    if (control == "method")
    {
        if ($("RebinMethod").selectedIndex==1)
        {
            $("RebinFWHMRow").style.display = "table-row";
            $("RebinNvectorRow").style.display = "table-row";
        }
        else if ($("RebinMethod").selectedIndex==0)
        {
            $("RebinFWHMRow").style.display = "none";
            $("RebinNvectorRow").style.display = "none";
        }
    }
    else if (control == "fwhm")
    {
        if (parseFloat($("RebinFWHM").value) == -1.0)
        {
            $("RebinFWHM").style.backgroundColor = colorRed;
        }
        else
        {
            $("RebinFWHM").style.backgroundColor = colorWhite;
        }
    }
    else if (control == "nvector")
    {
        if (parseFloat($("RebinNvector").value) == -1.0)
        {
            $("RebinNvector").style.backgroundColor = colorRed;
        }
        else
        {
            $("RebinNvector").style.backgroundColor = colorWhite;
        }
    }
  RebinCheck();
  }

function RebinCheck()
  {
  var rv = "rebin";
  var isok = 1;
    
  if ($("RebinCrop").checked)
    {
    rv = rv + ",c=1";
    }
    
  if ($("RebinRotate").checked)
    {
    rv = rv + ",u=0";
    }
  else
    {
    rv = rv + ",u=1";
    }
    
  if (parseFloat($("RebinScale").value) != 1.0)
    {
    if ($("RebinScale").value <= 0)
      isok = 0;
    rv = rv + ",scale=" + $("RebinScale").value;
    }
    
  if ($("RebinMethod").selectedIndex == 1)
    {
    rv = rv + ",method=gaussian";
    if (parseFloat($("RebinFWHM").value) <= 0)
      {
      isok = 0;
      }
    else
      {
      rv = rv + ",FWHM=" + $("RebinFWHM").value;
      }
        
    if (parseFloat($("RebinNvector").value) == -1.0)
      {
      isok = 0;
      }
    else
      {
      rv = rv + ",nvector=" +  $("RebinNvector").value;
      }
    }
  else
    {
    rv = rv + ",method=boxcar";
    }
    
  ExportProcessingOK = isok;
  CheckRediness();
  return (isok ? rv : "");
  }

// End Rebin

//
// Processing details for rebin using jsoc_resize
//

function ResizeInit(isActive)
  {
  if (!isActive)
    {
    // Clear out any old settings.
    $("ResizeScaleRow").style.display = "table-row";
    $("ResizeBicubic").checked = true;
    $("ResizeSunCenter").checked = true;
    $("ResizeCrop").checked = false;
    $("ResizeReplicate").checked = false;
    $("ResizeDoScale").checked = false;
    $("ResizeCdelt").value = "-1.0";
    $("ResizeCdelt").style.backgroundColor = colorRed;
    $("ResizeScaleRow").style.display = "none";
    }
  }

function ResizeSet(control)
  {
  if (control == "do_scale")
    {
    if ($("ResizeDoScale").checked)
      {
      $("ResizeScaleRow").style.display = "table-row";
      }
    else
      {
      $("ResizeScaleRow").style.display = "none";
      }
    }
  else if (control == "scale")
    {
    $("ResizeCdelt").style.backgroundColor = colorWhite;
    }
  // no specific action needed for "method", "register_to", "crop", or "replicate" 
  ResizeCheck();
  }

function ResizeCheck()
  {
  var rv = "resize";
  var isok = 1;
    
  if ($("ResizeBicubic").checked)
    rv += ",regrid=1";
  else
    rv += ",regrid=0";

  if ($("ResizeReplicate").checked)
    rv += ",do_stretchmarks=1";
  else
    rv += ",do_stretchmarks=0";
    
  if ($("ResizeSunCenter").checked)
    rv += ",center_to=0";
  else if ($("ResizeFirstImage").checked)
    rv += ",center_to=1";
  else
    rv += ",center_to=2";
    
  if ($("ResizeDoScale").checked)
    {
    rv = rv + ",rescale=1";
    rv = rv + ",scale_to=" + $("ResizeCdelt").value;
    if ($("ResizeCdelt").value <= 0)
      isok = 0;
    }
  else
    rv = rv + ",rescale=0";

  if ($("ResizeCrop").checked)
    {
    rv = rv + ",c=1";
    }
    
  ExportProcessingOK = isok;
  CheckRediness();
  return (isok ? rv : "");
  }

// End Resize

//
// Process ImPatch
//

var noaaColor;
var ImTracked = 0;

function ImPatchGetNoaa()
  {
  if ($("ImNOAA").value.strip().empty()) $("ImNOAA").value = "NotSpecified";
  if ($("ImNOAA").value != "NotSpecified")
    {
    var noaaNum = 1 * $("ImNOAA").value;
    if (noaaNum < 7000) 
      {
      noaaNum = noaaNum + 10000; // OK for times after 1996 Jan.
      $("ImNOAA").value = noaaNum + "";
      }
    $("AjaxBusy").innerHTML = Ajax.activeRequestCount;
    new Ajax.Request('http://' + Host + '/cgi-bin/ajax/jsoc_info_jsoc2',
      {
      method: 'get',
      parameters: {"op" : "rs_list", "ds": "su_rsb.NOAA_ActiveRegions[][" + noaaNum + "]", "key": "ObservationTime,LatitudeHG,LongitudeCM" },
      onSuccess: function(transport, json)
        {
        var response = transport.responseText || "no response text";
        var NOAA_rslist = response.evalJSON();
        try {if (NOAA_rslist.status > 0 || NOAA_rslist.count == 0) throw "noRecords";}
        catch(err) { $("ImNOAA").value = noaaNum + " " + err; return; }
        var minLong = 999, minLat, minTime;
        var irec, nrecs = NOAA_rslist.count;
        var thisTime, thisLong, thisLat;
        for (irec=0; irec<nrecs; irec++)
          {
          thisTime = NOAA_rslist.keywords[0].values[irec];
          thisLat = NOAA_rslist.keywords[1].values[irec];
          thisLong = NOAA_rslist.keywords[2].values[irec];
          if (Math.abs(thisLong) < Math.abs(minLong))
            {
            minLong = thisLong;
            minLat = thisLat;
            minTime = thisTime;
            }
          }
        try
          {
          if (minLong == 999) throw "noRegion";
          $("ImLocType").selectedIndex = 0;
          $("ImTRef").value = minTime;
          $("ImX").value = minLong + "";
          $("ImY").value = minLat + "";
	  $("ImNOAA").style.backgroundColor=colorOptionSet;
          noaaColor = colorPreset;
          $("ImLocType").style.backgroundColor = noaaColor;
          $("ImTRef").style.backgroundColor = noaaColor;
          $("ImX").style.backgroundColor = noaaColor;
          $("ImY").style.backgroundColor = noaaColor;
	  // CheckImPatch();
          }
        catch(err) { $("ImNOAA").value = noaaNum + " " + err; return; }
        },
      onFailure: function()
        {
        alert('Something went wrong with NOAA num data request');
        $("ImNOAA").value = "Not Found";
        },
      onComplete: function() { $("AjaxBusy").innerHTML = Ajax.activeRequestCount; }
      });
    }
  }
  
function ImPatchSet(param)
  {
  // Get first and last record info
  if (param == 1) defaultStartUsed = 0;
  if (param == 2) defaultStopUsed = 0;
  var needCheck = 1;
  if (!ImFirstRecord)
    {
    needCheck = 0;
    ImGetRecInfo(1);
    }
  if (!ImLastRecord)
    {
    needCheck = 0;
    ImGetRecInfo(-1);
    }
  if (needCheck)
    ImPatchCheck();
  }

var defaultStartUsed = 1;
var defaultStopUsed = 1;
function ImPatchCheck()
  {
  if (ImResetParams==2)
    {
    ImPatchInit(0);
    ImPatchSet(0);
    }
  else if (ImResetParams==3)
    {
    ImFirstRecord = null;
    ImLastRecord = null;
    ImResetParams = 0;
    ImPatchSet(0);
    }
  var isok = 1;
  var args = "im_patch";
  var ImLocOption;
// alert("ImPatchCheck, RecordCountNeeded="+RecordCountNeeded+", default start,stop="+defaultStartUsed+","+defaultStopUsed);
  if (RecordCountNeeded)
    {
    ImFirstRecord = null;
    ImLastRecord = null;
    if (defaultStartUsed) $("ImTStart").value = "NotSpecified";
    if (defaultStopUsed) $("ImTStop").value = "NotSpecified";
    ExportNewRS();
    return("");
    }

  if ( ($("ImTStart").value === "East Limb" && !$("ImEastLimb").checked) || ($("ImTStart").value.strip().empty()) )
     $("ImTStart").value = "NotSpecified";
  if ($("ImTStart").value == "NotSpecified")
    {
    if (ImFirstRecord)
      $("ImTStart").value = ImFirstRecord.keywords[0].values[0];
    }
  if ($("ImTStart").value == "NotSpecified")
    {
    $("ImTStart").style.backgroundColor=colorRed;
    isok = 0;
    }
  else
    {
    $("ImTStart").style.backgroundColor=colorWhite;
    if ($("ImEastLimb").checked)
      {
      $("ImTStart").value = "East Limb";
      args += ",t_start=" + "NotSpecified";
      defaultStartUsed = 0;
      }
    else
      args += ",t_start=" + $("ImTStart").value;
    }

  if ( ($("ImTStop").value === "West Limb" && !$("ImWestLimb").checked) || ($("ImTStop").value.strip().empty()) )
    $("ImTStop").value = "NotSpecified";
  if ($("ImTStop").value == "NotSpecified")
    {
    if (ImLastRecord)
      $("ImTStop").value = ImLastRecord.keywords[0].values[0];
    }
  if ($("ImTStop").value == "NotSpecified")
    {
    $("ImTStop").style.backgroundColor=colorRed;
    isok = 0;
    }
  else
    {
    $("ImTStop").style.backgroundColor=colorWhite;
    if ($("ImWestLimb").checked)
      {
      $("ImTStop").value = "West Limb";
      args += ",t_stop=" + "NotSpecified";
      defaultStopUsed = 0;
      }
    else
      args += ",t_stop=" + $("ImTStop").value;
    }

  if ($("ImTrack").checked) // checked for tracking, default
    {
    args += ",t=0";
    ImTracked = 1;
    }
  else
    {
    args += ",t=1";
    $("ImTRef").value = $("ImTStart").value;
    ImTracked = 0;
    }

  if ($("ImRegister").checked) // checked for registering by interpolation, default is no.
    args += ",r=1";
  else
    args += ",r=0";

  if ($("ImCrop").checked) // checked for crop to limb, default is no.
    args += ",c=1";
  else
    args += ",c=0";

  if ($("ImTDelta").value.strip().empty()) $("ImTDelta").value = "NotSpecified";
  if ($("ImTDelta").value == "NotSpecified" && ImFirstRecord)
    {
    var recset = $("ExportRecordSet").value;
    var posAt = recset.indexOf("@");
    if (posAt > 0)
      {
      var patt = /@[0-9]+[a-z]*/;
      var cads = $("ExportRecordSet").value.match(patt);
      $("ImTDelta").value = cads[0].substring(1);
      }
    else
      $("ImTDelta").value = ImFirstRecord.keywords[3].values[0] + "s";
    }
  if ($("ImTDelta").value == "NotSpecified")
    {
    $("ImTDelta").style.backgroundColor=colorWhite;
    }
  else
    {
    $("ImTDelta").style.backgroundColor=colorWhite;
    args += ",cadence=" + $("ImTDelta").value;
    }

  $("ImLocType").style.backgroundColor=noaaColor;
  ImLocOption = $("ImLocType").selectedIndex;
  // values are one of: stony, arcsec, pixels, carrlong
  args += ",locunits=" + $("ImLocType").options[ImLocOption].value;
  $("ImBoxType").style.backgroundColor=colorWhite;
  args += ",boxunits=" + $("ImBoxType").options[$("ImBoxType").selectedIndex].value;

  if ($("ImLocType").selectedIndex == 3 && ImTracked == 1)
    {
    $("ImCarrLi").style.display = "table-row";
    $("ImTRefLi").style.display = "none";
    $("ImTRef").value == "NotSpecified";
    }
  else
    {
    $("ImCarrLi").style.display = "none";
    $("ImTRefLi").style.display = "table-row";
    $("ImCarrot").value == "NotSpecified";
    }

  if ($("ImTRef").value.strip().empty()) $("ImTRef").value = "NotSpecified";
  if ($("ImTRef").value == "NotSpecified")
    {
    $("ImTRef").value = $("ImTStart").value;
    }
  if ($("ImTRef").value == "NotSpecified" && $("ImLocType").indexSelected != 3)
    {
    $("ImTRef").style.backgroundColor=colorRed;
    isok = 0;
    }
  else
    {
    $("ImTRef").style.backgroundColor=noaaColor;
    args += ",t_ref=" + $("ImTRef").value;
    }

  if ($("ImCarrot").value.strip().empty()) $("ImCarrot").value = "NotSpecified";
  if ($("ImLocType").selectedIndex == 3)
    {
    if ($("ImCarrot").value == "NotSpecified" && ImFirstRecord)
      $("ImCarrot").value = ImFirstRecord.keywords[1].values[0];
    if ($("ImCarrot").value == "NotSpecified" )
      {
      isok = 0;
      $("ImCarrot").style.backgroundColor=colorRed;
      }
    else
      {
      $("ImCarrot").style.backgroundColor=colorWhite;
      args += ",car_rot=" + $("ImCarrot").value;
      }
    }

  if ($("ImX").value.strip().empty()) $("ImX").value = "NotSpecified";
  if ($("ImLocType").selectedIndex == 3)
    {
    if ($("ImX").value == "NotSpecified" && ImFirstRecord)
      $("ImX").value = ImFirstRecord.keywords[2].values[0];
    }
  if ($("ImX").value == "NotSpecified")
    {
    isok = 0;
    $("ImX").style.backgroundColor=colorRed;
    }
  else
    {
    $("ImX").style.backgroundColor=noaaColor;
    args += ",x=" + $("ImX").value;
    }

  if ($("ImY").value.strip().empty()) $("ImY").value = "NotSpecified";
  if ($("ImY").value == "NotSpecified")
    {
    $("ImY").style.backgroundColor=colorRed;
    isok = 0;
    }
  else
    {
    $("ImY").style.backgroundColor=noaaColor;
    args += ",y=" + $("ImY").value;
    }

  if ($("ImWide").value.strip().empty()) $("ImWide").value = "NotSpecified";
  if ($("ImWide").value == "NotSpecified")
    {
    $("ImWide").style.backgroundColor=colorRed;
    isok = 0;
    }
  else
    {
    $("ImWide").style.backgroundColor=colorWhite;
    args += ",width=" + $("ImWide").value;
    }

  if ($("ImHigh").value.strip().empty()) $("ImHigh").value = "NotSpecified";
  if ($("ImHigh").value == "NotSpecified")
    {
    $("ImHigh").style.backgroundColor=colorRed;
    isok = 0;
    }
  else
    {
    $("ImHigh").style.backgroundColor=colorWhite;
    args += ",height=" + $("ImHigh").value;
    }

  if (isok)
    {
    $("ImVerify").innerHTML = "OK to submit";
    $("ImVerify").style.backgroundColor = colorWhite;
    }
  else
    {
    $("ImVerify").innerHTML = "Not Ready";
    $("ImVerify").style.backgroundColor = colorRed;
    }
  ExportProcessingOK = isok;
  CheckRediness();
  return (isok ? args : "");
  }

var ImFirstRecord = null;
var ImLastRecord = null;
function ImGetRecInfo(n)
  {
  // Get keywords for a single record.  n will be 1 for first record, -1 for last record.
  if (n==1)
    ImFirstRecord = null;
  else
    ImLastRecord = null;
  var timePrime = (firstTimePrime.length > 0 ? firstTimePrime : "T_REC");
// $("TESTMSG").innerHTML = timePrime;
  var keysneeded = timePrime+",CAR_ROT,CRLN_OBS,"+timePrime+"_step";
// alert("ImGetRecInfo("+n+") called");
  var recinfo;
  $("AjaxBusy").innerHTML = Ajax.activeRequestCount;
  var RecordSet = $("ExportRecordSet").value;
  $("ImRecordSet").innerHTML = RecordSet;
  new Ajax.Request('http://' + Host + '/cgi-bin/ajax/' + JSOC_INFO,
    {
    method: 'get',
    parameters: {"ds" : RecordSet, "op" : "rs_list", "n" : n, "key" : keysneeded },

    onSuccess: function(transport, json)
      {
      var thisN = ""+n;
      var response = transport.responseText || "no response text";
      var recinfo = response.evalJSON();

// $("TESTMSG").innerHTML += " thisN="+thisN;
      if (recinfo.status == 0)
        {
        if (thisN == "1")
          {
          ImFirstRecord = recinfo;
          if (ImLastRecord)
            ImPatchCheck();
          }
        if (thisN == "-1")
          {
          ImLastRecord = recinfo;
          if (ImFirstRecord)
            ImPatchCheck();
          }
        }
      else
        alert("ImPatch failed to get record info for n="+thisN+" of " + RecordSet);
      $("AjaxBusy").innerHTML = Ajax.activeRequestCount;
      },
    onFailure: function() { alert('Something went wrong...'); },
    onComplete: function() { $("AjaxBusy").innerHTML = Ajax.activeRequestCount; }
    });
  }

// Set display defaults

var ImResetParams = 1;
function ImPatchInit(isActive)
  {
  if (isActive == 0 && ImResetParams>0)
    {
    noaaColor = colorWhite;
    var requireColor = colorRed;
    if (ImResetParams==1)
        $("ProcessImPatch").style.display="none";
    $("ImTrack").checked = true;
    $("ImRegister").checked = false;
    $("ImCrop").checked = false;
    $("ImNOAA").style.backgroundColor=colorWhite; $("ImNOAA").value = "NotSpecified";
    $("ImTStart").value = "NotSpecified";
    $("ImTStop").value = "NotSpecified";
    $("ImTDelta").value = "NotSpecified";
    $("ImTRef").style.backgroundColor = requireColor; $("ImTRef").value = "NotSpecified";
    $("ImCarrot").style.backgroundColor = requireColor; $("ImCarrot").value = "NotSpecified";
    $("ImX").style.backgroundColor = requireColor; $("ImX").value = "NotSpecified";
    $("ImY").style.backgroundColor = requireColor; $("ImY").value = "NotSpecified";
    $("ImWide").style.backgroundColor = requireColor; $("ImWide").value = "NotSpecified";
    $("ImHigh").style.backgroundColor = requireColor; $("ImHigh").value = "NotSpecified";
    ImTracked = 1;
    ImFirstRecord = null;
    ImLastRecord = null;
    defaultStartUsed = 1;
    defaultStopUsed = 1;
    if (RecordCountNeeded == 1 || ImResetParams==2)
      ImResetParams = 0;
    }
  else
    {
    ImPatchSet(0);
    }
  }

// End of IM Patch code

//
// Maproj - map projections
//

// get NOAA AR info

var NOAA_TRef = "";
function MaprojGetNoaa()
  {
  if ($("MaprojNOAA").value.strip().empty()) $("MaprojNOAA").value = "NotSpecified";
  if ($("MaprojNOAA").value != "NotSpecified")
    {
    var noaaNum = 1 * $("MaprojNOAA").value;
    if (noaaNum < 7000) 
      {
      noaaNum = noaaNum + 10000; // OK for times after 1996 Jan.
      $("MaprojNOAA").value = noaaNum + "";
      }
    $("AjaxBusy").innerHTML = Ajax.activeRequestCount;
    new Ajax.Request('http://' + Host + '/cgi-bin/ajax/jsoc_info_jsoc2',
      {
      method: 'get',
      parameters: {"op" : "rs_list", "ds": "su_rsb.NOAA_ActiveRegions[][" + $("MaprojNOAA").value + "]", "key": "ObservationTime,LatitudeHG,LongitudeCM,LongitudeHG,LongitudinalExtent" },
      onSuccess: function(transport, json)
        {
        var response = transport.responseText || "no response text";
        var NOAA_rslist = response.evalJSON();
        if (NOAA_rslist.status == 0 && NOAA_rslist.count > 0) 
          {
          var minLong = 999, minRec=-1;
          var irec, nrecs = NOAA_rslist.count;
          var thisLong;
          for (irec=0; irec<nrecs; irec++)
            {
            thisLong = NOAA_rslist.keywords[2].values[irec];
            if (Math.abs(thisLong) < Math.abs(minLong))
              {
              minLong = thisLong;
              minRec = irec;
              }
            }
          if (minLong == 999)
            {
            $("MaprojNOAA").value = noaaNum + "Region not found"; 
            }
          else
            {
            NOAA_TRef =  NOAA_rslist.keywords[0].values[minRec];
            $("MaprojTRef").innerHTML = "at " + NOAA_TRef;
            $("MaprojX").value = NOAA_rslist.keywords[3].values[minRec];
            $("MaprojY").value = NOAA_rslist.keywords[1].values[minRec];
	    $("MaprojNOAA").style.backgroundColor=colorOptionSet;
            noaaColor = colorPreset;
            $("MaprojX").style.backgroundColor = noaaColor;
            $("MaprojY").style.backgroundColor = noaaColor;
	    // CheckMaproj();
            }
          }
        else
          {
          $("MaprojNOAA").value = noaaNum + "Region not found"; 
          }
        },
      onFailure: function()
        {
        alert('Something went wrong with NOAA num data request');
        $("MaprojNOAA").value = "Not Found";
        },
      onComplete: function() { $("AjaxBusy").innerHTML = Ajax.activeRequestCount; }
      });
    }
  }

function MaprojSet(param)
  {
  // Get first and last record info
  if (param == 1) defaultStartUsed = 0;
  if (param == 2) defaultStopUsed = 0;
  var needCheck = 1;
  if (!MaprojFirstRecord)
    {
    needCheck = 0;
    MaprojGetRecInfo(1);
    }
  if (!MaprojLastRecord)
    {
    needCheck = 0;
    MaprojGetRecInfo(-1);
    }
  if (needCheck)
    MaprojCheck();
  }

var defaultStartUsed = 1;
var defaultStopUsed = 1;
function MaprojCheck()
  {
  var isok = 1;
  var args = "Maproj";
  var MaprojLocOption;
// alert("MaprojCheck, RecordCountNeeded="+RecordCountNeeded+", default start,stop="+defaultStartUsed+","+defaultStopUsed);
  if (RecordCountNeeded)
    {
    ExportNewRS();
    return("");
    }

  args += ",map=" + $("MaprojProj").value;
  if ($("MaprojProj").selectedIndex == 1)
    {
    $("MaprojReflatLine").style.display="table-row";
    }
  else
    {
    $("MaprojReflatLine").style.display="none";
    }
 
  if ($("MaprojGrid").value != "none")
    {
    args += ",grid=" + $("MaprojGrid").value;
    }

  if (MaprojFirstRecord)
    {
    // var DegPerArcsec = (180.0/(696.0*3.14159))*(149640.0/(180.0*3600.0/3.14159)  // deg/Mm * Mm/arcsec 
    var DegPerArcsec = 215/3600.0;  // deg/RsunRadian  * AURadian/arcsec
    var DegPerPixel = new Number(MaprojFirstRecord.keywords[4].values[0] * DegPerArcsec);
    $("MaprojMaxScale").innerHTML = DegPerPixel.toPrecision(3);
    $("MaprojMaxScale").style.backgroundColor=colorWhite;
    if ($("MaprojScale").value == "NotSpecified") $("MaprojScale").value = $("MaprojMaxScale").innerHTML;
    }

  if ($("MaprojX").value.strip().empty()) $("MaprojX").value = "NotSpecified";
  if ($("MaprojX").value == "NotSpecified")
    {
    isok = 0;
    $("MaprojX").style.backgroundColor=colorRed;
    }
  else
    {
    $("MaprojX").style.backgroundColor=noaaColor;
    args += ",clon=" + $("MaprojX").value;
    }

  if ($("MaprojY").value.strip().empty()) $("MaprojY").value = "NotSpecified";
  if ($("MaprojY").value == "NotSpecified")
    {
    $("MaprojY").style.backgroundColor=colorRed;
    isok = 0;
    }
  else
    {
    $("MaprojY").style.backgroundColor=noaaColor;
    args += ",clat=" + $("MaprojY").value;
    }

  if ($("MaprojScale").value.strip().empty()) $("MaprojScale").value = "NotSpecified";
  if ($("MaprojScale").value == "NotSpecified")
    {
    $("MaprojScale").style.backgroundColor=colorRed;
    isok = 0;
    }
  else
    {
    $("MaprojScale").style.backgroundColor=colorWhite;
    args += ",scale=" + $("MaprojScale").value;
    }

  if ($("MaprojWide").value.strip().empty()) $("MaprojWide").value = "NotSpecified";
  if ($("MaprojWide").value == "NotSpecified")
    {
    $("MaprojWide").style.backgroundColor=colorRed;
    isok = 0;
    }
  else
    {
    $("MaprojWide").style.backgroundColor=colorWhite;
    args += ",cols=" + $("MaprojWide").value;
    }

  if ($("MaprojHigh").value.strip().empty())
    $("MaprojHigh").value = "NotSpecified";
  if ($("MaprojHigh").value == "NotSpecified")
    {
    $("MaprojHigh").style.backgroundColor=colorRed;
    isok = 0;
    }
  else
    {
    $("MaprojHigh").style.backgroundColor=colorWhite;
    args += ",rows=" + $("MaprojHigh").value;
    }

  if (isok)
    {
    $("MaprojVerify").innerHTML = "OK to submit";
    $("MaprojVerify").style.backgroundColor = colorWhite;
    }
  else
    {
    $("MaprojVerify").innerHTML = "Not Ready";
    $("MaprojVerify").style.backgroundColor = colorRed;
    }
  ExportProcessingOK = isok;
  CheckRediness();
  return (isok ? args : "");
  }

var MaprojFirstRecord = null;
var MaprojLastRecord = null;
function MaprojGetRecInfo(n)
  {
  // Get keywords for a single record.  n will be 1 for first record, -1 for last record.
  if (n==1)
    MaprojFirstRecord = null;
  else
    MaprojLastRecord = null;
  var timePrime = (firstTimePrime.length > 0 ? firstTimePrime : "T_REC");
// $("TESTMSG").innerHTML = timePrime;
  var keysneeded = timePrime+",CAR_ROT,CRLN_OBS,"+timePrime+"_step,CDELT1,CTYPE1";
// alert("MaprojGetRecInfo("+n+") called");
  var recinfo;
  $("AjaxBusy").innerHTML = Ajax.activeRequestCount;
  var RecordSet = $("ExportRecordSet").value;
  new Ajax.Request('http://' + Host + '/cgi-bin/ajax/' + JSOC_INFO,
    {
    method: 'get',
    parameters: {"ds" : RecordSet, "op" : "rs_list", "n" : n, "key" : keysneeded },

    onSuccess: function(transport, json)
      {
      var thisN = ""+n;
      var response = transport.responseText || "no response text";
      var recinfo = response.evalJSON();

// $("TESTMSG").innerHTML += " thisN="+thisN;
      if (recinfo.status == 0)
        {
        if (thisN == "1")
          {
          MaprojFirstRecord = recinfo;
          if (MaprojLastRecord)
            MaprojCheck();
          }
        if (thisN == "-1")
          {
          MaprojLastRecord = recinfo;
          if (MaprojFirstRecord)
            MaprojCheck();
          }
        }
      else
        alert("Maproj failed to get record info for n="+thisN+" of " + RecordSet);
      $("AjaxBusy").innerHTML = Ajax.activeRequestCount;
      },
    onFailure: function() { alert('Something went wrong...'); },
    onComplete: function() { $("AjaxBusy").innerHTML = Ajax.activeRequestCount; }
    });
  }

// Set display defaults

function MaprojInit(isActive)
  {
  if (!isActive)
    {
    noaaColor = colorWhite;
    var requireColor = colorRed;
    $("ProcessMaproj").style.display="none";
    $("MaprojReflatLine").style.display="none";
    $("MaprojNOAA").style.backgroundColor=colorWhite; $("MaprojNOAA").value = "NotSpecified";
    $("MaprojScale").style.backgroundColor = requireColor; $("MaprojScale").value = "NotSpecified";
    $("MaprojX").style.backgroundColor = requireColor; $("MaprojX").value = "NotSpecified";
    $("MaprojY").style.backgroundColor = requireColor; $("MaprojY").value = "NotSpecified";
    $("MaprojWide").style.backgroundColor = requireColor; $("MaprojWide").value = "NotSpecified";
    $("MaprojHigh").style.backgroundColor = requireColor; $("MaprojHigh").value = "NotSpecified";
    MaprojFirstRecord = null;
    MaprojLastRecord = null;
    defaultStartUsed = 1;
    defaultStopUsed = 1;
    }
  else
    {
    MaprojSet(0);
    }
  }

// End of MapProj code

//
// Processing details for all options
//

function ProcessingInit()
  {
  var ProcessingOptionsHTML = "";
  ExportProcessingOptions = new Array();
  var ExpOpt;
  var iOpt;

  // Add Select options for each processing type, in order that they may be used.
  iOpt = 0;
  ProcessingOptionsHTML += 
    '<input type="checkbox" checked="true" value="no_op" id="OptionNone" onChange="SetProcessing('+iOpt+');" />' +
    'no_op - none&nbsp;' +
    '<input id="ProcessingCheckboxHide" type="checkbox" checked="false" onChange="ProcessingEnabled();" />&nbsp;hide<br>';
  ExpOpt = new Object();
  ExpOpt.id="OptionNone";
  ExpOpt.rowid = "ExpSel_none";
  ExportProcessingOptions[iOpt] = ExpOpt;

  iOpt++;
  ProcessingOptionsHTML +=
    '<input type="checkbox" checked="false" value="aia_scale" id="OptionAiaScale" onChange="SetProcessing('+iOpt+');" /> ' +
    'aia_scale - Scale AIA lev1 to 0.6 arcsec/pixel<br>';
  ExpOpt = new Object();
  ExpOpt.id = "OptionAiaScale";
  ExpOpt.rowid = "ProcessAiaScale";
  ExpOpt.Init = AiaScaleInit;
  ExpOpt.Check = AiaScaleCheck;
  ExpOpt.Set = AiaScaleSet;
  ExportProcessingOptions[iOpt] = ExpOpt;

  iOpt++;
  ProcessingOptionsHTML += 
    '<input type="checkbox" checked="false" value="resize" id="OptionResize" onChange="SetProcessing('+iOpt+');" /> ' +
    'resize - Resize and rotate if needed, use sub-pixel registration<br>';
  ExpOpt = new Object();
  ExpOpt.id = "OptionResize";
  ExpOpt.rowid = "ProcessResize";
  ExpOpt.Init = ResizeInit;
  ExpOpt.Check = ResizeCheck;
  ExpOpt.Set = ResizeSet;
  ExportProcessingOptions[iOpt] = ExpOpt;

  iOpt++;
  ProcessingOptionsHTML += 
    '<input type="checkbox" checked="false" value="rebin" id="OptionRebin" onChange="SetProcessing('+iOpt+');" /> ' +
    'rebin - Rebin with boxcar or gaussian smoothing<br>';
  ExpOpt = new Object();
  ExpOpt.id = "OptionRebin";
  ExpOpt.rowid = "ProcessRebin";
  ExpOpt.Init = RebinInit;
  ExpOpt.Check = RebinCheck;
  ExpOpt.Set = RebinSet;
  ExportProcessingOptions[iOpt] = ExpOpt;

  iOpt++;
  ProcessingOptionsHTML += 
    '<input type="checkbox" checked="false" value="im_patch" id="OptionImPatch" onChange="SetProcessing('+iOpt+');" /> ' +
    'im_patch - Extract sub-frame<br>';
  ExpOpt = new Object();
  ExpOpt.id = "OptionImPatch";
  ExpOpt.rowid = "ProcessImPatch";
  ExpOpt.Init = ImPatchInit;
  ExpOpt.Check = ImPatchCheck;
  ExpOpt.Set = null;
  ExportProcessingOptions[iOpt] = ExpOpt;

  iOpt++;
  ProcessingOptionsHTML += 
    '<input type="checkbox" checked="false" value="maproj" id="OptionMaproj" onChange="SetProcessing('+iOpt+');" /> ' +
    'maproj - Extract a sub-frame and remap to a chosen projection.<br>';
  ExpOpt = new Object();
  ExpOpt.id = "OptionMaproj";
  ExpOpt.rowid = "ProcessMaproj";
  ExpOpt.Init = MaprojInit;
  ExpOpt.Check = MaprojCheck;
  ExpOpt.Set = null;
  ExportProcessingOptions[iOpt] = ExpOpt;

  $("ExportProcessing").innerHTML = ProcessingOptionsHTML;

  var nOpt = iOpt + 1;

for (iOpt=1; iOpt<nOpt; iOpt++)
    {
    var id;
    ExpOpt = ExportProcessingOptions[iOpt];
    ExpOpt.Init(0);
    $(ExpOpt.id).checked = false;
    $(ExpOpt.rowid).style.display = "none";
    }
  }

