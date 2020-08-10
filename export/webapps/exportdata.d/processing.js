// Functions for export options
//
// CheckXXXX creates the contents of ExportProcessingArgs from parameters for
// this export type
//
// Each processing type needs a common name, e.g. ImPatch, used for the first part of
// function and variable names associated with that processing type.
// For processing 'Xxxx' the following functions should be defined in this file.
//
//  XxxxInit()  - reset all visible elements to default state, reset associated global variables.
//                       XxxxInit() should be called from ProcessingInit() at the end of this file to initialize
//                       the html for this processing type.
//  XxxxCheck() - Examine all user set parameters for correctness and consistency, mark incorrect entries with color.
//                Build string to use in call to jsoc_fetch in the "process" parameter for this processing type.
//                this string will start with the process name, e.g. "hg_patch" and contains all needed command line
//                params with comma delimiters.  Flags should be e.g. "c=1" rather than "-c".
//                Returns the string if all is OK, or "" if user action is still required.
//                The Check function must update the SizeRatio global var with the size ratio resulting
//                from the chosen processing.
//  XxxxSet(xx) - Optional, may be used in the onChange or onClick events to accept user input,
//                or the XxxxCheckX function may be used for this purpose.
//  Other functions may be defined as needed and may be invoked in the XxxxSet() function where this processing type
//                is tested for completion, but only the XxxxInit and CheckXxxx funtions will be used elsewhere.
//
//  Finally, make entries in ProcessingInit() at the end of this file.
//

// Global vars for multiple processing options

var processingFirstRecord = null;
var processingLastRecord = null;
var defaultStartUsed = 1;
var defaultStopUsed = 1;
var keyTIME = 0;
var keySTEP = 1;
var keyCARROT = 2;
var keyCRLN = 3;
var keyCRLT = 4;
var keyCDELT = 5;
var keyCTYPE = 6;
var segDims = 2;

var keyNoaaTime = 0;
var keyNoaaLat = 1;
var keyNoaaLon = 2;
var keyNoaaLonX = 3;

// Support functions for multiple processing options
// This function gets the first or last record of a series.
// It gets called ONLY from the Set() functions. So, the argument to DoCheck() should be false.
function processingGetRecInfo(DoCheck, exportOption, n)
  {
    var parameters = null;
    var argString = null;
    var url = null;

  // Get keywords for a single record.  n will be 1 for first record, -1 for last record.
  if (n==1)
    processingFirstRecord = null;
  else
    processingLastRecord = null;
  var timePrime = (firstTimePrime.length > 0 ? firstTimePrime : "T_REC");
  var keysneeded = timePrime+","+timePrime+"_step,CAR_ROT,CRLN_OBS,CRLT_OBS,CDELT1,CTYPE1";
  var segsneeded = firstRealSegment;
  var recinfo;
  var RecordSet = $("ExportRecordSet").value;
  $("ImRecordSet").innerHTML = RecordSet;

    parameters = { "ds" : RecordSet, "op" : "rs_list", "n" : n, "key" : keysneeded, "seg" : segsneeded, "l" : 1 };

  var ajaxParameters = {
    method: 'get',
    parameters: parameters,

    onSuccess: function(transport, json)
      {
          var thisN = ""+n;
          var response = transport.responseText || "no response text";
          var recinfo = response.evalJSON();

          if (recinfo.status == 0)
            {
            if (thisN == "1")
              {
              processingFirstRecord = recinfo;

              }
            if (thisN == "-1")
              {
              processingLastRecord = recinfo;
              }

                if (exportOption)
                {
                    exportOption.argsReady = false;
                    exportOption.paramsValid = null;
                    exportOption.Check(false);
                }
            }
          else
          {
            alert("Processing setup code failed to get record info for n="+thisN+" of " + RecordSet);
          }

          if (exportOption)
          {
              exportOption.argsReady = true;
          }
      },
    onFailure: function() { alert('[ processingGetRecInfo ] Something went wrong...'); if (exportOption) { exportOption.argsReady = true; } },
    onComplete: function() { if (exportOption) { exportOption.argsReady = true; } }
    };

    argString = Object.keys(parameters).map(function(key) { return key + '=' + encodeURIComponent(parameters[key]); }).join('&');
    url = 'http://' + Host + '/cgi-bin/ajax/' + JSOC_INFO + '?' + argString;

  new Ajax.Request('http://' + Host + '/cgi-bin/ajax/' + JSOC_INFO, ajaxParameters);
  }

// Processing Options Code

//
// Processing details for HmiB2ptr --> HMI B to phi,theta,r
//

var HmiB2ptrOption;

function HmiB2ptrSet(param)
  {
    var checkRes = null;

    this.argsReady = false;
    this.paramsValid = null;

    checkRes = this.Check(false);

    if (checkRes === null)
    {
        this.argsReady = true;
        return 1;
    }
    else if (checkRes.length == 0)
    {
        this.argsReady = true;
        return 0;
    }

    this.argsReady = true;
    return 0;
  }

function HmiB2ptrInit(onLoad)
  {
    if (onLoad)
    {
        $("ProcessHmiB2ptr").style.display="none";
    }
  }

function HmiB2ptrCheck(fromSetProcessing)
  {
  var HmiB2ptrSizeRatio = 1.0;
  var isok = true;
  var args = "HmiB2ptr,l=1";

  if (!$("OptionHmiB2ptr").checked)
  {
    // If this processing option is not selected, then do not check parameter values.
    return '';
  }

    if (fromSetProcessing)
    {
        if (!this.argsReady)
        {
            // a Set() is pending
            return null;
        }
    }

    if (this.paramsValid !== null)
    {
        // Already checked.
        return this.paramsValid; // Empty string if the previous call determined the arguments to be invalid.
    }

    if (SeriesName.toUpperCase() === 'hmi.B_720s'.toUpperCase())
    {
        // $("ExportFilenameFmt").value = $("ExportFilenameFmt").value.replace('T_REC','T_OBS');
    }
    else
    {
        isok = false;
        alert("Error - HmiB2ptr can only be used for series hmi.B_720s");
        return 'error';
    }

  ExportProcessingOptions[HmiB2ptrOption].Size = HmiB2ptrSizeRatio;
  this.paramsValid = (isok ? args : "");
  CheckRediness();
  return (isok ? args : "");
  }

// End HmiB2ptr

//
// Processing details for AIA normalized Scaling
//

var AiaScaleOption;

function AiaScaleSet(param)
  {
    var checkRes = null;

    this.argsReady = false;
    this.paramsValid = null;

    if (param == 'usempt')
    {
        // if the user checked the use-mpt box, then display the mpt row element, else hide it
        if ($('AiaScaleUseMptCheckbox').checked)
        {
            $('AiaScaleMptRow').style.display = 'table-row';
        }
        else
        {
            $('AiaScaleMptRow').style.display = 'none';
        }
    }
    else if (param == 'mpt')
    {
        // nothing to do
    }
    else if (param = 'docutout')
    {
        // if the user checked the perform-cutout box, then display the perform-cutout tbody element, else hide it
        if ($('AiaScalePerformCutoutCheckbox').checked)
        {
            $$('tr.AiaScaleCutoutBody').each(function(elem) { elem.style.display = 'table-row'; });
        }
        else
        {
            $$('tr.AiaScaleCutoutBody').each(function(elem) { elem.style.display = 'none'; });
        }
    }
    else if (param == 'xc')
    {
        val = parseInt($('AiaScaleXc').value)
        if (val >= 4096)
        {
            $('AiaScaleXc').value = '4095'
        }
        else if (val <= -4096)
        {
            $('AiaScaleXc').value = '-4095'
        }
    }
    else if (param == 'yc')
    {
        val = parseInt($('AiaScaleYc').value)
        if (val >= 4096)
        {
            $('AiaScaleYc').value = '4095'
        }
        else if (val <= -4096)
        {
            $('AiaScaleYc').value = '-4095'
        }
    }
    else if (param == 'wide')
    {
        val = parseInt($('AiaScaleWide').value)
        if (val < 0)
        {
            $('AiaScaleWide').value = '0'
        }
        else if (val > 4096)
        {
            $('AiaScaleWide').value = '4096'
        }
    }
    else if (param == 'high')
    {
        val = parseInt($('AiaScaleHigh').value)
        if (val < 0)
        {
            $('AiaScaleHigh').value = '0'
        }
        else if (val > 4096)
        {
            $('AiaScaleHigh').value = '4096'
        }
    }
    else
    {
        alert('Error - invalid AiaScaleSet() argument ' + param)
        return 1;
    }

    checkRes = this.Check(false);

    if (checkRes === null)
    {
        // Error - uncheck the processing step.
        this.argsReady = true;
        return 1;
    }
    else if (checkRes.length == 0)
    {
        // The arguments are invalid, but do not uncheck the processing step.
        this.argsReady = true;
        return 0;
    }

    this.argsReady = true;
    return 0;
  }

function AiaScaleInit(onLoad)
{
    if (onLoad)
    {
        $("ProcessAiaScale").style.display = "none";
        $('AiaScaleUseMptCheckbox').checked = false;
        $('AiaScaleMptRow').style.display = 'none';
        $('AiaScaleMptSelect').value = 'aia.master_pointing3h';
        $('AiaScalePerformCutoutCheckbox').checked = false
        // $$('tr.AiaScaleCutoutBody').each(function(elem) { elem.setStyle({ 'display': 'none' }); })
        $$('tr.AiaScaleCutoutBody').each(function(elem) { elem.style.display = 'none'; })
        $('AiaScaleCutoutXc').value = '0.0';
        $('AiaScaleCutoutYc').value = '0.0';
        $('AiaScaleCutoutWide').value = '4096';
        $('AiaScaleCutoutHigh').value = '4096';
    }
}

function AiaScaleCheck(fromSetProcessing)
{
    var AiaScaleSizeRatio = 1.0;
    var isok = true;
    var args = null;
    var compatibleFormat = null;

    if (!$("OptionAiaScale").checked)
    {
        // If this processing option is not selected, then do not check parameter values.
        return '';
    }

    if (fromSetProcessing)
    {
        if (!this.argsReady)
        {
            // a Set() is pending
            return null;
        }
    }

    if (this.paramsValid !== null)
    {
        // already checked; even though we use a different set of arguments for cut-out vs. non-cut-out, a change
        // from cut-out to non-cut-out and vice versa always results in the Set() function to be called, which blows
        // away the cached arguments
        return this.paramsValid; // Empty string if the previous call determined the arguments to be invalid.
    }

    // if the series is aia.lev1 and we are producing non-cut-outs, then data go into aia.lev1p5 (indicated by
    // 'aia_scale' processing); otherwise, data go into SeriesName + '_mod' (indicated by 'aia_scale_mod' processing)

    if (IsAiaLev1)
    {
        // hack away! what we should do is read the jsoc.export_procs out field to determine what the output series is going to be;
        // however, it is time to move on; simply use the substitution s/lev1/lev1p5/ and if we end up with aia.lev1p5,
        // then change the filename format to one compatible with aia.lev1p5, and use the aia_scale_orig proc
        outputSeries = SeriesName.strip().replace(/lev1/i, 'lev1p5');

        if (outputSeries == AIA_LEV1P5)
        {
            args = 'aia_scale_orig';

            // the output series is aia.lev1p5, so let's use a filename format compatible with aia.lev1p5;
            // actually, we could do this for all processing since the filename format must always be compatible
            // with the output series, not the input series

            // save the original in case the user un-checks aia_scale
            $('ExportFilenameFmt').store({ 'originalFormat' : $('ExportFilenameFmt').value });
            compatibleFormat = GetCompatibleFormat(AIA_LEV1P5, AiaLev1AttributesGlobal, AiaLev1KeywordsGlobal);

            // now, replace aia.lev1p5 with the original series name, and remove the RequestID specification (output series
            // have this extra prime-key keyword)
            $("ExportFilenameFmt").value = compatibleFormat.replace(AIA_LEV1P5, SeriesName.strip()).replace(/[.]{requestid}/i, '');
        }
        else
        {
            args = 'aia_scale_aialev1';
        }
    }
    else
    {
        args = 'aia_scale_other';
    }

    // choose between cut-out and non-cut-out
    if (!$('AiaScalePerformCutoutCheckbox').checked)
    {
        // because the prime keys of the aia series are not consistent, there is no way to use aia.lev1p5 as the output series
        // for all aia-lev1 series
        if (0 && SeriesName.strip().toLowerCase().search(/aia[.]lev1/) == 0)
        {
            // scale aia-lev1 series, save results in aia.lev1p5
            args = 'aia_scale';
        }
    }
    else
    {
        // because the prime keys of the aia series are not consistent, there is no way to use aia.lev1p5 as the output series
        // for all aia-lev1 series
        if (0 && SeriesName.strip().toLowerCase().search(/aia[.]lev1/) == 0)
        {
            // scale aia-lev1 series and perform cut-out, use the save results in aia.lev1p5
            args = 'aia_scale';
        }

        val = parseInt($('AiaScaleCutoutXc').value)
        if (Math.abs(val) > 4096)
        {
            isok = false;
        }
        else
        {
            if (val != 0)
            {
                args = args + ',' + 'xc=' + $('AiaScaleCutoutXc').value;
            }

            val = parseInt($('AiaScaleCutoutYc').value);
            if (Math.abs(val) > 4096)
            {
                isok = false;
            }
            else
            {
                if (val != 0)
                {
                    args = args + ',' + 'yc=' + $('AiaScaleCutoutYc').value;
                }

                val = parseInt($('AiaScaleCutoutWide').value);
                if (val < 0 || val > 4096)
                {
                    isok = false;
                }
                else
                {
                    if (val != 4096)
                    {
                        args = args + ',' + 'wide=' + $('AiaScaleCutoutWide').value;
                    }

                    val = parseInt($('AiaScaleCutoutHigh').value);
                    if (val < 0 || val > 4096)
                    {
                        isok = false;
                    }
                    else
                    {
                        args = args + ',' + 'high=' + $('AiaScaleCutoutHigh').value;
                    }
                }
            }
        }
    }

    // provide MPT series if user has selected MPT
    if ($('AiaScaleUseMptCheckbox').checked)
    {
        args = args + ',' + 'mpt=' + $('AiaScaleMptSelect').value;
    }

    ExportProcessingOptions[AiaScaleOption].Size = AiaScaleSizeRatio;
    this.paramsValid = (isok ? args : "");
    CheckRediness();
    return (isok ? args : "");
}

// End AiaScale

//
// Processing details for rebin using jsoc_rebin
//

var RebinOption;

function RebinInit(onLoad)
  {
    if (onLoad)
    {
        // Clear out any old settings.
        $("RebinMethod").selectedIndex=0;
        $("RebinSegments").checked = false;
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
    }

    // Display the rebin options table.
  }

function RebinSet(control)
  {
    var checkRes = null;

    this.argsReady = false;
    this.paramsValid = null;

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

    checkRes = this.Check(false);

    if (checkRes === null)
    {
        this.argsReady = true;
        return 1;
    }
    else if (checkRes.length == 0)
    {
        this.argsReady = true;
        return 0;
    }

    this.argsReady = true;
    return 0;
  }

function RebinCheck(fromSetProcessing)
  {
  var RebinSizeRatio = 1.0;
  var rv = "rebin";
  var isok = true;

    if (!$("OptionRebin").checked)
    {
        // If this processing option is not selected, then do not check parameter values.
        return '';
    }

    if (fromSetProcessing)
    {
        if (!this.argsReady)
        {
            // a Set() is pending
            return null;
        }
    }

    if (this.paramsValid !== null)
    {
        // Already checked.
        return this.paramsValid; // Empty string if the previous call determined the arguments to be invalid.
    }

  if ($("RebinSegments").checked)
    {
    rv = rv + ",A=1";
    }

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
      isok = false;
    rv = rv + ",scale=" + $("RebinScale").value;
    RebinSizeRatio = parseFloat($("RebinScale").value);
    }

  if ($("RebinMethod").selectedIndex == 1)
    {
    rv = rv + ",method=gaussian";
    if (parseFloat($("RebinFWHM").value) <= 0)
      {
      isok = false;
      }
    else
      {
      rv = rv + ",FWHM=" + $("RebinFWHM").value;
      }

    if (parseFloat($("RebinNvector").value) == -1.0)
      {
      isok = false;
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

  ExportProcessingOptions[RebinOption].Size = RebinSizeRatio*RebinSizeRatio;
  this.paramsValid = (isok ? rv : "");
  CheckRediness();
  return (isok ? rv : "");
  }

// End Rebin

//
// Processing details for rebin using jsoc_resize
//

var ResizeOption;

function ResizeInit(onLoad)
  {
    if (onLoad)
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
    var checkRes = null;

    this.argsReady = false;
    this.paramsValid = null;

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

    processingGetRecInfo(this.Check, this, 1);
    return 0;
  }

// If the Check() functions are called without first calling the corresponding Set() functions, then
// it is possible that the Check() functions are operating on stale data because the Check() functions
// may use variables set by the Set() functions.
//
// We basically do not want to call the Check() functions from SetProcessing() until the the Set()
// functions have completed. The way to do that is to wait until this.paramsValid is !== null.
//
// The parameter, if true, says we are calling from SetProcessing(). If that is the case, then
// return this.paramsValid. SetProcessing() keeps calling the Check() function until
// this.paramsValid !== null.
function ResizeCheck(fromSetProcessing)
  {
  var ResizeSizeRatio = 1.0;
  var rv = "resize";
  var isok = true;

    if (!$("OptionResize").checked)
    {
        // If this processing option is not selected, then do not check parameter values.
        return '';
    }

    if (fromSetProcessing)
    {
        if (!this.argsReady)
        {
            // a Set() is pending
            return null;
        }
    }

    if (this.paramsValid !== null)
    {
        // Already checked.
        return this.paramsValid; // Empty string if the previous call determined the arguments to be invalid.
    }

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
      isok = false;
    if (processingFirstRecord)
      ResizeSizeRatio = processingFirstRecord.keywords[keyCDELT].values[0] / parseFloat($("ResizeCdelt").value);
    }
  else
    rv = rv + ",rescale=0";

  if ($("ResizeCrop").checked)
    {
    rv = rv + ",c=1";
    }

  ExportProcessingOptions[ResizeOption].Size = ResizeSizeRatio;
  this.paramsValid = (isok ? rv : "");
  CheckRediness();
  return (isok ? rv : "");
  }

// End Resize

//
// Process ImPatch
//

var ImPatchOption;
var noaaColor;
var ImTracked = 0;
var newDimX = 0;
var newDimY = 0;

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
          thisTime = NOAA_rslist.keywords[keyNoaaTime].values[irec];
          thisLat = NOAA_rslist.keywords[keyNoaaLat].values[irec];
          thisLong = NOAA_rslist.keywords[keyNoaaLon].values[irec];
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
      onComplete: function() { }
      });
    }
  }

// The values of the ImPatch inputs have changed. We need to call ImPatchCheck() to make sure the changes were valid.
// If there is no first or last record saved in global variables, then we need to get them. In this case, we have
// to put ImPatchCheck() in a callback, since we cannot actually do the check synchronously - we have to wait until
// a jsoc_info call has completed.
//
// If the processingFirstRecord or processingLastRecord is missing (no record-set query has been entered, or it hasn't
// been processed yet) OR if the record-set specification has changed, then the first and last record information will
// be fetched asynchronously. Then ImPatchCheck() will be run. If no first/last record information is fetched, then only
// ImPatchCheck() is run.
function ImPatchSet(param)
  {
    var checkRes = null;

    // Get first and last record info
    this.argsReady = false;
    this.paramsValid = null;
    if (param == 1) defaultStartUsed = 0;
    if (param == 2) defaultStopUsed = 0;

    if ($("ImRecordSet").innerHTML !== $("ExportRecordSet").value)
    {
        processingFirstRecord = null;
        processingLastRecord = null;
        $("ImTDelta").value = "NotSpecified";
    }

    // RecordCount cannot be null - you have to have a valid record-set before you can click on im_patch processing.
    if (RecordCount === null)
    {
        processingFirstRecord = null;
        processingLastRecord = null;
        if (defaultStartUsed) $("ImTStart").value = "NotSpecified";
        if (defaultStopUsed) $("ImTStop").value = "NotSpecified";
        ExportNewRS(); // This will gracefully fail if there is no record-set specification entered.
    }

    if ($("ImTrack").checked)
    {
        ImTracked = 1;
    }
    else
    {
        ImTracked = 0;
    }

    if (processingFirstRecord && processingLastRecord)
    {
        checkRes = this.Check(false);

        if (checkRes === null)
        {
            this.argsReady = true;
            return 1;
        }
        else if (checkRes.length == 0)
        {
            this.argsReady = true;
            return 0;
        }

        this.argsReady = true;
        return 0;
    }

    // This looks like it will call potentially ImPatchCheck() twice, which is not ideal.
    if (!processingFirstRecord)
    {
        processingGetRecInfo(this.Check, this, 1);
    }
    if (!processingLastRecord)
    {
        processingGetRecInfo(this.Check, this, -1);
    }

    return 0;
  }

// ImResetParams values:
//  0 set in ImPatchCheck (do not reset parameters)
//  1 set onload (reset parameters to default)
//  2 set by Reset params button (reset parameters to default)
//  3 set by Update RecordSet Times button (reset parameters to default)

function ImPatchCheck(fromSetProcessing)
  {
    var args = "im_patch";
    var isok = true;
    var localImPatchSizeRatio = 1.0;

    if (!$("OptionImPatch").checked)
    {
        // If this processing option is not selected, then do not check parameter values.
        return '';
    }

    // If this.argsReady == false, we know that some Set() call was made and we are possibly waiting for
    // resolution of asynchronous calls. If this Check() was called by Set(), then we need to continue
    // and finalize the Set() call. But if this Check() was called by SetProcessing(), then we need to
    // check if there is a pending Set() call. We do that by looking at this.argsReady. If it is
    // false, then we know that there is a pending Set(), and SetProcessing() must wait. If it is true,
    // then there is no pending Set(), and SetProcessing() can get an answer and continue.
    if (fromSetProcessing)
    {
        if (!this.argsReady)
        {
            // a Set() is pending
            return null;
        }
    }

    if (this.paramsValid !== null)
    {
        // Already checked.
        return this.paramsValid; // Empty string if the previous call determined the arguments to be invalid.
    }

  if (RecordCount === null)
    {
        return("");
    }

  if ( ($("ImTStart").value === "East Limb" && !$("ImEastLimb").checked) || ($("ImTStart").value.strip().empty()) )
     $("ImTStart").value = "NotSpecified";
  if ($("ImTStart").value == "NotSpecified")
    {
    if (processingFirstRecord)
      $("ImTStart").value = processingFirstRecord.keywords[keyTIME].values[0];
    }
  if ($("ImTStart").value == "NotSpecified")
    {
    $("ImTStart").style.backgroundColor=colorRed;
    isok = false;
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
    if (processingLastRecord)
      $("ImTStop").value = processingLastRecord.keywords[keyTIME].values[0];
    }
  if ($("ImTStop").value == "NotSpecified")
    {
    $("ImTStop").style.backgroundColor=colorRed;
    isok = false;
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
    }
  else
    {
    args += ",t=1";
    $("ImTRef").value = $("ImTStart").value;
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
  if ($("ImTDelta").value == "NotSpecified" && processingFirstRecord)
    {
    var recset = $("ExportRecordSet").value; // This can be undefined if user has changed record set and has not waited for the
                                             // record count to complete.
    var posAt = recset.indexOf("@");
    if (posAt > 0)
      {
      var patt = /@[0-9]+[a-z]*/;
      var cads = $("ExportRecordSet").value.match(patt);
      $("ImTDelta").value = cads[0].substring(1);
      }
    else
      $("ImTDelta").value = processingFirstRecord.keywords[keySTEP].values[0] + "s";
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
  // values are one of: stony, arcsec, pixels, carrlong
  args += ",locunits=" + $("ImLocType").options[$("ImLocType").selectedIndex].value;
  $("ImBoxType").style.backgroundColor=colorWhite;
  args += ",boxunits=" + $("ImBoxType").options[$("ImBoxType").selectedIndex].value;

  if ($("ImLocType").selectedIndex == 3 && ImTracked == 1)
    {
    $("ImCarrLo").style.display = "table-row";
    $("ImTRefLi").style.display = "none";
    $("ImTRef").value == "NotSpecified";
    }
  else
    {
    $("ImCarrLo").style.display = "none";
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
    isok = false;
    }
  else
    {
    $("ImTRef").style.backgroundColor=noaaColor;
    args += ",t_ref=" + $("ImTRef").value;
    }

  if ($("ImCarrot").value.strip().empty()) $("ImCarrot").value = "NotSpecified";
  if ($("ImLocType").selectedIndex == 3)
    {
    if ($("ImCarrot").value == "NotSpecified" && processingFirstRecord)
      $("ImCarrot").value = processingFirstRecord.keywords[keyCARROT].values[0];
    if ($("ImCarrot").value == "NotSpecified" )
      {
      isok = false;
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
    if ($("ImX").value == "NotSpecified" && processingFirstRecord)
      $("ImX").value = processingFirstRecord.keywords[keyCRLN].values[0];
    }
  if ($("ImX").value == "NotSpecified")
    {
    isok = false;
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
    isok = false;
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
    isok = false;
    }
  else
    {
    $("ImWide").style.backgroundColor=colorWhite;
    args += ",width=" + $("ImWide").value;
    newDimX = $("ImWide").value * 1;
    }

  if ($("ImHigh").value.strip().empty()) $("ImHigh").value = "NotSpecified";
  if ($("ImHigh").value == "NotSpecified")
    {
    $("ImHigh").style.backgroundColor=colorRed;
    isok = false;
    }
  else
    {
    $("ImHigh").style.backgroundColor=colorWhite;
    args += ",height=" + $("ImHigh").value;
    newDimY = $("ImHigh").value * 1;
    }

  // Special check to set the All Segments flag if the HmiB2ptr processing is selected
  if ($("OptionHmiB2ptr").checked)
    args += ",A=1";

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
// XXXXXXX need to set sizeratio based on high and wide, can use keyCDELT if needed

  if (processingFirstRecord)
    {
    var dims = processingFirstRecord.segments[0].dims[0];
    var oldDimX = dims.substring(0,dims.indexOf("x"));
    var oldDimY = dims.substring(dims.indexOf("x")+1);
    localImPatchSizeRatio = (newDimX * newDimY) / (oldDimX * oldDimY);
    }
  ExportProcessingOptions[ImPatchOption].Size = localImPatchSizeRatio;
  this.paramsValid = (isok ? args : "");
  CheckRediness();
  return (isok ? args : "");
  }


// Set display defaults
function ImPatchInit(onLoad)
  {
    if (onLoad)
    {
        // Truly reset ImPatch parameters to default values.
        noaaColor = colorWhite;
        var requireColor = colorRed;

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
        processingFirstRecord = null;
        processingLastRecord = null;
        defaultStartUsed = 1;
        defaultStopUsed = 1;
    }
    else
    {
        // set first and last record for impatch

        // blow away arguments cache since we are changing the values of processingFirstRecord and processingLastRecord
        processingGetRecInfo(this.Check, this, 1); // changes this.argsReady to true
        processingGetRecInfo(this.Check, this, -1);
    }
  }

// End of IM Patch code

//
// Maproj - map projections
//

var MaprojOption;
var newDimX = 0;
var newDimY = 0;

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
            thisLong = NOAA_rslist.keywords[keyNoaaLon].values[irec];
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
            NOAA_TRef =  NOAA_rslist.keywords[keyNoaaTime].values[minRec];
            $("MaprojTRef").innerHTML = "at " + NOAA_TRef;
            $("MaprojX").value = NOAA_rslist.keywords[keyNoaaLonX].values[minRec];
            $("MaprojY").value = NOAA_rslist.keywords[keyNoaaLat].values[minRec];
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
      onComplete: function() { }
      });
    }
  }

function MaprojSet(param)
  {
    var checkRes = null;

    this.argsReady = false;
    this.paramsValid = null;

    // Get first and last record info
    if (param == 1) defaultStartUsed = 0;
    if (param == 2) defaultStopUsed = 0;

    if (processingFirstRecord && processingLastRecord)
    {
        checkRes = this.Check(false);

        if (checkRes === null)
        {
            this.argsReady = true;
            return 1;
        }
        else if (checkRes.length == 0)
        {
            this.argsReady = true;
            return 0;
        }

        this.argsReady = true;
        return 0;
    }

    // This looks like it will call potentially MaprojCheck() twice, which is not ideal.
    if (!processingFirstRecord)
    {
        processingGetRecInfo(this.Check, this, 1);
    }
    if (!processingLastRecord)
    {
        processingGetRecInfo(this.Check, this, -1);
    }

    return 0;
  }

function MaprojCheck(fromSetProcessing)
  {
  var MaprojSizeRatio = 1.0;
  var isok = true;
  var args = "Maproj";
  var MaprojLocOption;

    if (!$("OptionMaproj").checked)
    {
        // If this processing option is not selected, then do not check parameter values.
        return '';
    }

    if (fromSetProcessing)
    {
        if (!this.argsReady)
        {
            // a Set() is pending
            return null;
        }
    }

    if (this.paramsValid !== null)
    {
        // Already checked.
        return this.paramsValid; // Empty string if the previous call determined the arguments to be invalid.
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

  if (processingFirstRecord)
    {
    // var DegPerArcsec = (180.0/(696.0*3.14159))*(149640.0/(180.0*3600.0/3.14159)  // deg/Mm * Mm/arcsec
    var DegPerArcsec = 215/3600.0;  // deg/RsunRadian  * AURadian/arcsec
    var DegPerPixel = new Number(processingFirstRecord.keywords[keyCDELT].values[0] * DegPerArcsec);
    $("MaprojMaxScale").innerHTML = DegPerPixel.toPrecision(3);
    $("MaprojMaxScale").style.backgroundColor=colorWhite;
    if ($("MaprojScale").value == "NotSpecified") $("MaprojScale").value = $("MaprojMaxScale").innerHTML;
    }

  if ($("MaprojLnObs").checked) $("MaprojX").value = processingFirstRecord.keywords[keyCRLN].values[0];
  if ($("MaprojX").value.strip().empty()) $("MaprojX").value = "NotSpecified";
  if ($("MaprojX").value == "NotSpecified")
    {
    isok = false;
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
    isok = false;
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
    isok = false;
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
    isok = false;
    }
  else
    {
    $("MaprojWide").style.backgroundColor=colorWhite;
    args += ",cols=" + $("MaprojWide").value;
    newDimX = $("MaprojWide").value * 1;
    }

  if ($("MaprojHigh").value.strip().empty())
    $("MaprojHigh").value = "NotSpecified";
  if ($("MaprojHigh").value == "NotSpecified")
    {
    $("MaprojHigh").style.backgroundColor=colorRed;
    isok = false;
    }
  else
    {
    $("MaprojHigh").style.backgroundColor=colorWhite;
    args += ",rows=" + $("MaprojHigh").value;
    newDimY = $("MaprojHigh").value * 1;
    }

  // Special check to set the All Segments flag if the HmiB2ptr processing is selected
  if ($("OptionHmiB2ptr").checked)
    args += ",A=1";

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

  if (processingFirstRecord)
    {
    var dims = processingFirstRecord.segments[0].dims[0];
    var oldDimX = dims.substring(0,dims.indexOf("x"));
    var oldDimY = dims.substring(dims.indexOf("x")+1);
    MaprojSizeRatio = (newDimX * newDimY) / (oldDimX * oldDimY);
    }

  ExportProcessingOptions[MaprojOption].Size = MaprojSizeRatio;
  this.paramsValid = (isok ? args : "");
  CheckRediness();
  return (isok ? args : "");
  }

// Set display defaults

function MaprojInit(onLoad)
  {
    if (onLoad)
    {
        noaaColor = colorWhite;
        var requireColor = colorRed;
        $("ProcessMaproj").style.display="none";
        $("MaprojReflatLine").style.display="none";
        $("MaprojNOAA").style.backgroundColor=colorWhite; $("MaprojNOAA").value = "NotSpecified";
        $("MaprojScale").style.backgroundColor = requireColor; $("MaprojScale").value = "NotSpecified";
        $("MaprojX").style.backgroundColor = requireColor; $("MaprojX").value = "NotSpecified";
        $("MaprojLnObs").checked = 0;
        $("MaprojY").style.backgroundColor = requireColor; $("MaprojY").value = "NotSpecified";
        $("MaprojWide").style.backgroundColor = requireColor; $("MaprojWide").value = "NotSpecified";
        $("MaprojHigh").style.backgroundColor = requireColor; $("MaprojHigh").value = "NotSpecified";
        processingFirstRecord = null;
        processingLastRecord = null;
        defaultStartUsed = 1;
        defaultStopUsed = 1;
    }
    else
    {
        processingGetRecInfo(this.Check, this, 1);
        processingGetRecInfo(this.Check, this, -1);
    }
  }

// End of MapProj code

// HTML element event handlers

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
  ExpOpt.Size = 1.0;
  ExpOpt.argsReady = true; // no Set() call pending
  ExpOpt.paramsValid = null; // cache last Check() call
  ExportProcessingOptions[iOpt] = ExpOpt;

  iOpt++;
  ProcessingOptionsHTML +=
    '<input type="checkbox" checked="false" value="aia_scale" id="OptionAiaScale" onChange="SetProcessing('+iOpt+');" /> ' +
    '<label id="OptionAiaScaleLabel" for="OptionAiaScale">aia_scale - Scale image to 0.6 arcsec/pixel</label><br>';
  ExpOpt = new Object();
  ExpOpt.id = "OptionAiaScale";
  ExpOpt.rowid = "ProcessAiaScale";
  ExpOpt.Init = AiaScaleInit;
  ExpOpt.Check = AiaScaleCheck;
  ExpOpt.Set = AiaScaleSet;
  ExpOpt.Size = 1.0;
    ExpOpt.argsReady = true;
    ExpOpt.paramsValid = null;
    // to determine if checkbox is disabled, must check record-set keyword values
    ExpOpt.disabled_state_rec_dep = true;
  ExportProcessingOptions[iOpt] = ExpOpt;
  AiaScaleOption = iOpt;

  iOpt++;
  ProcessingOptionsHTML +=
    '<input type="checkbox" checked="false" value="HmiB2ptr" id="OptionHmiB2ptr" onChange="SetProcessing('+iOpt+');" /> ' +
    'HmiB2ptr - Convert HMI B_720s to Phi,Theta,R coordinates<br>'
  ExpOpt = new Object();
  ExpOpt.id = "OptionHmiB2ptr";
  ExpOpt.rowid = "ProcessHmiB2ptr";
  ExpOpt.Init = HmiB2ptrInit;
  ExpOpt.Check = HmiB2ptrCheck;
  ExpOpt.Set = HmiB2ptrSet;
  ExpOpt.Size = 1.0;
    ExpOpt.argsReady = true;
    ExpOpt.paramsValid = null;
  ExportProcessingOptions[iOpt] = ExpOpt;
  HmiB2ptrOption = iOpt;

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
  ExpOpt.Size = 1.0;
    ExpOpt.argsReady = true;
    ExpOpt.paramsValid = null;
  ExportProcessingOptions[iOpt] = ExpOpt;
  ResizeOption = iOpt;

  iOpt++;
  ProcessingOptionsHTML +=
    '<input type="checkbox" checked="false" value="im_patch" id="OptionImPatch" onChange="SetProcessing('+iOpt+');" /> ' +
    'im_patch - Extract sub-frame<br>';
  ExpOpt = new Object();
  ExpOpt.id = "OptionImPatch";
  ExpOpt.rowid = "ProcessImPatch";
  ExpOpt.Init = ImPatchInit;
  ExpOpt.Check = ImPatchCheck;
  ExpOpt.Set = ImPatchSet;
  ExpOpt.Size = 1.0;
    ExpOpt.argsReady = true;
    ExpOpt.paramsValid = null;
  ExportProcessingOptions[iOpt] = ExpOpt;
  ImPatchOption = iOpt;

    // add processing HTML element event handlers
    $('ResetImPatchButton').addEventListener('click', function() { ImPatchInit(1); ExportProcessingOptions[ImPatchOption].Set(0); $("ProcessImPatch").style.display="table-row";}, false);


  iOpt++;
  ProcessingOptionsHTML +=
    '<input type="checkbox" checked="false" value="maproj" id="OptionMaproj" onChange="SetProcessing('+iOpt+');" /> ' +
    'maproj - Extract a sub-frame and remap to a chosen projection.<br>';
  ExpOpt = new Object();
  ExpOpt.id = "OptionMaproj";
  ExpOpt.rowid = "ProcessMaproj";
  ExpOpt.Init = MaprojInit;
  ExpOpt.Check = MaprojCheck;
  ExpOpt.Set = MaprojSet;
  ExpOpt.Size = 1.0;
    ExpOpt.argsReady = true;
    ExpOpt.paramsValid = null;
  ExportProcessingOptions[iOpt] = ExpOpt;
  MaprojOption = iOpt;

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
  ExpOpt.Size = 1.0;
    ExpOpt.argsReady = true;
    ExpOpt.paramsValid = null;
  ExportProcessingOptions[iOpt] = ExpOpt;
  RebinOption = iOpt;

  $("ExportProcessing").innerHTML = ProcessingOptionsHTML;

  var nOpt = iOpt + 1;

    for (iOpt=1; iOpt<nOpt; iOpt++)
    {
        var id;

        ExpOpt = ExportProcessingOptions[iOpt];
        ExpOpt.Init(1);
        $(ExpOpt.id).checked = false;
        $(ExpOpt.rowid).style.display = "none";

        if (ExpOpt.hasOwnProperty('disabled_state_rec_dep'))
        {
            $(ExpOpt.id).store({ 'disabled_state_rec_dep' : ExpOpt.disabled_state_rec_dep });
        }
        else
        {
            $(ExpOpt.id).store({ 'disabled_state_rec_dep' : false });
        }
    }
}

function ProcessingOK()
{
    var processingIsOk = true;
    var iOpt;

    // Skip the no-op processing step (index 0).
    for (iOpt = 1; iOpt < ExportProcessingOptions.length; iOpt++)
    {
        ExpOpt = ExportProcessingOptions[iOpt];

        if ($(ExpOpt.id).checked && (ExpOpt.paramsValid === null || ExpOpt.paramsValid.length == 0))
        {
            processingIsOk = false;
            break;
        }
    }

    return processingIsOk;
}

function ExportProcessingArgsReady()
{
    var argsReady = true;
    var iOpt;

    // Skip the no-op processing step (index 0).
    for (iOpt = 1; iOpt < ExportProcessingOptions.length; iOpt++)
    {
        ExpOpt = ExportProcessingOptions[iOpt];

        if ($(ExpOpt.id).checked && !ExpOpt.argsReady)
        {
            argsReady = false;
            break;
        }
    }

    return argsReady;
}
