// Functions for export options
//
// CheckXXXX creates the contents of ExportProcessingArgs from parameters for
// this export type
//
// Each processing type needs a common name, e.g. HgPatch, used for the first part of
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
    $("RebinFWHM").style.backgroundColor = "#D88080";
    $("RebinNvector").style.backgroundColor = "#D88080";
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
            $("RebinFWHM").style.backgroundColor = "#D88080";
        }
        else
        {
            $("RebinFWHM").style.backgroundColor = "#FFFFFF";
        }
    }
    else if (control == "nvector")
    {
        if (parseFloat($("RebinNvector").value) == -1.0)
        {
            $("RebinNvector").style.backgroundColor = "#D88080";
        }
        else
        {
            $("RebinNvector").style.backgroundColor = "#FFFFFF";
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
    $("ResizeCdelt").style.backgroundColor = "#D88080";
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
    $("ResizeCdelt").style.backgroundColor = "#FFFFFF";
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
  else
    rv += ",center_to=1";
    
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
    
  return (isok ? rv : "");
  }

// End Resize

//
// Process HgPatch
//

var noaaColor;

function HgPatchGetNoaa()
  {
  if ($("HgNOAA").value.strip().empty()) $("HgNOAA").value = "NotSpecified";
  if ($("HgNOAA").value != "NotSpecified")
    {
    var noaaNum = 1 * $("HgNOAA").value;
    if (noaaNum < 7000) 
      {
      noaaNum = noaaNum + 10000; // OK for times after 1996 Jan.
      $("HgNOAA").value = noaaNum + "";
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
        catch(err) { $("HgNOAA").value = noaaNum + " " + err; return; }
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
          $("HgLocType").selectedIndex = 0;
          $("HgTRef").value = minTime;
          $("HgX").value = minLong + "";
          $("HgY").value = minLat + "";
	  $("HgNOAA").style.backgroundColor="#FFCC66";
          noaaColor = "#D8D8D8";
          $("HgLocType").style.backgroundColor = noaaColor;
          $("HgTRef").style.backgroundColor = noaaColor;
          $("HgX").style.backgroundColor = noaaColor;
          $("HgY").style.backgroundColor = noaaColor;
	  // CheckHgPatch();
          }
        catch(err) { $("HgNOAA").value = noaaNum + " " + err; return; }
        },
      onFailure: function()
        {
        alert('Something went wrong with NOAA num data request');
        $("HgNOAA").value = "Not Found";
        },
      onComplete: function() { $("AjaxBusy").innerHTML = Ajax.activeRequestCount; }
      });
    }
  }
  
function HgPatchCheck()
  {
  var isok = 0;
  var args = "hg_patch";
  var HgLocOption;
  CheckHgRecordSet(0);
  if ($("HgTrack").checked) // checked for tracking, default
    {
    args += ",t=0";
    HgTracked = 1;
    }
  else
    {
    args += ",t=1";
    $("HgTRef").value = $("HgTStart").value;
    HgTracked = 0;
    }
  if ($("HgTStart").value.strip().empty()) $("HgTStart").value = "NotSpecified";
  if ($("HgTStart").value == "NotSpecified")
    {
    $("HgTStart").style.backgroundColor="#FFFFFF";
    isok = 1;
    }
  else
    {
    $("HgTStart").style.backgroundColor="#FFFFFF";
    args += ",t_start=" + $("HgTStart").value;
    isok += 1;
    }
  if ($("HgTStop").value.strip().empty()) $("HgTStop").value = "NotSpecified";
  if ($("HgTStop").value == "NotSpecified")
    {
    $("HgTStop").style.backgroundColor="#FFFFFF";
    isok = 1;
    }
  else
    {
    $("HgTStop").style.backgroundColor="#FFFFFF";
    args += ",t_stop=" + $("HgTStop").value;
    isok += 1;
    }
  if ($("HgTDelta").value.strip().empty()) $("HgTDelta").value = "NotSpecified";
  if ($("HgTDelta").value == "NotSpecified")
    {
    $("HgTDelta").style.backgroundColor="#FFFFFF";
    isok = 1;
    }
  else
    {
    $("HgTDelta").style.backgroundColor="#FFFFFF";
    args += ",cadence=" + $("HgTDelta").value;
    isok += 1;
    }

  $("HgLocType").style.backgroundColor=noaaColor;
  HgLocOption = $("HgLocType").selectedIndex;
  args += ",locunits=" + $("HgLocType").options[HgLocOption].value;
  $("HgBoxType").style.backgroundColor="#FFFFFF";
  args += ",boxunits=" + $("HgBoxType").options[$("HgBoxType").selectedIndex].value;

  if ($("HgLocType").options[HgLocOption].value == "carrlong" && HgTracked == 1)
    {
    $("HgCarrLi").style.display = "table-row";
    $("HgTRefLi").style.display = "none";
    $("HgTRef").value == "NotSpecified";
    }
  else
    {
    $("HgCarrLi").style.display = "none";
    $("HgTRefLi").style.display = "table-row";
    $("HgCarrot").value == "NotSpecified";
    }

  if ($("HgCarrot").value.strip().empty()) $("HgCarrot").value = "NotSpecified";
  if ($("HgCarrot").value == "NotSpecified" )
    {
    $("HgCarrot").style.backgroundColor="#D88080";
    isok = 0;
    }
  else
    {
    $("HgCarrot").style.backgroundColor="#FFFFFF";
    args += ",car_rot=" + $("HgCarrot").value;
    isok += 1;
    }
  if ($("HgTRef").value.strip().empty()) $("HgTRef").value = "NotSpecified";
  if ($("HgTRef").value == "NotSpecified")
    {
    $("HgTRef").style.backgroundColor="#D88080";
    isok = 0;
    }
  else
    {
    $("HgTRef").style.backgroundColor=noaaColor;
    args += ",t_ref=" + $("HgTRef").value;
    isok += 1;
    }
  if ($("HgX").value.strip().empty()) $("HgX").value = "NotSpecified";
  if ($("HgX").value == "NotSpecified")
    {
    $("HgX").style.backgroundColor="#D88080";
    isok = 0;
    }
  else
    {
    $("HgX").style.backgroundColor=noaaColor;
    args += ",x=" + $("HgX").value;
    isok += 1;
    }
  if ($("HgY").value.strip().empty()) $("HgY").value = "NotSpecified";
  if ($("HgY").value == "NotSpecified")
    {
    $("HgY").style.backgroundColor="#D88080";
    isok = 0;
    }
  else
    {
    $("HgY").style.backgroundColor=noaaColor;
    args += ",y=" + $("HgY").value;
    isok += 1;
    }
  if ($("HgWide").value.strip().empty()) $("HgWide").value = "NotSpecified";
  if ($("HgWide").value == "NotSpecified")
    {
    $("HgWide").style.backgroundColor="#D88080";
    isok = 0;
    }
  else
    {
    $("HgWide").style.backgroundColor="#FFFFFF";
    args += ",width=" + $("HgWide").value;
    isok += 1;
    }
  if ($("HgHigh").value.strip().empty()) $("HgHigh").value = "NotSpecified";
  if ($("HgHigh").value == "NotSpecified")
    {
    $("HgHigh").style.backgroundColor="#D88080";
    isok = 0;
    }
  else
    {
    $("HgHigh").style.backgroundColor="#FFFFFF";
    args += ",height=" + $("HgHigh").value;
    isok += 1;
    }

  return (isok ? args : "");
  }

// Special action to modify input RecordSet if no records spec present
// Also gets list of available series.

function HgSeriesSelect()
  {
  HgSeriesSelected = 1;
  CheckHgRecordSet(1);
  //if (CheckHgRecordSet())
  //ExportNewRS();
  }

function CheckHgRecordSet(clicked) // If HG Patch selected, convert empty record selector to last record.
  {
  var posbracket = $("ExportRecordSet").value.indexOf("[");
  var currentSeries = $("ExportRecordSet").value.substring(0,posbracket);
  var currentSpec = (posbracket < 0 ? "[$]" : $("ExportRecordSet").value.substring(posbracket));
  var selectedSeries;
  var n = $("HgSerList").length;
  if (HgSeriesSelected)
    selectedSeries = $("HgSerList").selectedIndex;
  else
    {
    for (selectedSeries=1; selectedSeries<n; selectedSeries++)
      if (currentSeries == $("HgSerList").options[selectedSeries].value) break;
    }
  if ((clicked && selectedSeries == 0) || selectedSeries == n)
    {
    HgSeriesSelected = 0;
    $("HgSerList").options[0].style.backgroundColor = "Red";
    alert("Please select valid XXX hg_patch available series from list.");
    return(0);
    }
  HgSeriesSelected = 1;
  $("HgSerList").selectedIndex = selectedSeries;
  currentSeries = $("HgSerList").options[selectedSeries].value;
  $("ExportRecordSet").value = currentSeries + currentSpec;
  if (posbracket < 0)
    {
    ExportNewRS();
    return(0);
    }
  return(1);
  }

// Set display defaults

function HgPatchInit(isActive)
  {
  if (!isActive)
    {
    noaaColor = "#FFFFFF";
    var requireColor = "#D88080";
    $("ProcessHgPatch").style.display="none";
    $("HgTrack").checked = true;
    $("HgNOAA").style.backgroundColor="#FFFFFF"; $("HgNOAA").value = "NotSpecified";
    $("HgTStart").value = "NotSpecified";
    $("HgTStop").value = "NotSpecified";
    $("HgTDelta").value = "NotSpecified";
    $("HgTRef").style.backgroundColor = requireColor; $("HgTRef").value = "NotSpecified";
    $("HgCarrot").style.backgroundColor = requireColor; $("HgCarrot").value = "NotSpecified";
    $("HgX").style.backgroundColor = requireColor; $("HgX").value = "NotSpecified";
    $("HgY").style.backgroundColor = requireColor; $("HgY").value = "NotSpecified";
    $("HgWide").style.backgroundColor = requireColor; $("HgWide").value = "NotSpecified";
    $("HgHigh").style.backgroundColor = requireColor; $("HgHigh").value = "NotSpecified";
    HgGetSeriesList();
    HgTracked = 1;
    }
  HgSeriesSelected = 0;
  }

function HgGetSeriesList()
  {
  var HgPatchSerList = new Ajax.Request('http://' + Host + '/cgi-bin/ajax/show_hgpatch',
    {
    method: 'get',
    onSuccess: function(transport, json)
      {
      var response = transport.responseText || "no HgPatch series";
      HgSeriesList = response.evalJSON();
      for (var i=$("HgSerList").length; i > 0; i--)
        $("HgSerList").remove(i-1);
      var n = HgSeriesList.n;
      if (n < 1) alert("WARNING: No _hgpatch series found.\n"+response);
      insertOption("HgSerList","Not Selected Yet", "");
      for (var i=0; i<n; i++)
         insertOption("HgSerList",HgSeriesList.list[i], "");
      },
    onFailure: function() { alert('Failed to get HgPatch Series List'); },
    onComplete: function() { $("AjaxBusy").innerHTML = Ajax.activeRequestCount; }
    });
  return;
  }

// End of HG Patch code


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
    'no_op - none<br>';
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
    '<input type="checkbox" checked="false" value="hg_patch" id="OptionHgPatch" onChange="SetProcessing('+iOpt+');" /> ' +
    'hg_patch - Extract sub-frame<br>';
  ExpOpt = new Object();
  ExpOpt.id = "OptionHgPatch";
  ExpOpt.rowid = "ProcessHgPatch";
  ExpOpt.Init = HgPatchInit;
  ExpOpt.Check = HgPatchCheck;
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

