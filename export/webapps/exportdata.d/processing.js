// Functions for export options
//
// CheckXXXX creates the contents of ExportProcessingArgs from parameters for
// this export type

//
// Processing details for Patch in Heliographic Coords
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
          $("HgLocType").value = "stony";
          $("HgTRef").value = minTime;
          $("HgX").value = minLong + "";
          $("HgY").value = minLat + "";
	  $("HgNOAA").style.backgroundColor="#FFCC66";
          noaaColor = "#D8D8D8";
          $("HgLocType").style.backgroundColor = noaaColor;
          $("HgTRef").style.backgroundColor = noaaColor;
          $("HgX").style.backgroundColor = noaaColor;
          $("HgY").style.backgroundColor = noaaColor;
	  CheckHgPatch();
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
  
function CheckHgPatch()
  {
  var isok = 0;
  CheckHgRecordSet(0);
  ExportProcessingArgs = "hg_patch,";
  if ($("HgTrack").checked) // checked for tracking, default
    {
    ExportProcessingArgs = ExportProcessingArgs + "t=0,";
    HgTracked = 1;
    }
  else
    {
    ExportProcessingArgs = ExportProcessingArgs + "t=1,";
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
    ExportProcessingArgs = ExportProcessingArgs + "t_start=" + $("HgTStart").value + ",";
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
    ExportProcessingArgs = ExportProcessingArgs + "t_stop=" + $("HgTStop").value + ",";
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
    ExportProcessingArgs = ExportProcessingArgs + "cadence=" + $("HgTDelta").value + ",";
    isok += 1;
    }

  $("HgLocType").style.backgroundColor=noaaColor;
  ExportProcessingArgs = ExportProcessingArgs + "locunits=" + $("HgLocType").value + ",";
  $("HgBoxType").style.backgroundColor="#FFFFFF";
  ExportProcessingArgs = ExportProcessingArgs + "boxunits=" + $("HgBoxType").value + ",";

  if ($("HgLocType").value == "carrlong" && HgTracked == 1)
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
    ExportProcessingArgs = ExportProcessingArgs + "car_rot=" + $("HgCarrot").value + ",";
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
    ExportProcessingArgs = ExportProcessingArgs + "t_ref=" + $("HgTRef").value + ",";
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
    ExportProcessingArgs = ExportProcessingArgs + "x=" + $("HgX").value + ",";
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
    ExportProcessingArgs = ExportProcessingArgs + "y=" + $("HgY").value + ",";
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
    ExportProcessingArgs = ExportProcessingArgs + "width=" + $("HgWide").value + ",";
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
    ExportProcessingArgs = ExportProcessingArgs + "height=" + $("HgHigh").value + ",";
    isok += 1;
    }
  return (isok > 0);
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
    alert("Please select valid hg_patch available series from list.");
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

function HgPatchInit()
  {
  noaaColor = "#FFFFFF";
  var requireColor = "#D88080";
  $("ProcessHgPatch").style.display="none";
  $("HgTrack").checked = 1;
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
  HgPatchActive = 0;
  HgGetSeriesList();
  HgSeriesSelected = 0;
  HgTracked = 1;
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
// Processing details for normalized Scaling
//

function AiaScaleInit()
  {
  $("ProcessAiaScale").style.display="none";
  }

//
// Processing details for all options
//

function ProcessingInit()
{
    AiaScaleInit();
    HgPatchInit();
    $("ProcessRebin").style.display="none";
}

function RebinInit()
{
    // Clear out any old settings.
    $("RebinMethod").value = "boxcar";
    $("RebinCrop").checked = false;
    $("RebinRotate").checked = false;
    $("RebinScale").value = "1.0";
    
    $("RebinFWHMRow").style.display = "none";
    $("RebinNvectorRow").style.display = "none";
    
    // Display the rebin options table.
    $("ProcessRebin").style.display="table-row";   
}

function CheckRebin(control)
{
    if (control == "method")
    {
        if ($("RebinMethod").value == "gaussian")
        {
            $("RebinFWHMRow").style.display = "table-row";
            $("RebinNvectorRow").style.display = "table-row";
        }
        else if ($("RebinMethod").value == "boxcar")
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
}

function GetRebinArgs()
{
    var rv = "rebin";
    
    if ($("RebinCrop").checked)
    {
        rv = rv + ",";
        rv = rv + "-c";
    }
    
    if ($("RebinRotate").checked)
    {
        rv = rv + ",";
        rv = rv + "-u";
    }
    
    if (parseFloat($("RebinScale").value) != 1.0)
    {
        rv = rv + ",";
        rv = rv + "scale=" + $("RebinScale").value;
    }
    
    if ($("RebinMethod").value == "gaussian")
    {
        if (parseFloat($("RebinFWHM").value) != -1.0)
        {
            rv = rv + ",";
            rv = rv + "FHWM=" + $("RebinFWHM").value;
        }
        
        if (parseFloat($("RebinNvector").value) != -1.0)
        {
            rv = rv + ",";
            rv = rv + "nvector=" +  $("RebinNvector").value;
        }
    }
    
    if ($("RebinMethod").value != "boxcar")
    {
        rv = rv + ",";
        rv = rv + "method=" + $("RebinMethod").value;
    }
    
    return rv;
}