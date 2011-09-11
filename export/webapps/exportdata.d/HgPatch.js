// Functions for export options
//
// CheckXXXX creates the contents of ExportProcessingArgs from parameters for
// this export type
//
// Processing details for Patch in Heliographic Coords
//

//function HgPatchGetNoaa()
  //{
  // if $("HgNOAA").value is a number
  //   { // lookup NOAA number
  //   do lookdata rs list fetch of su_rsb.NOAA_Regions key=....
  //   find instance nearest CM, get long and lat and time, use stonyhurst since do not know rot
  //   populate $("HgLocType").value, $("HgX").value, $("HgY").value, $("HgTRef").value
  //   call CheckHgPatch()
  //   return
  //   }
  //}

var noaaColor;

function HgPatchGetNoaa()
  {
  var noaaNum = 1 * $("HgNOAA").value;
  if (noaaNum < 7000) 
    {
    noaaNum = noaaNum + 10000; // OK for times after 1996 Jan.
    $("HgNOAA").value = noaaNum + "";
    }
  $("AjaxBusy").innerHTML = Ajax.activeRequestCount;
  new Ajax.Request('http://jsoc2.stanford.edu/cgi-bin/ajax/' + JSOC_INFO,
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

function CheckHgPatch()
  {
  var isok = 0;
  CheckHgRecordSet();
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

  if ($("HgCarrot").value == "NotSpecified" || $("HgCarrot").value.length == 0)
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
  if ($("HgTRef").value == "NotSpecified" || $("HgTRef").value.length == 0)
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
  if (CheckHgRecordSet())
    ExportNewRS();
  }

function CheckHgRecordSet() // If HG Patch selected, convert empty record selector to last record.
  {
  var posbracket = $("ExportRecordSet").value.indexOf("[");
  var currentSeries = $("ExportRecordSet").value.substring(0,posbracket);
  var currentSpec = (posbracket < 0 ? "[$]" : $("ExportRecordSet").value.substring(posbracket));
  var selectedSeries;
  if (HgSeriesSelected)
    selectedSeries = $("HgSerList").selectedIndex;
  else
    {
    var n = $("HgSerList").length;
    var i;
    for (i=0; i<n; i++)
      if (currentSeries == $("HgSerList").options[i].value) break;
    if (i == n)
      {
      alert("Please select valid hg_patch available series from list.");
      return;
      }
    HgSeriesSelected = 1;
    $("HgSerList").selectedIndex = i;
    $("ExportRecordSet").value = currentSeries + currentSpec;
    selectedSeries = i;
    }
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
      var n = HgSeriesList.n;
      if (n < 1) alert("WARNING: No _hgpatch series found.\n"+response);
      for (var i=0; i<n; i++)
         insertOption("HgSerList",HgSeriesList.list[i], "");
      },
    onFailure: function() { alert('Failed to get HgPatch Series List'); },
    onComplete: function() { $("AjaxBusy").innerHTML = Ajax.activeRequestCount; }
    });
  return;
  }


// End of HG Patch code
