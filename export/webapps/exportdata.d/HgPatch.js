// Functions for export options
//
// CheckXXXX creates the contents of ExportProcessingArgs from parameters for
// this export type
//
// Processing details for Patch in Heliographic Coords
//

function CheckHgPatch()
  {
  var isok = 0;
  ExportProcessingArgs = "hg_patch,";
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

  $("HgLocType").style.backgroundColor="#FFFFFF";
  ExportProcessingArgs = ExportProcessingArgs + "locunits=" + $("HgLocType").value + ",";
  $("HgBoxType").style.backgroundColor="#FFFFFF";
  ExportProcessingArgs = ExportProcessingArgs + "boxunits=" + $("HgBoxType").value;

  if ($("HgLocType").value == "carrlong")
    {
    $("HgCarrLi").style.display = "inline";
    $("HgTRefLi").style.display = "none";
    $("HgTRef").value == "NotSpecified";
    }
  else
    {
    $("HgCarrLi").style.display = "none";
    $("HgTRefLi").style.display = "inline";
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
    $("HgTRef").style.backgroundColor="#FFFFFF";
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
    $("HgX").style.backgroundColor="#FFFFFF";
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
    $("HgY").style.backgroundColor="#FFFFFF";
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

function CheckHgRecordset() // If HG Patch selected, convert empty record selector to last record.
  {
  var posbracket = $("ExportRecordSet").value.indexOf("[");
  if (posbracket == -1)
    {
    $("ExportRecordSet").value = $("ExportRecordSet").value + "[$]";
    ExportNewRS();
    }
  }

// Set display defaults

function HgPatchInit()
  {
  $("ProcessHgPatch").style.display="none";
  $("HgTStart").value = "NotSpecified";
  $("HgTStop").value = "NotSpecified";
  $("HgTDelta").value = "NotSpecified";
  $("HgTRef").style.backgroundColor = "#D88080"; $("HgTRef").value = "NotSpecified";
  $("HgCarrot").style.backgroundColor = "#D88080"; $("HgCarrot").value = "NotSpecified";
  $("HgX").style.backgroundColor = "#D88080"; $("HgX").value = "NotSpecified";
  $("HgY").style.backgroundColor = "#D88080"; $("HgY").value = "NotSpecified";
  $("HgWide").style.backgroundColor = "#D88080"; $("HgWide").value = "NotSpecified";
  $("HgHigh").style.backgroundColor = "#D88080"; $("HgHigh").value = "NotSpecified";
  HgPatchActive = 0;
  }

// End of HG Patch code
