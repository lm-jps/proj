// Functions for export protocol options
//
// CheckXXXX creates the contents of ExportProtocolArgs from parameters for
// this export type
//
// Protocol details for Images

// Global vars
var ColorInfo;
var IsCT = 0;
var Protocol_ctnames = {};
var Protocol_scalenames = {};
var WaveLength = 0;

function ProtocolGetImageDefaults()
  {
  // For now this gets ONLY the value of WAVELNTH for the first record in the export recordset
  $("AjaxBusy").innerHTML = Ajax.activeRequestCount;
  // $("RSCountPlace").innerHTML = "Getting count - wait...";
  WaveLength = 0;
  var RecordSet = $("ExportRecordSet").value;
  var SampleRecord;
  new Ajax.Request('http://' + Host + '/cgi-bin/ajax/' + JSOC_INFO,
    {
    method: 'get',
    parameters: {"ds" : RecordSet, "op" : "rs_list", "key" : "WAVELNTH", "n" : "1" },

    onSuccess: function(transport, json)
      {
      var response = transport.responseText || "no response text";
      SampleRecord = response.evalJSON();
      if (SampleRecord.status == 0)
        {
        var nkeys = SampleRecord.keywords.length;
        for (var ikey=0; ikey<nkeys; ikey++)
          {
          if (SampleRecord.keywords[ikey].name == "WAVELNTH")
            {
            WaveLength = Math.round(SampleRecord.keywords[ikey].values[0]);
            ProtocolImageDefaults();
            break;
            }
          }
        }
      else
        ProtocolImageDefaults();
      $("AjaxBusy").innerHTML = Ajax.activeRequestCount;
      },
    onFailure: function() { alert('Something went wrong...'); },
    onComplete: function() { $("AjaxBusy").innerHTML = Ajax.activeRequestCount; }
    });
  }
 
function ProtocolImageDefaults()
  {
  var foundCT = 0;
  if (IsCT == 1)
    {
    for (var i=0; i < ColorInfo.count; i++)
      {
      if (ColorInfo.keywords[0].values[i].toLowerCase() == SeriesName.toLowerCase() &&
          ColorInfo.keywords[1].values[i] == WaveLength)
        {          
        $("ImageSer").innerHTML = ColorInfo.keywords[0].values[i];               
        $("ImageCT").selectedIndex = Protocol_ctnames[ColorInfo.keywords[2].values[i]];
        $("ImageMin").value = ColorInfo.keywords[3].values[i];                 
        $("ImageMax").value = ColorInfo.keywords[4].values[i];
        $("ImageScl").selectedIndex = Protocol_scalenames[ColorInfo.keywords[5].values[i]];
        foundCT = 1;
        break;
        }            
      }   //end of for loop     
    }
  if (foundCT == 0)
    {
    $("defaultvals").innerHTML="No colortable specified for this series, use defaults or change:  Min/Max will be calculated from data"
    $("ImageMin").value = "NOT_SPECIFIED"; 
    $("ImageMax").value = "NOT_SPECIFIED";
    $("ImageSer").innerHTML = SeriesName;
    $("ImageScl").selectedIndex= Protocol_scalenames["MINMAX"];
    $("ImageCT").selectedIndex= $("ImageCT").length - 1;  
    }                         
  }

// Get color table info for the series at hand
function ProtocolGetImageInfo()
  { 
  $("AjaxBusy").innerHTML = Ajax.activeRequestCount;
  if (IsCT == 0)
    { // Not initialized yet, find colortables list then set defaults
    new Ajax.Request('http://' + Host + '/cgi-bin/ajax/jsoc_info_jsoc2',
      {   
      method: 'get',
      parameters: {"op" : "rs_list", "ds": "jsoc.Color_Tables[][]", "key": "Inseries,WAVELNTH,CT_name,scalemin,scalemax,scaling" },
      onSuccess: function(transport, json)
        {
        var response = transport.responseText || "no response text";
        ColorInfo = response.evalJSON(); 
        if($("ImageCT").length > 0)
          {
          $("ImageCT").length=0;
          $("ImageScl").length=0;
          } 
        if (ColorInfo.status > 0 || ColorInfo.count == 0) 
          {
          insertOption("ImageCT", "No Colortables Found");
          IsCT = -1;
          }
        else
          {                   
          Protocol_ctnames = {};
          Protocol_scalenames = {};
          var ctindex=0, scaleindex=0;
          // Add new CT name and scale name to drop down options, and set default for each ColorInfo record.
          for (var i=0; i<ColorInfo.count; i++)
            {              
            var ctname=ColorInfo.keywords[2].values[i];
            var scalename=ColorInfo.keywords[5].values[i];
            if (Protocol_ctnames[ctname] == undefined)
              {
              insertOption("ImageCT",ctname, "");                                 
              Protocol_ctnames[ctname] = ctindex++;
              }
            if (Protocol_scalenames[scalename] == undefined)
              {
              insertOption("ImageScl",scalename, "");  
              Protocol_scalenames[scalename] = scaleindex++;
              }
            }           
          // Make sure final CT is good default
          insertOption("ImageCT","grey.sao", "");                                 
          IsCT = 1;
          }
        ProtocolGetImageDefaults();
        },
      onFailure: function()
        {
        alert('Something went wrong with color table data request');
        $("ImageCT").value = "Not Found";
        },
      onComplete: function() { $("AjaxBusy").innerHTML = Ajax.activeRequestCount;}    
      });
    }  
  else // set the defaults
    ProtocolGetImageDefaults();
  };

function ProtocolOptionsInit()
  {
  $("ProtocolImageOptions").style.display="none";
  ProtocolOptionsSet = 0;
  }


function ProtocolImageUserSet(opt)
  {
  if (opt == 1)
    {
    var optindex = $("ImageCT").selectedIndex;
    var optname = $("ImageCT").options[optindex].value;
    for (var i=0; i < ColorInfo.count; i++)
      {
      if (ColorInfo.keywords[2].values[i] == optname)
        {          
        $("ImageCTTxt").innerHTML = ColorInfo.keywords[0].values[i];
        $("ImageMinTxt").innerHTML = ColorInfo.keywords[3].values[i];                 
        $("ImageMaxTxt").innerHTML = ColorInfo.keywords[4].values[i];
        $("ImageSclTxt").innerHTML = ColorInfo.keywords[5].values[i];
        break;
        }
      }
    }
  ProtocolOptionsSet = 2;
  }

function ProtocolImageInit(imageType)
  {
  if (ProtocolOptionsSet == 0)
    {
    ProtocolGetImageInfo();
    ProtocolOptionsSet = 1;
    }
  $("ProtocolImageOptions").style.display="table-row";   
  }

function ProtocolImageCheck()
  {
  var protocol_base = $("ExportProtocolHidden").value.substr(0,3);
  $("ExportProtocolHidden").value = protocol_base;
  if (ProtocolOptionsSet > 0) // Always provide full set pf parameters.
    {
    var min= parseFloat($("ImageMin").value);
    var max= parseFloat($("ImageMax").value);
    var protocol_args = ",min="+$("ImageMin").value + ",max="+$("ImageMax").value +
       ",CT="+ $("ImageCT").value + ",scaling=" +$("ImageScl").value + 
       ",size=" +$("ImageSize").value.substr(4);
    $("ExportProtocolHidden").value += protocol_args;
    }
  ExportProtocolArgsOK=1;       
  }

// End of Protocol options code
