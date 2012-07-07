// Functions for export protocol options
//
// CheckXXXX creates the contents of ExportProtocolArgs from parameters for
// this export type
//
// Protocol details for Images
//

var ColorInfo;

function ProtocolsGetImageInfo()
  {
  $("AjaxBusy").innerHTML = Ajax.activeRequestCount;
  new Ajax.Request('http://' + Host + '/cgi-bin/ajax/jsoc_info_jsoc2',
      {
      method: 'get',
      parameters: {"op" : "rs_list", "ds": "su_priya.ColorTables[][]", "key": "InSeries,WAVELNTH,CT_name,min,max,scaling" },
      onSuccess: function(transport, json)
        {
        var response = transport.responseText || "no response text";
        ColorInfo = response.evalJSON();
        if (ColorInfo.status > 0 || ColorInfo.count == 0) 
          {
          insertOption("ImageCT", "No Colortables Found");
          }
        else
          {
          for (var i=0; i<ColorInfo.count; i++)
            {
            insertOption("ImageCT",ColorInfo.keywords[2].values[i],"");
            if (ColorInfo.keywords[0].values[i] == SeriesName)
              {
              $("ImageMin").value = ColorInfo.keywords[3].values[i];
              $("ImageMax").value = ColorInfo.keywords[4].values[i];
              $("ImageSer").innerHTML = ColorInfo.keywords[0].values[i];
              $("ImageCT").selectedIndex = i;
              }
            }
          }
        },
      onFailure: function()
        {
        alert('Something went wrong with color table data request');
        $("ImageCT").value = "Not Found";
        },
      onComplete: function() { $("AjaxBusy").innerHTML = Ajax.activeRequestCount; }
      });
  }

function ProtocolOptionsInit()
  {
  // $("ProtocolOptions").style.display="none";
  $("ImagingOptions").style.display="none";
  }

function ProtocolImageInit(imageType)
  {
  ProtocolsGetImageInfo();
  ProtocolOptionsSet = 0;
  // $("ProtocolOptions").style.display="table-row";
  $("ImagingOptions").style.display="table-row";
  }


function ProtocolSetMinMax()
  {
  var protocol_args = "min="+$("ImageMin").value + ",max="+$("ImageMax").value;
  var protocol_base = $("ExportProtocolHidden").value.substr(0,3);
  $("ExportProtocolHidden").value = protocol_base + "," + protocol_args;
  ProtocolOptionsSet = 1;
  }

function ProtocolSetColorTable()
  {
  var i = $("ImageCT").selectedIndex;
  $("ImageMin").value = ColorInfo.keywords[3].values[i];
  $("ImageMax").value = ColorInfo.keywords[4].values[i];
  $("ImageSer").innerHTML = ColorInfo.keywords[0].values[i];
  ProtocolSetMinMax();
  }

// End of Protocol options code
