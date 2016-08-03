#!/home/jsoc/anaconda3/bin/python3

import cgi
import math
from datetime import datetime
import sys
import os
import mpld3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.time import Time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../../base/libs/py'))
from drmsCmdl import CmdlParser
import drms_json as drms
from matplotlib.widgets import Button
import json
import matplotlib.colors as mcol
from mpl_toolkits.axes_grid1 import AxesGrid

if sys.version_info < (3, 0):
    raise Exception("You must run the 3.0 release, or a more recent release, of Python.")
    
# JSON return status.

JSON_RV_SUCCESS = 0
JSON_RV_BADARGS = 1
JSON_RV_CREATEIMGFAIL = 2

#radii is hard coded for the test. need to have the radii inside the data
slice = 0
radii = [1.001,1.050,1.100,1.150,1.250,1.500,1.750,2.000,2.250,2.499]
manifest = {}
carrot = None

class Arguments(object):

    def __init__(self, parser):
        # This could raise in a few places. Let the caller handle these exceptions.
        self.parser = parser
        
        # Parse the arguments.
        self.parse()
        
        # Set all args.
        self.setAllArgs()
        
    def parse(self):
        try:
            self.parsedArgs = self.parser.parse_args()      
        except Exception as exc:
            if len(exc.args) == 2:
                type, msg = exc
                  
                if type != 'CmdlParser-ArgUnrecognized' and type != 'CmdlParser-ArgBadformat':
                    raise # Re-raise

                raise Exception('args', msg)
            else:
                raise # Re-raise

    def setArg(self, name, value):
        if not hasattr(self, name):
            # Since Arguments is a new-style class, it has a __dict__, so we can
            # set attributes directly in the Arguments instance.
            setattr(self, name, value)
        else:
            raise Exception('args', 'Attempt to set an argument that already exists: ' + name + '.')

    def setAllArgs(self):
        for key,val in list(vars(self.parsedArgs).items()):
            self.setArg(key, val)
        
    def getArg(self, name):
        try:
            return getattr(self, name)
        except AttributeError as exc:
            raise Exception('args', 'Unknown argument: ' + name + '.')



def getData(car,series):
    #getData() takes 4 different fits files from the disk at a given carrington rotaion and returns them
    global carrot 
    c = drms.Client()
    link = series + "[" + str(car) + "]"
    k, s = c.get(link, key=drms.const.all, seg='synop,chmap,logq,brq')
    #this code will only work on a singe record, but 's' may countain multiple records
    carrot = int(k.CAR_ROT.tolist()[0])
    synop = fits.open(s.loc[0,"synop"])
    chmap = fits.open(s.loc[0,"chmap"])
    logq = fits.open(s.loc[0,"logq"])
    brq = fits.open(s.loc[0,"brq"])
    return logq, brq,chmap,synop

def getImageData():
    #getImageData() takes the 3D numpy array data from the logq and brq files and transforms them into 1 slogq dataset
    brq_image = brqData[0].data
    logq_image = logqData[0].data
    brq_sign = np.sign(brq_image)
    Q = np.power(10,logq_image)
    Q[Q < 2] = 2
    logEq = (Q/2) + np.sqrt(((np.power(Q,2))/4)-1)
    slogq_image = brq_sign * np.log10(logEq)
    return slogq_image

def createImages(wr):
    #the first part here takes the carrington rotation and tranforms it into the day month and year that the carrot began
    global carrot
    JD = 27.2753*(carrot - 1.0) + 2398167.0
    t = Time(JD,format = 'jd')
    realDate = t.datetime.strftime('%B %d, %Y')
    #creates a manifest that the javascript and html reference in order to find the images
    manifest['slogfilenames'] = {}
    manifest['carrot'] = carrot
    manifest['date'] = realDate
    #saves the 10 heights of the slogq maps as seperate files labeled by their solar radii in the public_html/img/ directory
    colors = np.loadtxt(os.path.dirname(os.path.realpath(__file__)) + '/../data/qmap-colormap.txt')
    normalcolors = colors/255
    colormap_customized = mcol.ListedColormap(normalcolors)
    
    for slice in range(len(slogq_imageData[:,0,0])):
        plt.cla()
        slogq_slice = slogq_imageData[slice,:,:]
        absslogq_slice = np.abs(slogq_slice)
        colormax = np.nanmax(absslogq_slice)
        plt.imshow(slogq_slice, cmap = colormap_customized, origin = 'lower', extent = [0,360,-90,90], vmax = 3, vmin = -3)
        plt.title('S-logQ Map at %.3f Solar Radii' % radii[slice], y = 1.09)
        axisLines()
        SlogfileName = '/img/' + str(carrot) + 'slogqmap%.3f.html'  %radii[slice]
        with open(wr + SlogfileName, 'w') as mpld3image:
            mpld3.save_html(fig,mpld3image)   
        manifest['slogfilenames']["%.3f" %radii[slice]] = SlogfileName
        
    #saves the coronal hole map and the synoptic map to the same directory as the slogqmaps 
    plt.cla()
    chmap_slice = np.asarray(chmapData[1].data, dtype=np.float)
    plt.imshow(chmap_slice, origin = 'lower', extent = [0,360,-90,90], cmap = plt.cm.get_cmap('seismic_r'))
    plt.title('Coronal Hole Map', y = 1.09)
    axisLines()
    ChmapfileName = '/img/' + str(carrot) + 'chmap.html'
    with open(wr + ChmapfileName, 'w') as mpld3image:
        mpld3.save_html(fig,mpld3image)
    manifest['chfilename'] = ChmapfileName
        
    plt.cla()
    synop_slice = synopData[0].data
    plt.imshow(synop_slice, origin = 'lower', extent = [0,360,-90,90], cmap = plt.cm.get_cmap('Greys'))
    plt.title('Synoptic Map', y = 1.09)
    axisLines()
    SynopfileName = '/img/' + str(carrot) + 'synop.html'
    with open(wr + SynopfileName, 'w') as mpld3image:
        mpld3.save_html(fig,mpld3image)        
    manifest['synopfilename'] = SynopfileName
    
    #prints the manifest into a text file that is json serializable    
    jsonString = json.dumps(manifest)
    with open(wr + '/js/jsonString.js','w') as jsonFile:
        print(jsonString,file = jsonFile)
        
def axisLines():
    #creates the titles for the axis and the lines that surround each img
    plt.hlines(0,0,360,linewidth = .5)
    plt.vlines(360,-90,90,linewidth = 1)
    plt.hlines(90,0,360,linewidth = 1)
    plt.hlines(-90,0,360, linewidth = 1)
    plt.vlines(0,-90,90, linewidth = 1)
    ax.set_xlabel("Degrees Longitude(°)", y = 1.2)
    ax.set_ylabel("Degrees Latitude(°)", x = 1.2)
    
if __name__ == "__main__":
    #runs all of the functions
    error = False
    webRun = False
    
    if os.getenv('REQUEST_URI'):
        webRun = True
        jsonRoot = {}
        arguments = Arguments(None)
        
        try:
            # Try to get arguments with the cgi module. If that doesn't work, then fetch them from the command line.
            args = cgi.FieldStorage()
        
            if args:
                for key in args.keys():
                    val = args.getvalue(key)
                    if key in ('series'):
                        arguments.setArg('series', val)
                    elif key in ('carrot'):
                        arguments.setArg('carrot', int(val))
                    elif key in ('webroot'):
                        arguments.setArg('webroot', val)
                    else:
                        raise ValueError
        except ValueError:
            jsonRoot['status'] = JSON_RV_BADARGS
            jsonRoot['errMsg'] = 'Invalid CGI URL.\nUsage:\n  http://jsoc.stanford.edu/cgi-bin/qmapviewer?series=<DRMS qmap series>&carrot=<Carrington Rotation>&webroot=<web root directory>'
            error = True

    else:    
        parser = CmdlParser(usage='%(prog)s [ -h ] series=<DRMS series containing qmap images> carrot=<Carrington Rotation>')
        parser.add_argument('series', '--series', help='The DRMS data series that contains the input qmap images.', metavar='<qmap series>', dest='series', required=True)
        parser.add_argument('carrot', '--carrot', help='The Carrington Rotation for the qmap image.', metavar='<Carrington Rotation>', dest='carrot', type=int, required=True)
        parser.add_argument('webroot', '--webroot', help='The directory in which the images will be placed.', metavar = '<Web Root>', dest='webroot', required=True)
        arguments = Arguments(parser)

    if not error:
        try:
            wr = arguments.webroot
            car = arguments.carrot
            ds = arguments.series

            logqData, brqData, chmapData, synopData = getData(car,ds)
            slogq_imageData = getImageData()
            np.nanmax(slogq_imageData)
            fig, ax = plt.subplots(1,1)
            createImages(wr)

            if webRun:
                jsonRoot['status'] = JSON_RV_SUCCESS
        except:
            if webRun:
                jsonRoot['status'] = JSON_RV_CREATEIMGFAIL

    if webRun:
        print('Content-type: application/json\n')
        print(json.dumps(jsonRoot))

        sys.exit(0)
