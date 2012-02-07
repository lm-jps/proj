function [x0,y0,r,b0,p0]=hmidisk(trec,mode)
%hmidisk	find disk parameters from HMI T_REC index
% 
% [x0,y0,r,b0,p0]=hmidisk(trec,mode)
% * Given a T_REC index, finds the disk parameters which are assumed 
% to lie in the X0, Y0, and R_SUN keywords.  Returns these in "Matlab"
% coordinates, which means indexed from 1, so a center at the first
% pixel would be returned as x0 = y0 = 1.  (I.e., 1 is added to the
% x0,y0 keyword values.)
% * If five outputs are given, also finds SOLAR_B0 and SOLAR_P0, the
% solar tilt parameters.
% * This is done by querying the JSOC database.
% * A typical T_REC is: 2010.07.03_13:12:00_TAI .  If trec is given 
% as a matlab datenum (a double), this is converted to a T_REC string 
% before the JSOC database is queried.
% * Given a list of T_REC's or datenums, returns the disk parameters 
% for each one.  If given as T_REC's, the list of strings can be a 
% cell array or a string array, with the strings "stacked vertically" 
% in each case; see cellstr for example.  If a cell array of strings 
% is used, it must be nf-by-1.
% * If only one output is to be returned, it is returned as an nf-by-5
% matrix containing all values.  If mode = 'disk', only the first three
% disk parameters are returned.
%
% * There is a near-duplicate of this routine that could be harmonized
% with this routine; see hmi/wcs2center.m
% 
% Inputs:
%   string or real trec(nf);   -- a valid time index
%   opt char mode = 'geom'     -- 'geom' or 'disk'
% 
% Outputs:
%   real x0(nf);
%   real y0(nf);
%   real r(nf);
%   real b0(nf);
%   real p0(nf);
% 
% See Also:

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 31 Aug 2010
% Copyright (c) 2010.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [1 2]), error ('Bad input arg number'); end;
if all(nargout ~= [0 1 2 3 5]), error ('Bad output arg number'); end;
if nargin < 2, mode = 'geom'; end;

%
% Computation
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WCS transformations

% Nested functions for WCS transformations.  
% ASSUMES: 
%   crpix1, crpix2 = CRPIX1, CRPIX2
%   sina,cosa = sin and cos of CROTA2 resp.
%   crvalx and crvaly are CRVAL1 and CRVAL2, 
%   cdelt = CDELT1 == CDELT2
% PIX_X and PIX_Y are CCD pixel addresses, 
% WX and WY are arc-sec W and N on the Sun from disk center.

DTOR = pi/180;
RAD2ARCSEC = 648000/pi;

% nested function
function px = PIX_X(wx,wy)
px = ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1);
return;
end

% nested function
function py = PIX_Y(wx,wy) 
py = ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2);
return;
end

% (end WCS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert input from double (datestr) to char (t_rec)
if isnumeric(trec),
  trec = datestr(trec, 'yyyy.mm.dd_HH:MM:SS_TAI');
end;

% convert char input to cell array
if ischar(trec),
  trec = cellstr(trec);
end;

% initialize sizes
nf = size(trec,1);
x0 = zeros(nf,1); 
y0 = x0;
r  = x0;
b0 = x0;
p0 = x0;

% try to find disk metadata in these data series
Parents = { 'hmi.M_720s', 'hmi.M_720s_nrt' };

% keys needed for disk params
keys = {'CRVAL1', ... % disc center in arcsec
        'CRVAL2', ...
        'CDELT1', ... % arcsec, assumimg dx=dy
        'CRPIX1', ... % disk center on ccd
        'CRPIX2', ...
        'CROTA2', ... % rotation
        'CRLT_OBS', ... % b angle
        'RSUN_REF', ...
        'DSUN_OBS'};

key_query = ['key=' keys{1}];
for k = keys(2:end),
  key_query = strcat(key_query, ',', k{1});
end;

% begin getting keywords
for f = 1:nf,
  for parent = Parents,
    query = sprintf('%s[%s]&%s', parent{1}, trec{f}, key_query);
    a = rs_list(query, 'web_access');
    if a.count > 0,
      break; % found the WCS
    end;
  end;
  if a.count == 0,
    error('Did not find disk metadata for %s', trec{f});
  end;
  wcs = jsoc_cell2struct_keys(a.keywords);

  % this is not comprehensive, but I think it catches all cases
  if iscell(wcs.crota2) || iscell(wcs.crpix1),
    error('Missing values in WCS keywords using %s', query);
  end;

  % Ephemeris courtesy of Phil's (nested functions) above
  % I have preserved the somewhat over-wrought nature of these
  % functions so that these lines are almost the same as the
  % original C macros I copied.
  crvalx   = wcs.crval1;
  crvaly   = wcs.crval2;
  cdelt    = wcs.cdelt1;
  crpix1   = wcs.crpix1;
  crpix2   = wcs.crpix2;
  crota2   = wcs.crota2;
  dsun_obs = wcs.dsun_obs;
  rsun_ref = wcs.rsun_ref;
  % rsun_ref = 6.96e8;

  sina = sin(crota2 * DTOR); 
  cosa = cos(crota2 * DTOR); 

  Vx0 = PIX_X(0.0,0.0) - 1.0; % zero-based coordinates
  Vy0 = PIX_Y(0.0,0.0) - 1.0;
  Vr  = asin(rsun_ref / dsun_obs) * RAD2ARCSEC / cdelt;
  Vb0 = wcs.crlt_obs; % sign is correct
  Vp0 = -crota2;
    
  % insert 
  x0(f) = Vx0; y0(f) = Vy0; r(f) = Vr; b0(f) = Vb0; p0(f) = Vp0;
end;
x0 = x0+1; y0 = y0+1; % convert to Matlab coordinates
% combine into a disk or geom if it looks right
if nargout <= 1,
  if strcmp(mode, 'geom'),
    x0 = [x0 y0 r b0 p0];
  else,
    x0 = [x0 y0 r];
  end;
end;
return;
end
