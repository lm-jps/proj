function varargout=hmidisk(trec,mode)
%hmidisk	find disk parameters from HMI T_REC index
% 
% geom = hmidisk(trec,mode)            -- vector, error
% [geom,msg] = hmidisk(trec,mode)      -- vector, quiet
% [x0,y0,r,b0,p0] = hmidisk(trec,mode) -- geom, scalar
% [x0,y0,r] = hmidisk(trec,mode)       -- disk, scalar
% * Given a T_REC index and an output mode, finds the disk parameters
% by querying JSOC.
% * A typical T_REC is: 2010.07.03_13:12:00_TAI .  If trec is given 
% as a matlab datenum (a double), this is converted to a T_REC string 
% before the JSOC database is queried.
% * Several output modes are allowed.  
% * The first is `geom', `disk', or `wcs'.  There are 5 `geom' parameters, 
% the position stored in the X0, Y0, and R_SUN keywords, and the tilt
% angles stored in the SOLAR_B0 and SOLAR_P0 keywords.  Say 'geom' for
% all 5, 'disk' for just the first 3, or 'wcs' for a structure that 
% contains all the WCS fields.
% * The second is 'vector' vs. 'scalar'.  If vector is used, a 1x3 or
% 1x5 vector is returned, else, corresponding scalars are returned.
% (This need not be given for 'wcs'.)
% * The last is 'error' vs. 'quiet'.  If 'error' is given, DB query
% errors cause Matlab errors.  Otherwise, a descriptive message is returned
% as the fourth or sixth output.  If the message is empty, the query succeeded.
% * If mode is not given, we use: geom,vector,error.  But, if anything
% is specified, all must be specified.
% * X0 and Y0 are returned in Matlab coordinates, which means indexed 
% from 1, so a center at the first pixel would be returned as 
% x0 = y0 = 1.  (Specifically, 1 is added to the x0,y0 keyword values.)
% 
% Inputs:
%   string or real trec;  -- a valid time index
%   opt string mode = 'geom,vector,error'
%                         -- {geom, disk, wcs} X {vector, scalar} X {error,quiet}
% 
% Outputs:
%   real x0
%   real y0
%   real r
%   real b0
%   real p0
%  -or-
%   real geom(3) or (5)
%  -or-
%   struct wcs
%  -plus optionally-
%   string msg
% 
% See Also: wcs2center

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 31 Aug 2010
% Copyright (c) 2010.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [1 2]), error ('Bad input arg number'); end;
if all(nargout ~= [0:6]), error ('Bad output arg number'); end;
if nargin < 2, mode = 'geom,vector,error'; end;

% set up output mode
if ~isempty(strfind(mode, 'geom')),
  geom_mode = 1;
elseif ~isempty(strfind(mode, 'disk')),
  geom_mode = 0;
elseif ~isempty(strfind(mode, 'wcs')),
  geom_mode = 2;
else,
  error('Need valid output mode (geom/disk/wcs)');
end;

% set up output mode
if geom_mode == 2,
  % (not actually used)
  N_out = 1;
elseif ~isempty(strfind(mode, 'vector')),
  vec_mode = 1;
  N_out = 1;
elseif ~isempty(strfind(mode, 'scalar')),
  vec_mode = 0;
  N_out = 3 + (geom_mode == 1) * 2;
else,
  error('Need valid output mode (vector/scalar)');
end;

% set up error handling
if ~isempty(strfind(mode, 'error')),
  do_err_out = 1;
elseif ~isempty(strfind(mode, 'quiet')),
  do_err_out = 0;
else,
  error('Need valid error-handling mode (error/quiet)');
end;

% outputs, each initialized to []
varargout = cell(1, N_out + (1-do_err_out));
% check that the numbers agree; if no outs asked for, let it go
if nargout > 0 && nargout ~= length(varargout),
  error('Mode did not match number of output arguments requested');
end;

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

% look for disk metadata in this data series
nrt_mode = hmi_property('get', 'nrt_mode');
if nrt_mode,
  parent = 'hmi.M_720s_nrt';
else,
  parent = 'hmi.M_720s';
end

% keys needed for disk params
keys = {'CRVAL1', ... % disc center in arcsec
        'CRVAL2', ...
        'CDELT1', ... % arcsec, assumimg dx=dy
        'CRPIX1', ... % disk center on ccd
        'CRPIX2', ...
        'CROTA2', ... % rotation
        'CRLT_OBS', ... % b angle
        'CRLN_OBS', ... 
        'RSUN_REF', ...
        'DSUN_OBS'};

keylist = sprintf('%s,', keys{:});
keylist(end) = []; % kill last comma

% get keywords
%  note: the camera>0 clause guards against duplicate entries in M_720s
%  with missing data -- it forces only valid records to be returned.
series = sprintf('%s[%s][?camera>0?]', parent, trec);
if do_err_out,
  % do not ask for msg argument -- will error out if WCS is not present
  a = rs_list(series, {'key', keylist});
else,
  % ask for msg argument
  [a,msg] = rs_list(series, {'key', keylist});
  if ~isempty(msg),
    varargout{end} = msg;
    return;
  end;
end;

% insist on exactly one response
if a.count ~= 1,
  msg = sprintf('Response for %s for WCS had count = %d, needed 1', series, a.count);
  if do_err_out,
    error(msg);
  else,
    varargout{end} = msg;
    return;
  end;
end;

% get wcs from key structure
wcs = jsoc_cell2struct_keys(a.keywords);

if ~isfield(wcs, 'crota2') || ~isfield(wcs, 'crpix1'),
  msg = sprintf('Keyword response for %s for WCS not present', series);
  if do_err_out,
    error(msg);
  else,
    varargout{end} = msg; % record error
    return;
  end;
end;
% this is not comprehensive, but I think it catches all cases
if iscell(wcs.crota2) || iscell(wcs.crpix1),
  msg = sprintf('Missing values in WCS keywords using %s', series);
  if do_err_out,
    error(msg);
  else,
    varargout{end} = msg; % record error
    return;
  end;
end;

% can leave this routine now
if geom_mode == 2,
  varargout{1} = wcs;
  return;
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
    
% set outputs (convert x,y to Matlab coordinates)
x0 = Vx0 + 1; 
y0 = Vy0 + 1; 
r  = Vr; 
b0 = Vb0; 
p0 = Vp0;

% plug outputs in
if geom_mode == 1,
  if vec_mode,
    varargout{1} = [x0 y0 r b0 p0];
  else,
    varargout{1} = x0;
    varargout{2} = y0;
    varargout{3} = r;
    varargout{4} = b0;
    varargout{5} = p0;
  end;
else,
  if vec_mode,
    varargout{1} = [x0 y0 r];
  else,
    varargout{1} = x0;
    varargout{2} = y0;
    varargout{3} = r;
  end;
end;

return;
end
