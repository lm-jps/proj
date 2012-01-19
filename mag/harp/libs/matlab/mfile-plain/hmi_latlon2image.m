function [x,y,z]=hmi_latlon2image(lat,lon,geom)
%hmi_latlon2image	transform lat/lon coordinates to hmi image coordinates
% 
% [x,y,z]=hmi_latlon2image(lat,lon,geom)
% * Given (lat,lon) pairs and a disk geometry, transforms into image
% coordinates.
% * The lat and lon are assumed to be in degrees.
% * The solar tilt parameters B0 and P0 within geom are assumed to be 
% in degrees (as returned by hmidisk).
% 
% Inputs:
%   real lat(nf);
%   real lon(nf);
%   real geom(5)
% 
% Outputs:
%   real x(nf);
%   real y(nf);
%   real z(nf);
% 
% See Also:

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 31 Aug 2011
% Copyright (c) 2010.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [3]), error ('Bad input arg number'); end;
% if all(nargout ~= [2]), error ('Bad output arg number'); end;

%
% Computation
% 

% initialize sizes
nf = length(lat);
if length(lon) ~= length(lat),
  error('Incommensurate lat/lon lengths');
end;

% unpack geometry
xcen = geom(1);
ycen = geom(2);
rsun = geom(3);
beta = geom(4);
pang = geom(5);

% initial coordinates in (x,y,z)
% lat/lon in degrees
p1x = cosd(lat) .* sind(lon);
p1y = sind(lat);
p1z = cosd(lat) .* cosd(lon);

% column of (x,y,z) rows
p1 = [p1x(:) p1y(:) p1z(:)];

% rotation matrix: beta
%   sign check:  
%     when beta > 0, north pole should be visible (have z > 0)
%     north pole = NP = [0 1 0]
%     NP * Rbeta = [0 cos(beta) sin(beta)] -> z = sin(beta), correct.
Rbeta = zeros(3);
Rbeta(1,1) = 1; % beta (tip) leaves "x" fixed
Rbeta(2:3,2:3) = [cosd(beta) sind(beta); -sind(beta) cosd(beta)];

% rotation matrix: p
%   sign check:
%   when p > 0, north pole should move eastward (have x < 0)
%   NP * Rp = [-sin(p) cos(p) 0] -> x = -sin(p), correct.
Rp = zeros(3);
Rp(3,3) = 1; % p-angle (twist) leaves "z" fixed
Rp(1:2,1:2) = [cosd(pang) sind(pang); -sind(pang) cosd(pang)];

% rotate p0 into p1
p2 = p1 * (Rbeta * Rp);

% scale into image coordinates
% (note, xcen and ycen are already too large by 1)
p3 = rsun * p2 + repmat([xcen ycen 0], [nf 1]);

% unpack into x,y,z if desired
if nargout > 1,
  x = p3(:,1);
  y = p3(:,2);
  z = p3(:,3);
else,
  x = p3;
end;

return;
end
