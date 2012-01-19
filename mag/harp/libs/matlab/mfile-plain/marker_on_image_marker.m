function [x,y]=marker_on_image_marker(name, scale, thick)
%marker_on_image_marker	return mask for image marker given by name
% 
% [x,y]=marker_on_image_marker(name, scale, thick)
% * Return a list of point locations, (x,y), for a marker given
% by `name' at the given scale.  Valid names are:
%  + = cross
%  x = ex
%  * = star
%  s = square
% * The point locations are confined to the box [1..scale] x [1..scale].
% You can transform them to the coordinates of a given image im using:
%   inx = sub2ind(size(im), x, y);
% and then assign the given locations to a value by: 
%   im(inx) = 255
% If you need the marker to be centered about some (x0,y0) point, 
% you will need to add (x0,y0) and subtract off (scale-1)/2.
% * If thick is given, it is used to widen the marker shape.
% The box size stays the same.
% 
% Inputs:
%  string name
%  opt int scale = 11
%  opt int thick = 1
% 
% Outputs:
%  int x(n)
%  int y(n)
% 
% See Also:  marker_on_image

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 25 Sep 09.
% Copyright (c) 2009.  All rights reserved.

% Error checking
% 
if all(nargin  ~= [1 2 3]), error ('Bad input arg number'); end
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end  
if nargin < 2, scale = 11; end;
if nargin < 3, thick = 1; end;

%
% Computation
% 

% index sets
Npt = 21;
inx  = linspace(0, 1, Npt);
inx0 = zeros(size(inx));
inx1 = ones (size(inx));

if strcmp(name, '+'),
  % + symbol
  % draw order: horizontal, vertical
  rpts = [0.5*inx1 inx     ];
  cpts = [inx      0.5*inx1];
elseif strcmp(name, 's'),
  % square symbol
  % draw order: bottom, right, top, left
  rpts = [inx0 inx  inx1  inx ];
  cpts = [inx  inx1 1-inx inx0];
elseif strcmp(name, 'x'),
  % x
  % draw order: SW->NE, NW->SE
  rpts = [inx inx];
  cpts = [inx 1-inx-eps]; % avoid doubled line
elseif strcmp(name, '.'),
  % .
  % draw order: SW->NE, NW->SE
  rpts = [0.5];
  cpts = [0.5];
elseif strcmp(name, '*'),
  % * (four lines)
  % first, do +: draw order: horizontal, vertical
  rpts = [0.5*inx1 inx     ];
  cpts = [inx      0.5*inx1];
  % then, do x: draw order: SW->NE, NW->SE
  rpts = [rpts inx inx];
  cpts = [cpts inx 1-inx-eps]; % avoid doubled line
  % finally, strip extreme points
  edges = hypot(rpts-0.5, cpts-0.5) > 0.6;
  rpts(edges) = [];
  cpts(edges) = [];
else,
  error('Did not recognize symbol %s', name);
end;

% transform by scale
x = round(1 + (scale-1) * cpts);
y = round(1 + (scale-1) * rpts);

% clip to box [1..scale] x [1..scale]
bad = (x < 1) | (x > scale) | (y < 1) | (y > scale);
x(bad) = [];
y(bad) = [];

% widen by thick 
ker = ones(1, thick);
% round-trip to image
% (filters identical indexes, even if ker = 1)
z = zeros(scale, scale);
inx = sub2ind(size(z), x, y);
z(inx) = 1;
z = conv2(ker, ker, z, 'same');
[xt,yt] = find(z);
x = xt';
y = yt';

return;
end

