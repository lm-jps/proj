function [image,err] = loadfits(file, block, type)
%loadfits	load array from file in FITS format (transposed)
%
% [image,err] = loadfits(file, block, type)
% * This is the transposed version of loadfitsa, equivalent to
% im = loadfitsa(file, block, type)';
% See that file for help.
%
% Inputs:
%   string file
%   opt int block = 1
%   opt string type = 'double'
%
% Outputs:
%   real image(...)
%   opt int err
%
% See Also: loadfitsa, savefits, savefits_old

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 19 Aug 2002
% Copyright (c) 2002.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [1:2]), error('Bad input arg number'); end;
% if all(nargout ~= [0 1]), error('Bad output arg number'); end;
if (nargin < 2), block = 1; end;
if (nargin < 3), type = 'double'; end;

%
% Computation
% 
% preserve nargout, because it affects error indications.
if nargout > 2,
  [imp,err] = loadfitsa(file, block, type);
else,
  imp = loadfitsa(file, block, type);
end;
 
% transpose (annoyingly expensive)
image = permute(imp, [ndims(imp):-1:1]); 
return;
