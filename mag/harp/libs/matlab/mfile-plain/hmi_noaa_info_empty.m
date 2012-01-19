function ar=hmi_noaa_info_empty()
%hmi_noaa_info_empty	empty struct for noaa ARs
% 
% ar=hmi_noaa_info_empty()
% * Return empty (0x0) struct containing correct fields for NOAA AR.
% 
% Inputs:
%   (none)
% 
% Outputs:
%   struct ars(0);
% 
% See Also: hmi_noaa_info_interp

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 26 Sep 2011
% Copyright (c) 2011.  All rights reserved.

% 
% Error checking
% 
% if all(nargin  ~= [1]), error ('Bad input arg number'); end;
% if all(nargout ~= [0 1 2]), error ('Bad output arg number'); end;

%
% Computation
% 
ar = struct('observationtime', {}, ...
            'regionnumber', {}, ...
            'zurichclass', {}, ...
            'magnetictype', {}, ...
            'spotcount', {}, ...
            'area', {}, ...
            'latitudehg', {}, ...
            'longitudehg', {}, ...
            'longitudecm', {}, ...
            'longitudinalextent', {}, ...
            't', {});
return;
end


