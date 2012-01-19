function hmi_init_run_local(mode, hooks, rgn)
%hmi_init_run_local	tracker file-output begin/end hook
% 
% hmi_init_run_local(mode, hooks, rgn)
% * Initializes (or finalizes) a tracker run, depending on 
% the sign of the mode input.
% * This version is for debugging purposes at JPL.
% 
% Inputs:
%   (see hmi_init_run_common)
% 
% Outputs:
%   (to filesystem)
% 
% See Also:

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 31 Dec 2010
% Copyright (c) 2010.  All rights reserved.

scenario = 'diagnostic';
hmi_init_run_common(scenario, mode, hooks, rgn);
return;
end

