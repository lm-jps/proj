function hmi_init_run_production(mode, hooks, rgn)
%hmi_init_run_production	tracker file-output begin/end hook
% 
% hmi_init_run_production(mode, hooks, rgn)
% * Initializes (or finalizes) a tracker run, depending on 
% the sign of the mode input.
% * This version is for HMI production, so it produces metafiles 
% for ingestion into JSOC.
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

scenario = 'jsoc,diagnostic';
hmi_init_run_common(scenario, mode, hooks, rgn);
return;
end
