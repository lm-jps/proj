function err=hmi_save_track_local(rois, tid, birth, status, run_name, tmpl)
%hmi_save_track_local	save track information to disk
% 
% err=hmi_save_track_local(rois, tid, birth, status, run_name, tmpl)
% * Saves a set of rois representing a track.
% 
% Inputs:
%   struct rois(nt)
%   int tid
%   real birth
%   string status
%   string run_name
%   string tmpl
% 
% Outputs:
%   int err
% 
% See Also:

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 31 Aug 2010
% Copyright (c) 2010.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [6]), error ('Bad input arg number'); end;
%if all(nargout ~= [0 1 2 3 5]), error ('Bad output arg number'); end;

%
% Computation
% 
fprintf('closing track %d of %d chips [%s]\n', tid, length(rois), status);
% save roi-list to a mat-file
fn = sprintf(tmpl, 'mask', '-mask.', num2str(tid, '%06d'), 'mat');
[junk1,junk2] = mkdir(fileparts(fn)); % suppress warning if exists
save(fn, 'rois', 'tid');
err = 0;
return;
end
