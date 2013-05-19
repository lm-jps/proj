function tids_unused=hmi_unused_harpnums(tid_now, tids)
%hmi_gapfill_update_tid	produce updated track id of bogey track
% 
% tids_unused=hmi_unused_harpnums(tid_now, tids)
% * Returns a list of so-far-unused track IDs (HARPNUM's).  These HARPNUMs
% are "behind" some given current track, tid_now, and moreover do not 
% overlap the currently-known track IDs (tids).
% * FIXME: this does not account for tracks that have been written out
% to intermediate files, but not ingested.
% 
% Inputs:
%   int tid_now
%   int tids(n)
% 
% Outputs:
%   int tids_unused(m)
% 
% See Also: 

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 13 Mar 2013
% Copyright (c) 2013.  All rights reserved.

% 
% Error checking
% 
% if all(nargin  ~= [4]), error ('Bad input arg number'); end;
%if all(nargout ~= [0 1 2 3 5]), error ('Bad output arg number'); end;

%
% Computation
% 
gapwidth = 50;
% pool of possible track id's
%  (includes all tracks in above query except tid_now)
tids_pool = [max(1,tid_now-gapwidth):max(1,tid_now-1)];
if isempty(tids_pool),
  % not an error
  tids_unused = [];
  return;
end;

% determine which tracks have been ingested already
query = sprintf('hmi.Mharp_720s[%d-%d][][?T_REC=T_FRST1?]', tids_pool(1), tids_pool(end));
[res,msg] = rs_list(query, {'key','HARPNUM'});
if ~isempty(msg),
  fprintf('%s: Trouble: JSOC query to get unused HARPs failed (%s).\n', mfilename, msg);
  fprintf('%s: Proceeding without unused HARPs.\n', mfilename);
  tids_unused = [];
  return;
end;
if res.count < 0,
  fprintf('%s: Trouble: JSOC query to get unused HARPs failed.\n', mfilename);
  fprintf('%s: Proceeding without unused HARPs.\n', mfilename);
  tids_unused = [];
  return;
end;
% tids_used: have used these track IDs already
if res.count > 0,
  tids_used = str2double(res.keywords{1}.values);
else,
  % found no used tids: res.keywords will not exist
  tids_used = [];
end;

% remove already-ingested and in-progress tracks
tids_unused = setdiff(tids_pool, tids_used);
tids_unused = setdiff(tids_unused, tids);

% FIXME: need to also setdiff out tracks that have been written to
% ingestor-ready files, but not actually ingested.

% take largest first
tids_unused = sort(tids_unused, 'descend');

return;
end

