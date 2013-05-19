function [tid_new,t_add]=hmi_gapfill_update_tid(tid, tids, currtime, stats)
%hmi_gapfill_update_tid	produce updated track id of bogey track
% 
% [tid_new,t_add]=hmi_gapfill_update_tid(tid, tids, currtime, stats)
% * Returns the updated track ID of a bogey track.  This is a track id
% which tags a track structure that was just created, and we're
% not sure if it corresponds to some future track, or a new track.
% Upon creation, it is coded as a tid field with: ROI_t(inx).tid < 0.
% * Also returns a 0/1 variable saying whether the track count,
% ROI_rgn.Ntot, should be incremented (1 if it should).  If tid_new
% is a recycled track ID, Nt should not be incremented.
% * The usage would be:
%    [tid_new,t_add] = hmi_gapfill_update_tid(...)
%    ROI_t(inx).tid = tid_new;
%    ROI_rgn.Ntot = ROI_rgn.Ntot + t_add;
% * There are three cases:
%  -- the track does correspond to a future track; then we return
%     the track id of the future track
%  -- the track doesn't correspond to a future track, but some 
%     recycled track IDs exist; then we return such a recycled ID.
%  -- otherwise, we negate the given track ID and use it.  In gap-filling
%     mode, this is an error condition, because it means post-gap tracks 
%     will be offset by one.
% 
% Inputs:
%   int tid
%   int tids(n)
%   real currtime -- datenum
%   real stats(1,NS)
% 
% Outputs:
%   int tid_new
%   int t_add    -- 0/1
% 
% See Also: track_hmi_loop

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 13 Mar 2013
% Copyright (c) 2013.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [4]), error ('Bad input arg number'); end;
%if all(nargout ~= [0 1 2 3 5]), error ('Bad output arg number'); end;

%
% Computation
% 

% This updates "bogey" tracks by finding the correct track id.  It cannot
% be done by track_hmi_loop>tracker_new_track because the location of the 
% ROI is not known when tracker_new_track is called.  
% So, that routine assigns a "bogey" track id (< 0), which is updated by 
% this routine based on future HARPs.
% This only makes sense when filling in gaps.

% from roi_stats_mag_defs.h
RS_rgn_min_lat = 4; % lower corner of (lat,lon) bounding box for region
RS_rgn_min_lon = 5;
RS_rgn_max_lat = 6; % upper corner of (lat,lon) bounding box for region  
RS_rgn_max_lon = 7;

% fill in harp.*
%   need time and lat/lon bounding box
harp.t_rec    = currtime;
harp.latdtmin = stats(RS_rgn_min_lat);
harp.londtmin = stats(RS_rgn_min_lon);
harp.latdtmax = stats(RS_rgn_max_lat);
harp.londtmax = stats(RS_rgn_max_lon);

fprintf('  Gapfill: Updating bogey TID = %d\n', tid);

% does a future harp overlap the current harp?
%   (returns best match, or 0 if no acceptable match)
[harpid,harps,crits] = hmi_overlap_future_harps(harp, tids);
if harpid > 0,
  fprintf('  Gapfill: Got a matching TID = %d (crit = %.2f)\n', ...
          harpid, max(crits));
  tid_new = harpid;
  % this track looks new to us, but it was already seen in the gappy run,
  % so it's OK to add to the track count
  t_add = 1; 
  return;
end;

% otherwise: can we recycle a past harpnum that was unused?
%  Note: abs(tid) is the non-bogey value
% FIXME: can anything in tids() be negative?
% What if there is more than one bogey at a given time?
tids_unused = hmi_unused_harpnums(abs(tid), tids(tids>0));
if ~isempty(tids_unused),
  fprintf('  Gapfill: Got %d recyled TIDs, using %d\n', ...
          length(tids_unused), tids_unused(1));
  % (it returns a list of all of them)
  tid_new = tids_unused(1);
  % this should not count as a new track; we want to hide its 
  % presence, so don't disturb the track count
  t_add = 0; 
  return;
end;

% otherwise: negate the given track id
%  it will be as if a new regular track was allocated
% Note, in gapfill mode, this is an error condition, because
% we can't be adding new track IDs.
tid_new = -1 * tid;
t_add = 1;

fprintf('  Gapfill: ERROR: had to update as regular track, TID = %d\n', tid_new);

return;
end

