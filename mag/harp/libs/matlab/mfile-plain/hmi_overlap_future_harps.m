function [harpid,harps,crits]=hmi_overlap_future_harps(harp, harpids_now)
%hmi_overlap_future_harps	identify future HARPs overlapping a current one
%
% [harpid,harps,crits]=hmi_overlap_future_harps(harp, harpids_now)
% * Given a patch's trec (time) and geometry, load a set of future harps
% and identify which, if any, overlap the given patch.
% * Also requires the list of current track IDs, to filter them from the
% matching to future HARPs.  We already know no current track ID matches,
% because it did not turn up as a match in the more careful calculation 
% using bitmap overlap.
% * Plan is to call this from track_hmi_loop in gap-filling mode.
% * Input structure contains:
%      t_rec, a T_REC value
%      {lat,lon}dt{min,max}
% * Return value is the best harp id (or 0 if none matched), a struct of 
% future harp geometries, and a vector of overlap criteria (in [0,1]).
% 
% Inputs:
%   struct harp
%   int harpids_now(n)
% 
% Outputs:
%   int harpid
%   struct harps(m)
%   int crits(m)
% 
% See Also:  track_hmi_loop, hmi_gapfill_update_tid

% Written by Michael Turmon (turmon@jpl.nasa.gov) in 2013
% Copyright (c) 2013.  All rights reserved.  

% 
% Error checking
% 

if all(nargin  ~= [2]), error ('Bad input arg number'); end;
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end;

%
% Computation
% 
crit_good_enough = 0.8;
jsoc_time_format = 'yyyy.mm.dd_HH:MM_TAI';

% get time as a string, and as a datenum
if isreal(harp.t_rec),
  trec1_S = datestr(harp.t_rec, jsoc_time_format);
  trec1_R = harp.t_rec;
else,
  trec1_S = harp.t_rec;
  trec1_R = datenum(strrep(harp.t_rec, '_TAI', ''), 'yyyy.mm.dd_HH:MM:SS');
end;

trange = '4d@6h'; % could be run-dependent?
harps_all = hmi_compile_future_harps(trec1_S, trange);

% strip current harps from future harps
[junk,inx] = setdiff(harps_all.harpnum, harpids_now);
for f = fieldnames(harps_all)',
  fchar = f{1};
  x = harps_all.(fchar);
  x = x(:)'; % force into a row vector
  % restrict all the fields according to inx
  harps.(fchar) = x(inx);
end;
% note: future harp info now in "harps"

% time offset in days, per-harp (if in future, is > 0)
dt = harps.t_rec - trec1_R;
% longitude offset in degrees, per-harp (if in future, is > 0)
off = harps.omega_dt .* dt;
% lat/lons, offset
lat0 = harps.latdtmin;
lat1 = harps.latdtmax;
lon0 = harps.londtmin - off;
lon1 = harps.londtmax - off;

% area of new harp (scalar, deg^2)
A_harp = ...
    max(0, harp.latdtmax - harp.latdtmin) .* ...
    max(0, harp.londtmax - harp.londtmin);
% area of intersection: new_harp INT future_harps (vector)
% note: rectangle INT rectangle is still a rectangle
A_int = ...
    max(0, min(harp.latdtmax, lat1) - max(harp.latdtmin, lat0)) .* ...
    max(0, min(harp.londtmax, lon1) - max(harp.londtmin, lon0));
% overlap criterion, 0 <= crit <= 1, 1 is best (vector)
crits = A_int ./ A_harp;

% find the best, and make sure it's good enough
[crit_max, inx_max] = max(crits);
if crit_max >= crit_good_enough,
  harpid = harps.harpnum(inx_max);
else,
  harpid = 0;
end;

return;
end


