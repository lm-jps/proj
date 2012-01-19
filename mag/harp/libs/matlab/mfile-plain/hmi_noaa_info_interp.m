function ars=hmi_noaa_info_interp(t)
%hmi_noaa_info_interp	interpolate NOAA AR info at datetime t
% 
% ars=hmi_noaa_info_interp(t)
% * Given a matlab datenum t, finds all NOAA AR information at 
% bracketing times t0 <= t < t1, and returns a list of structs 
% containing the AR info -- one struct per AR.
% * In fact, t0 = floor(t), and t1 = t0 + 1 day.
% * Typically the information will be interpolated.  Discrete
% quantities (like AR Zurich type) will use nearest-neighbor, 
% and continuous ones will use linear interpolation from the 
% bracketing times.
% * The fields returned are:
%    t                  -- a datenum
%    observationtime    -- trec string for t, in TAI
%    regionnumber       -- NOAA ID
%    zurichclass        -- string
%    magnetictype       -- string
%    spotcount          -- integer
%    area
%    latitudehg
%    longitudehg
%    longitudecm
%    longitudinalextent
% 
% Inputs:
%   real t;        -- a datenum
% 
% Outputs:
%   struct ars(nr);
% 
% See Also: hmi_noaa_info_raw

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 26 Sep 2011
% Copyright (c) 2011.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [1]), error ('Bad input arg number'); end;
% if all(nargout ~= [0 1 2]), error ('Bad output arg number'); end;

%
% Computation
% 
date_fmt = 'yyyy.mm.dd'; % date part 
trec_fmt = '%s_TAI/1d'; % query for a 1-day range in TAI
% get bracketing times
t1 = floor(t); % nearest day, rounding down
trec1 = sprintf(trec_fmt, datestr(t1,   date_fmt));
trec2 = sprintf(trec_fmt, datestr(t1+1, date_fmt));
% these indirectly query JSOC
ars1 = hmi_noaa_info_raw(trec1);
ars2 = hmi_noaa_info_raw(trec2);

% allow for empty AR lists
if isempty(ars1) && isempty(ars2),
  % no ARs -- not an error
  % ars = struct([]); % should be an empty struct?
  ars = hmi_noaa_info_empty();
  return;
end;

% fields to interpolate (others are manual, or nearest-neighbor)
interpolated_fields = {'area', ...
                    'latitudehg', 'longitudehg', ...
                    'longitudecm', 'longitudinalextent'};

% cache region IDs
ids1 = [ars1.regionnumber];
ids2 = [ars2.regionnumber];
ids = union(ids1, ids2); % unique AR id's, sorted
Nid = length(ids);

% set up each AR in result, 1:Nid
ars = hmi_noaa_info_empty();
trec = sprintf('%s_TAI', datestr(t, 'yyyy.mm.dd_HH:MM:SS'));
N_ar = 1; % number of ARs created so far
for i = 1:Nid,
  id = ids(i);
  % one of the below may be empty, but not both
  inx1 = find(ids1 == id);
  inx2 = find(ids2 == id);
  if isempty(inx1),
    % not in earlier ar list
    wt1 = 0; wt2 = 1;
    ar1 = ars2(inx2); % fake it (TODO: extrapolate longitudecm?)
    ar2 = ars2(inx2);
    if abs(t - ar1.t) > 0.25, 
      continue; % too far away to fake
    end;
  elseif isempty(inx2),
    % not in later ar list
    wt1 = 1; wt2 = 0;
    ar1 = ars1(inx1);
    ar2 = ars1(inx1); % fake it
    if abs(t - ar2.t) > 0.25, 
      continue; % too far away to fake
    end;
  else,
    % typical case: appears in both
    dt1 = t - ars1(inx1).t; % >= 0
    dt2 = ars2(inx2).t - t; % >= 0
    % allow for dt1 == dt2 == 0, even though it should not happen
    wt1 = dt2 / (dt1 + dt2 + eps);
    wt2 = 1 - wt1;
    clear dt1 dt2
    ar1 = ars1(inx1);
    ar2 = ars2(inx2);
  end;
  % nearest neighbor is default
  if wt1 > wt2,
    ars(N_ar) = ar1;
  else,
    ars(N_ar) = ar2;
  end;
  % but, interpolate some fields
  for f = interpolated_fields,
    f1 = f{1};
    x = wt1 * ar1.(f1) + wt2 * ar2.(f1);
    ars(N_ar).(f1) = x;
  end;
  % finally, set up time
  ars(N_ar).t = t;
  ars(N_ar).observationtime = trec;
  % this finishes the AR
  N_ar = N_ar + 1;
end;
return;
end
