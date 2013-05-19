function ars=hmi_noaa_info_interp(t,mode)
%hmi_noaa_info_interp	interpolate NOAA AR info at datetime t
% 
% ars=hmi_noaa_info_interp(t,mode)
% * Given a matlab datenum t, finds all NOAA AR information at 
% bracketing times t0 <= t < t1, and returns a list of structs 
% containing the AR info -- one struct per AR.
% * In fact, t0 = floor(t), and t1 = t0 + 1 day.
% * Typically the returned information will be interpolated.  
% Discrete quantities (like AR Zurich type) will use nearest-
% neighbor, and continuous ones will use linear interpolation 
% from the bracketing times.
% * If mode contains 'extrap', the information can be extrapolated
% for up to 3 days.  This is used for near-real-time results, and
% the latest time where information is present is automatically found.
% In this case, only longitudecm is altered, using a mean solar 
% rotation of the Carrington rate.  If the result is greater
% than 100 degrees, it is deleted from the list.
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
%   opt string mode = ''
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
if all(nargin  ~= [1 2]), error ('Bad input arg number'); end;
% if all(nargout ~= [0 1 2]), error ('Bad output arg number'); end;
if nargin < 2, mode = ''; end;

% allow extrapolation?
extrap = ~isempty(strfind(mode, 'extrap'));
if extrap, mode_raw = 'quiet'; else, mode_raw = ''; end;

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
[ars1,ok1] = hmi_noaa_info_raw(trec1, mode_raw);
[ars2,ok2] = hmi_noaa_info_raw(trec2, mode_raw);

% allow for empty AR lists
% FIXME: can/should we weaken this to && rather than || ?
% extrap_done: are we extrapolating beyond the last T_REC?
extrap_done = false;
extrap_max = 3; % max extrap time, in days
if ~ok1 || ~ok2,
  % no ARs found
  if extrap,
    % chance for extrapolation
    [arsX,okX] = hmi_noaa_info_raw('$', mode_raw);
    if okX && length(arsX) > 0, 
      tX = arsX(1).t;
    else,
      tX = NaN;
    end;
    if okX && (t > tX) && (t < (tX+extrap_max)),
      % ok, and t not too far ahead of noaa[$]: extrapolate to t
      ars1 = noaa_ar_extrap(t, arsX);
      ars2 = hmi_noaa_info_empty();
      % unique case where extrapolation is done
      extrap_done = true;
    else,
      % did not work
      okX = false; 
    end;
  else,
    % not allowed to extrapolate
    okX = false;
  end;
  % possible extrapolation failed too
  if ~okX,
    ars = hmi_noaa_info_empty();
    return;
  end;
end;

% keyboard

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
  if extrap_done,
    % extrapolating: always uses ars1 (see above)
    % note: the max time gap was set above
    wt1 = 1; wt2 = 0;
    ar1 = ars1(inx1);
    ar2 = ars1(inx1); % fake it
  elseif isempty(inx1),
    % not in earlier ar list
    wt1 = 0; wt2 = 1;
    ar1 = ars2(inx2);
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
    wt1 = dt2 / max(dt1 + dt2, eps);
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

% extrapolate arsX to time t
%   (expect t ahead of arsX)
function ars = noaa_ar_extrap(t, arsX)

  ars = arsX; % this almost does it all
  % eliminate edge condition
  if length(ars) == 0, return; end;
  % requires length > 0
  dt = t - ars(1).t; % in days, typically > 0
  % update longitudecm, allow for ar-dependent rate
  for i = 1:length(ars),
    rate = 360/27.2753; % degrees/day -- using Carrington rotation
    ars(i).longitudecm = ars(i).longitudecm + (dt * rate);
  end;
  % delete ARs that rotated off
  inx = ([ars.longitudecm] > 100);
  % note, this does leave field names intact, 
  % even if it empties the AR list out
  ars(inx) = [];
  return
  end
