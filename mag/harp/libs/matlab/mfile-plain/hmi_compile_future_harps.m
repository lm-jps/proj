function harps=hmi_compile_future_harps(trec, trange)
%hmi_compile_future_harps	compile a list of future HARPs and their geometry
%
% harps=hmi_compile_future_harps(trec, trange)
% * Use given trec (time) and trange to load HARP metadata at a given time
% so that metadata can be saved and a HARP/NOAA match made.
% * Plan is to call this from track_hmi_loop in gap-filling mode.
% * Input trec is a T_REC value, e.g., corresponding to a mask image. 
% * Input trange is an interval, e.g., 4d@6h
% * The output is a structure that contains:
%   harpnum, t_rec, omegadt, {lat,lon}dt{min,max}
% 
% Inputs:
%   string trec;    -- trec for HARP series
%   string trange;  -- matlab time number
% 
% Outputs:
%   int res
% 
% See Also:  track_hmi_loop

% Written by Michael Turmon (turmon@jpl.nasa.gov) in 2013
% Copyright (c) 2010.  All rights reserved.  

% 
% Error checking
% 

if all(nargin  ~= [2]), error ('Bad input arg number'); end;
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end;

%
% Computation
% 

% fields in returned structure
fields = {'t_rec',    ...
          'harpnum',  ...
          'latdtmin', ...
          'londtmin', ...
          'latdtmax', ...
          'londtmax', ...
          'omega_dt', ...
         };

% parameters to get from each HARP
fields_comma = sprintf('%s,', fields{:});
params = {'key', fields_comma(1:end-1)}; % strip trailing ,

ds_harp = 'hmi.Mharp_720s';

% get HARP numbers and HARP metadata
series = sprintf('%s[][%s/%s]', ds_harp, trec, trange);
[reply,msg] = rs_list(series, params);

% make empty harps in case of early return
Cs = [fields; cell(size(fields))]; % alternate fieldname with []
harps = struct(Cs{:});
if ~isempty(msg),
  return;
end;
if reply.count == 0,
  return;
end;

harps_all = jsoc_cell2struct_keys(reply.keywords);

% de-duplicate -- don't use fancy R2012 features
[junk,inx] = unique(harps_all.harpnum);

% we take all the fields
for f = fieldnames(harps_all)',
  fchar = f{1};
  x = harps_all.(fchar);
  x = x(:)'; % force into a row vector
  harps.(fchar) = x(inx);
end;

return;
end


