function ars=hmi_noaa_info_raw(trec)
%hmi_noaa_info_raw	find NOAA AR numbers from HMI T_REC index
% 
% ars=hmi_noaa_info_raw(trec)
% * Given a T_REC index, finds all NOAA AR information and returns
% a list of structs containing the AR info -- one struct per AR per 
% day.  
% * Note that the NOAA series uses UT (Zulu time), and are sampled 
% at 23:59:59 on the given day, so it can be useful to use queries 
% like "YYYY.MM.DD_UT/1d" to get the AR info for a given day.
% * Could be extended to take Matlab datetime's.
% * Given a list of T_REC's, returns the AR parameter for each one.
% The list of strings can be a cell array or a string array, with the
% strings "stacked vertically" in each case; see cellstr for example.
% If a cell array is used, it should be nf-by-1.
% 
% Inputs:
%   string trec(nf);        -- a valid time index
% 
% Outputs:
%   struct ars(nf);
%   bool ok(nf);
% 
% See Also: hmi_noaa_info_interp

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

% convert input to cell array
if ischar(trec),
  trec = cellstr(trec);
end;

% initialize sizes
nf = length(trec);

% source for NOAA metadata, maintained by Rick's scripts
source = 'su_rsb.NOAA_ActiveRegions';
% which keys to get
key_query = 'key=**ALL**';
% date string format (omit trailing 'Z')
date_fmt = 'yyyy.mm.dd_HH:MM:SS';

% begin getting keywords
ars = [];
for f = 1:nf,
  query = sprintf('%s[%s]&%s', source, trec{f}, key_query);
  a = rs_list(query, 'web_access');
  if a.count == 0,
    ar1 = hmi_noaa_info_empty();
    warning('No NOAA info for %s, leaving empty, continuing.', trec{f});
  else,
    % does some of the conversion work
    s1 = jsoc_cell2struct_keys(a.keywords);
    % put into a list of structs
    ar1 = [];
    for i = 1:length(s1.regionnumber),
      for f = fieldnames(s1)',
        f1 = f{1};
        x = s1.(f1)(i);
        if iscell(x), x = x{1}; end; % unbox strings
        ar1(i).(f1) = x;
      end;
    end;
  end;
  ars = cat(2, ars, ar1);
end;

% add in the obs-time as a matlab datenum
for i = 1:length(ars),
  t1 = datenum(ars(i).observationtime(1:end-1), date_fmt);
  % record this time as TAI -- it came in Z = UTC
  ars(i).t = t1 + 34/(3600*24);
end;

return;
end


    