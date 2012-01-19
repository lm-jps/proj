function pn=hmipatchnum(trec)
%hmipatchnum	find patchnum parameter from HMI T_REC index
% 
% pn=hmipatchnum(trec)
% * Given a T_REC index, finds the number-of-patches parameter which lives
% in the PATCHNUM keyword.
% * A typical T_REC is: 2010.07.03_13:12:00_TAI
% * Given a list of T_REC's, returns the quality parameter for each one.
% The list of strings can be a cell array or a string array, with the
% strings "stacked vertically" in each case; see cellstr for example.
% If a cell array is used, it must be nf-by-1.
% 
% Inputs:
%   string trec(nf);        -- a valid time index
% 
% Outputs:
%   real pn(nf);
% 
% See Also:

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 17 Feb 2011
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
nf = size(trec,1);
pn = cell(nf,1); 

% try to find disk metadata in these data series
Parents = { 'hmi.Mpatch_720s' };

key_query = 'key=PATCHNUM';

% begin getting keywords
for f = 1:nf,
  for parent = Parents,
    % extra [1] extracts first patch
    query = sprintf('%s[%s][1]&%s', parent{1}, trec{f}, key_query);
    a = rs_list(query, 'web_access');
    if a.count > 0,
      break; % found the PATCHNUM
    end;
  end;
  if a.count == 0,
    error('Did not find quality bits for %s', trec{f});
  end;
  p1 = jsoc_cell2struct_keys(a.keywords);
  pnnum = p1.patchnum; % the string

  % insert the number(s)
  pn{f} = pnnum;
end;

% might want to reshape pn?
% or pair it with its trec's?

return;
end
