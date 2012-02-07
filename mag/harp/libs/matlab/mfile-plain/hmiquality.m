function [q,ok]=hmiquality(trec)
%hmiquality	find quality parameter from HMI T_REC index
% 
% [q,ok]=hmiquality(trec)
% * Given a T_REC index, finds the quality parameter which lives
% in the QUALITY keyword.
% * If two outputs are given, also finds ok, a boolean value which
% tells if the quality is sufficient for most operations.
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
%   real q(nf);
%   bool ok(nf);
% 
% See Also:

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 17 Feb 2011
% Copyright (c) 2011.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [1]), error ('Bad input arg number'); end;
if all(nargout ~= [0 1 2]), error ('Bad output arg number'); end;

%
% Computation
% 

% convert input to cell array
if ischar(trec),
  trec = cellstr(trec);
end;

% initialize sizes
nf = size(trec,1);
q  = zeros(nf,1); 
ok = false(nf,1);

% try to find disk metadata in these data series
Parents = { 'hmi.M_720s', 'hmi.M_720s_nrt' };

key_query = 'key=QUALITY';

% begin getting keywords
for f = 1:nf,
  for parent = Parents,
    query = sprintf('%s[%s]&%s', parent{1}, trec{f}, key_query);
    a = rs_list(query, 'web_access');
    if a.count > 0,
      break; % found the QUALITY
    end;
  end;
  if a.count == 0,
    error('Did not find quality bits for %s', trec{f});
  end;
  q1 = jsoc_cell2struct_keys(a.keywords);
  qstr = q1.quality{1}; % the string
  qchr = qstr(3:end); % mask off the leading 0x
  qnum = hex2dec(qchr); % actually a double

  % insert the numeric quality (a 32-bit integer stored as a double)
  q(f) = qnum;
  % insert the "ok for most purposes" flag
  ok(f) = (qnum < 2^31);
end;

return;
end
