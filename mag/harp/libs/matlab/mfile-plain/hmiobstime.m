function t=hmiobstime(trec)
%hmiobstime	find observation time from HMI T_REC index
% 
% t=hmiobstime(trec)
% * Given a T_REC index, return the observation time as "Matlab"
% dates, a la datenum.
% * Given a list of indexes, returns the time for each one.
% The list of strings can be a cell array or a string array, with the
% strings "stacked vertically" in each case; see cellstr for example.
% If a cell array is used, it must be nf-by-1.
% 
% Inputs:
%   string trec(nf);  -- valid T_REC indexes
% 
% Outputs:
%   real t(nf);
% 
% See Also:

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 31 Aug 2010
% Copyright (c) 2010.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [1]), error ('Bad input arg number'); end;
if all(nargout ~= [0 1]), error ('Bad output arg number'); end;

%
% Computation
% 

% convert input to cell array
if ischar(trec),
  trec = cellstr(trec);
end;

% initialize sizes
nf = size(trec,1);
t  = zeros(nf,1); 

% extract times
for f = 1:nf,
  % get the f'th one
  trec1 = trec{f};
  % use all but the _TAI on the end
  t(f) = datenum(trec1(:,1:end-4),'yyyy.mm.dd_HH:MM:SS');
end;
return;
