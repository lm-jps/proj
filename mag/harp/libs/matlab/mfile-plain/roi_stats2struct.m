function s2=roi_stats2struct(s1)
%roi_stats2struct	convert roi_stats matrix into a structure
% 
% s2=roi_stats2struct(s1)
% * Given a roi_stats matrix containing nb rows, one for each ROI,
% and ns statistics, as defined in roi_stats_mag_defs.h, convert
% it into a length-nb structure array.
% * So for instance, s1(3,12) might contain the minimum longitude of
% the third ROI; this would be referred to using s2 as s2(3).arminlon.
% A column vector of all such longitudes would be [s2.arminlon]'.
% * This calls roi_stats_mag to get the name-index correspondence.
% 
% Inputs:
%   real s1(nb,ns)
% 
% Outputs:
%   struct s2(nb)
% 
% See Also: roi_stats_mag

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
[nb,ns] = size(s1);
% get the field names in order of index
[junk,names] = roi_stats_mag([],[],[],zeros(1,5),0,'sesw');
if ns ~= size(names,1),
  error('Incommensurate #statistics between input s1 and roi_stats_mag');
end;

% num2cell converts to nb-by-ns cell array of individual numbers
% cellstr converts char array of names to a cell array of names 
%  (and strips spaces)
s2 = cell2struct(num2cell(s1), cellstr(names), 2);

return;
end
