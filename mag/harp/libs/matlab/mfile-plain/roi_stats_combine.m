function sc=roi_stats_combine(s)
%roi_stats_combine	combine roi_stats_mag vectors into one
% 
% sc=roi_stats_combine(s)
% * Given a list of roi_stats vectors containing nb rows, one 
% for each ROI, and ns statistics, as defined in roi_stats_mag_defs.h, 
% construct a single vector with ns statistics that corresponds to
% the union of all the vectors.
% * This is done by averaging, summing, taking the min/max, or recomputing
% based on other statistics, as appropriate for that statistic.  The 
% correct operation is summarized by roi_stats_mag.
% * These combinations are essentially exact, except for daysgone and
% daysback.  These two are inexact because the velocity used to compute
% return times depends on the AR centroid.  The deviation is small,
% these parameters are inexact anyway, and we did not want to introduce
% the dependencies needed to make the computation exact.
% 
% Inputs:
%   real s(nb,ns)
% 
% Outputs:
%   struct sc(1,ns)
% 
% See Also: roi_stats_mag

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 22 Feb 2011
% Copyright (c) 2011.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [1]), error ('Bad input arg number'); end;
if all(nargout ~= [0 1]), error ('Bad output arg number'); end;

%
% Computation
% 
[nb,ns] = size(s);
% get the field names in order of index
[junk,names,combo] = roi_stats_mag([],[],[],zeros(1,5),0,'sesw');
if ns ~= size(names,1),
  error('Incommensurate #statistics between input s and roi_stats_mag');
end;
% convert to cells to make life easier
names = cellstr(names);
combo = cellstr(combo);

% the combined statistic
sc = zeros(1,ns);

% set divide-by-zero warnings to off, save state
swarn = warning('off', 'Matlab:divideByZero');

% fill in sc
for i = 1:ns,
  % determine operation
  name = names{i};
  opname = combo{i};
  % perform operation
  if any(strcmp(opname, {'min', 'max', 'sum'})),
    op = str2func(opname);
    sc(i) = op(s(:,i));
  elseif ~isempty(strfind(opname, 'recomp')),
    % retrieve the (non-central) moments and the normalizer
    %   below, let N = rgn_area ("normalizer")
    % all outputs are scalars
    [S_1,S_2,S_3,S_4,N] = stat_retrieve(sc, names);
    % let mu = first moment (extracted from sc)
    mu = S_1/N; % this is used throughout
    % let M_2 = variance = second central moment
    M_2 = (S_2/N) - mu.^2;
    % recomputation is case-by-case
    if strcmp(name, 'rgn_bmean'),
      % mean = rgn_bsum1 / N = mu
      sc(i) = mu;
    elseif strcmp(name, 'rgn_bsdev'),
      % sdev = sqrt( (rgn_bsum2 / N) - rgn_bmean^2 )
      sc(i) = sqrt(M_2);
    elseif strcmp(name, 'rgn_bskew'),
      % skewness = M_3 / (M_2 ^ 3/2) where algebra reveals
      % M_3 = rgn_bsum3/N - 3 rgn_bmean (rgn_bsum2/N) + 2 rgn_bmean^3
      M_3 = (S_3/N) - 3*mu*(S_2/N) + 2*mu^3; % third central moment
      sc(i) = M_3 / (M_2 * sqrt(M_2));
    elseif strcmp(name, 'rgn_bkurt'),
      % kurtosis = M_4 / M_2^2 - 3 where algebra reveals:
      % M_4 = (rgn_bsum4/N) - 4 rgn_bmean (rgn_bsum3/N) + ...
      %       6 rgn_bmean^2 (rgn_bsum2/N) - 3 rgn_bmean^4
      M_4 = (S_4/N) - 4 * mu * (S_3/N) + 6 * mu^2 * (S_2/N) - 3 * mu^4;
      sc(i) = (M_4 / (M_2*M_2)) - 3;
    end;
  elseif ~isempty(strfind(opname, 'avg.')),
    wtname = opname(5:end);
    wtinx = find(strcmp(wtname, names));
    if length(wtinx) ~= 1,
      warning(swarn); % restore warning state
      error('No unique weight %s for field %s (op-type %s)', ...
            wtname, name, opname);
    end;
    wts = s(:,wtinx);
    if any(wts < 0),
      warning(swarn); % restore warning state
      error('Negative weight %s for field %s (op-type %s)', ...
            wtname, name, opname);
    end;
    if sum(wts) > 0,
      wts = wts/sum(wts);
    else
      % this can happen legally, e.g., when every flux-weighted 
      % center, across all merged regions, had zero weight.
      % the following suppresses the /0 warning, but still
      % sets the resulting average to NaN
      wts = wts + NaN;
    end;
    % When averaging, need to make zero-weight * NaN contribute
    % zero, because the NaN arises from a 0/0 situation anyway.
    % Basically 0 * (0/0) should go to 0, and if we did not do
    % this, it will go to NaN instead, forcing the whole average 
    % to NaN.
    s1 = s(:,i);
    s1(wts == 0) = 0;
    sc(i) = dot(wts, s1);
  else,
    warning(swarn); % restore warning state
    error('Did not recognize op %s for field %s', opname, names{i});
  end;
end;
warning(swarn); % restore warning state
return;
end

% stat_retrieve: helper function to extract some fields for
% moments from the structure

function [m1,m2,m3,m4,normizer] = stat_retrieve(s, names)

inx = find(strcmp('rgn_bsum1', names));
m1 = s(:,inx);
inx = find(strcmp('rgn_bsum2', names));
m2 = s(:,inx);
inx = find(strcmp('rgn_bsum3', names));
m3 = s(:,inx);
inx = find(strcmp('rgn_bsum4', names));
m4 = s(:,inx);
% FIXME: assuming no area normalization
inx = find(strcmp('rgn_num', names));
normizer = s(:,inx);
return;
end

