function roisp=hmi_rois_truncate(rois,retain)
%hmi_rois_truncate	curtail history contained in ROI_s structures
% 
% roisp=hmi_rois_truncate(rois,retain)
% * Truncate the struct arrays contained in the cells of rois
% to `retain' entries.  Output the truncated cell array as roisp.
% * You can let retain = inf, in which case no truncation
% is done and roisp = rois.  You can let retain = 0, but I'm 
% not sure if an empty scroll works with the current tracker code.
% * Typical usage is to truncate past history of all chips
% to just the last entry:
% 
% >> roisp = hmi_rois_truncate(ROI_s, 1)
% >> save ROI_t ROI_rgn roisp
%
% The second line saves the truncated history as a starting point
% for next time.
%
% Inputs:
%  cell rois{nr} of struct
%  int retain;   -- or inf
% 
% Outputs:
%  cell roisp{nr} of struct
% 
% See Also:  

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 9 jul 2011
% Copyright (c) 2011.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [2]), error ('Bad input arg number'); end;
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end;
if ~iscell(rois) || ~isstruct(rois{1}),
  error('Did not get a valid rois cell array');
end;
if retain < 0 || floor(retain) ~= retain,
  % this is ok for inf
  error('Cannot retain a negative or fractional number of chips')
end;

%
% Computation
% 
if isinf(retain),
  % easist to special-case this one
  roisp = rois;
elseif false,
  % complete truncation
  roisp = cell(size(rois)); % fills roisp with []
  for r = 1:length(rois),
    % if rois{r} is full, it is a 1xnc struct array 
    % if rois{r} is empty, it is a 0x0 struct array with no
    %   fields, which still works below.
    % either way, the max() works with large "retain", even inf.
    roisp{r} = rois{r}(:,max(1,end-retain+1):end);
  end;
else,
  % just sparsify "chip" field within the struct
  roisp = rois;
    for r = 1:length(roisp),
      for c = 1:max(0,length(roisp{r})-retain),
      % sparse zeros, retains only size information
      C1 = roisp{r}(c).chip;
      roisp{r}(c).chip = sparse(size(C1, 1), size(C1, 2));
    end;
  end
end

return

