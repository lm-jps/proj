function [bb,roimap,roi_ok,im_crit]=hmi_find_roi(s,geom,tp)
%hmi_find_roi	find sunspot groups in HMI labeling
% 
% [bb,roimap,roi_ok,im_crit]=hmi_find_roi(s,geom,tp)
% * Process labeling s to identify active region groups, returning them
% in a bounding-box set bb and an image-size map roimap, plus a 0/1 
% image-size indicator of below-threshold activity that is allowable.
% * The free parameters used in the identification are passed in 
% via the tp structure.
% * This API is used by tracker2_loop for sunspot-group identification
% components.
% * This function uses the newer hmi_patch omnibus mex-file to do
% the identification.
% * For debugging and diagnostic purposes, the pre-thresholded im_crit
% map can also be passed out.  This is not used by the tracker_loop
% API.
%
% Inputs:
%  real s(m,n);
%  real geom(5);  -- [c_x c_y r_sun b0 p0]
%  struct tp
% 
% Outputs:
%  real bb(nbox,4);
%  int  roimap(m,n);
%  bool roi_ok(m,n)
%  real im_crit(m,n);
% 
% See Also:  tracker2_loop, sunspotgroupfind

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 9 sep 2010
% Copyright (c) 2010.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [3]), error ('Bad input arg number'); end;
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end;
if ~isstruct(tp) || ~isfield(tp, 'tau'),
  error('Did not get a valid tp structure');
end;

%
% Computation
% 

%% find kernel
% (note, these definitions agree with hmi_patch)
nker = 256; % # entries in kernel lookup table
hi_ker = 0.015; % in units of distance on unit sphere
v1inx = linspace(0, hi_ker, nker); % just a dummy index
% the ker points are found using t2inx
ker = exp(-0.5*v1inx/((.0325*tp.kwid)^2));
% normalize it
ker = ker/(pi*(hi_ker/length(ker))*sum(ker));

% kernel weighting function
kwt  = [tp.klat 1 1];

% empty m-gram argument -> statistics not computed
[bb,stats,roimap,im_crit] = hmi_patch(s, [], geom, ...
                                      tp.active, ker, kwt, tp.tau);
% (due to empty m-gram argument, stats is junk)

% from crit, get a mask for sites that are allowed to be part of an ROI
% this is more permissive than roimap, which is (im_crit > tp.tau).
roi_ok = (im_crit > (tp.tau * tp.tau2));

return;
