function res=track_hmi_loop(fn, s_time, im, hooks)
%track_hmi_loop	identify sunspot groups in images and link in series
%
% res=track_hmi_loop(fn, s_time, im, hooks)
% * Identifies sunspot groups in a segmentation and processes 
% the list of ROI's by entering the ROI's in certain global variables.
% * jsoc_trec_loop should control this function.  The first three
% arguments are supplied automatically when invoked by jsoc_trec_loop;
% the last argument (hooks) is a catchall containing parameters and
% function hooks that is supplied to jsoc_trec_loop, and then passed
% downwards to this function.
% * Input fn is a T_REC value corresponding to a mask image. (fn is
% used to get image metadata.)  The matlab tile
% * Usage: Call the image-map function
%  jsoc_trec_loop('hmi.Marmask_720s[2010.07.01/1d]', [-1 1 1],...
%                 'track_hmi_loop', hooks);
% Using skip=[-1 1 1] here indicates we want the masks loaded
% fully, don't want to skip any masks, and want begin/end hook
% calls made to this function by jsoc_trec_loop.
% * If skip(1) is negative, the image-loader will use http to get the
% images; if positive, it will look in the filesystem.  See 
% jsoc_trec_loop and jsoc_imload.
% * `hooks' is a catchall length-1 struct containing these fields:
%   - first_track: 0 -> restart from saved state; >0 -> TID of first track
%   - retain_history: inf -> retain all; integer -> retain n chips
%   - filename: filename template info (struct)
%   - load_meta: loads a mask's geometry metadata
%   - load_quality: loads a mask's quality metadata
%   - save_track: saves a track (list-of-chips) to disk, DB, etc.
%   - loop_hook: script called near end of track match process
%   - find_roi: given segmentation, finds ROIs
%   - tracker_params: see find_roi for the parameters needed
% * If res is returned as empty, this indicates the file was skipped.
% 
% Inputs:
%   string fn;     -- tag naming a mask
%   real s_time;   -- matlab time number
%   real im(m,n);  -- the contents of fn
%   struct hooks(1);
% 
% Outputs:
%   struct res; -- contains several fields
% 
% See Also:  

% Written by Michael Turmon (turmon@jpl.nasa.gov) in 2010 and 2011.
% Copyright (c) 2010.  All rights reserved.  

% these global variables hold the scrolls of identified spot-groups, 
% and the list of regions in the scrolls
global ROI_t ROI_s ROI_rgn

% Naming conventions for track numbering schemes:
% "track id" - number from 1 up, uniquely identifying each and every
%    track that was found in the entire image sequence.
% "track slot" - number from 1 to length(ROI_t), identifying which
%    entry in ROI_t is storing a particular track's state information.
%    track slots are recycled once a track is written out.
% "track index" - index into a list of track slots.  Often goes
%    from 1..Nt where Nt is the number of occupied track slots,
%    but can index into shorter lists as well.  It is important to 
%    distinguish these numbers from track slot numbers.

% 
% Error checking
% 

if all(nargin  ~= [4]), error ('Bad input arg number'); end;
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end;
if isnumeric(fn),
  res = 0; % set it now, for consistency (non-empty)
  % grab the hooks for initialization and clean-up
  if fn == -Inf, tracker_loop_begin(hooks); return; end;
  if fn == +Inf, tracker_loop_end(hooks);   return; end;
  % it's not a hook we know
  error('Bad input args: fn is numeric');
end;
% allow im to be a cell of {mask, mag}
if iscell(im)
  mag = im{2};
  im = im{1};
else,
  mag = []; % this disables some metadata in roi_stats
end;

%-----------------------------------------------------------------
%
% Computation
% 
% Check that we have not seen this image before
if ROI_rgn.currtime == s_time,
  % overlap by one -- assume a fencepost issue
  res = sprintf('Image at %s just before last time %s. Skipping.\n', ...
                fn, datestr(ROI_rgn.currtime));
  % skip it
  return;
elseif ROI_rgn.currtime > s_time,
  % (alternatively, could error out here)
  res = sprintf('Image at %s before last time %s. Skipping.\n', ...
                fn, datestr(ROI_rgn.currtime));
  % skip it
  return;
end;

% image is new: update current progress info
ROI_rgn.currfile = fn;
ROI_rgn.currtime = s_time;
if length(ROI_rgn.firstfile) == 0,
  ROI_rgn.firstfile = fn;
end;

% check the QUALITY at the current time
[q,q_ok,msg] = hooks.load_quality(fn);
if ~isempty(msg),
  res = msg; % will continue
  return;
elseif ~q_ok,
  res = sprintf('Bad quality 0x%08x.  Skipping the image at %s', q, fn);
  return;
end;
clear q q_ok msg

% get geometry 5-tuple, msg is error message or ''
[s_geom,msg] = hooks.load_meta(fn);
if ~isempty(msg),
  res = msg; % will continue
  return;
end;
clear msg

% image will be processed by this tracker: update its number
ROI_rgn.imagenum = ROI_rgn.imagenum + 1; % have done one more image

%% tracker parameters
% map parameters
tkp_maptype = 'uint8';  % ROI_t.map is stored as this type
tkp_dbl2map = str2func(tkp_maptype); % type converter function
tkp_mapmax  = double(intmax(tkp_maptype)); % map is in [0..mapmax]
% wt = weight for old map (= 1 - "weight for current map")
% formula: wt = exp(-maprate * dt [days])
% e.g., 120 log(2) gives wt = 1/2 for dt = 720sec
tkp_maprate = hooks.tracker_params.maprate;
% track finalization parameters
tkp_final_num = ...
    hooks.tracker_params.final_num; % delay in frames before declaring final (<0)
tkp_final_time = ...
    hooks.tracker_params.final_time;% delay in days before declaring final
% pixels to pad placeholder ROIs for unseen tracks
tkp_pixpad = 4;

% get observation time for catfiles
s_timeS = [datestr(s_time,10) '/' datestr(s_time,6) ' ' datestr(s_time,13) ];
% demote the active tracks by one image
tlist = find(([ROI_t(:).tid] >= 0) & ([ROI_t(:).state] <= 0));
for t = tlist, ROI_t(t).state = ROI_t(t).state-1; end;

% Identify ROIs in current labeling
% b is ROI bounding boxes; roimap maps pixels to regions (1..Nr or 0)
[b,roimap,roi_ok] = hooks.find_roi(im, s_geom, hooks.tracker_params);
% set up image number of this segmenatation as a string
s_imagenumS = sprintf(hooks.filename.imagenum, ROI_rgn.imagenum);
% 
% Find correspondence matrix
% 
% initialize correspondence matrix
Nt0 = ROI_rgn.Nt0; % number of active tracks
Nr = size(b,1);  % number of new ROIs
clear b; % we find this again later, don't try to keep track of it now
rho = zeros(Nt0, Nr); % correspondence matrix (tracks to roi's)
% list of active track slots within ROI_t
tlist = find(([ROI_t(:).tid] >= 0) & ([ROI_t(:).state] <= 0)); 
% loop over tracks, determining their correspondence to new ROIs
Nr_orig = Nr; % just the range of the original ROIs, excluding grandfathers
tinx = 1; % silly, maintain the iteration number (track index)
t_harder = []; % remember which track slots we tried harder to continue
clear map_save; % placeholder for old maps: avoid repeatedly rotating them forward
[map_save(1:length(ROI_t)).map] = deal([]); % map_save(:).map = [] => uninitialized
for t = tlist,
  % note: t is a track slot, tinx is a track index
  % it is much faster to find bb1 here than in the mexfile - R2010b, mar 2011
  % (and it takes very little time here, so bb1 is not worth keeping track of)
  [i,j] = find(ROI_t(t).map); bb1 = [min(i) min(j) max(i) max(j)];
  % extrapolate roi(t).map, at prior time, to roi1 at current time
  % roi(t).map is a uint8.  *24 converts days -> hours.
  % roi1 lives between 0 and tkp_mapmax
  roi1 = rotate(double(ROI_t(t).map), ...
                ROI_t(t).mapgeom, s_geom, ...
                (ROI_t(t).maptime-s_time)*24, 'sparse,sesw', bb1);
  clear i j bb1;
  % save old map in case there's no match for this track
  % (only use the memory if the track was not matched last time)
  if ROI_t(t).state < -1,
    map_save(t).map     = ROI_t(t).map;
    map_save(t).mapgeom = ROI_t(t).mapgeom;
    map_save(t).maptime = ROI_t(t).maptime;
  end;
  % save map at (s_time, s_geom) now, for re-use later
  % (note, ROI_t(t).time != ROI_t(t).maptime)
  ROI_t(t).map = tkp_dbl2map(roi1); % converts nans to 0
  ROI_t(t).mapgeom = s_geom;
  ROI_t(t).maptime = s_time;

  % correlate roi1 and current roi map (only with original ROIs)
  % 1:Nr covers any new `try harder' ROIs also
  rho(tinx, 1:Nr) = cstat(roi1, roimap, Nr, 'unscaled');
  % try harder if a current track is unmatched
  if sum(rho(tinx,:)) == 0,
    mask = (roi1 > 0) & roi_ok;
    if nnz(mask) > 0,
      t_harder(end+1) = t; % list of `try harder' tracks
      % there was something present - call it a new ROI
      Nr = Nr + 1; % the new number of ROIs
      roimap(mask) = Nr; % plug in the current ROI
      rho(:,end+1) = 0; % append ROI to correspondence matrix
      % rho(tinx, Nr) = sum(znan(roi1(mask))); 
      % update rho(1:tinx,Nr)
      % the maps of tracks tlist(1:tinx-1) may overlap the new ROI
      % note, rho(tinx+1:end,Nr) is done as a matter of course in later loops
      for tinx2 = 1:tinx,
        % roi2 is the map for tlist(tinx2), analogous to roi1 for tlist(tinx)
        % note, this map was already rotated to match the coordinates of mask
        roi2 = ROI_t(tlist(tinx2)).map;
        % add agreeing pixels, in 0..mapmax
        rho(tinx2,Nr) = sum(znan(double(roi2(mask))));
      end;
    end;
  end;
  tinx = tinx + 1;
end;
% rescale rho to count agreeing pixels, rather than map values
rho = rho / tkp_mapmax; 
clear roi1 tinx mask roi2 tinx2

%
% detect merging tracks
%
% tbuds is a map on tlist x tlist. 
% if i ~= j and tbuds(i,j) = 1, then tlist(i) will merge with tlist(j)
tbuds = expm(double((rho*rho')>0)) > 0; % track index of matching tracks
t_merged_away = []; % remember the merged-out track *indexes*
t_merged_into = []; % merged-into track *slots*
for t = 1:Nt0,
  % t is a track index, tlist(t) is a track slot for ROI_t
  tbuds1 = find(tbuds(t,:)); % tbuds1 are track indexes
  if length(tbuds1) <= 1, continue; end; 
  if tbuds1(1) ~= t, continue; end; % only merge the first time we hit the group
  % perform the merge
  fprintf('merge:%s', sprintf(' %d', tlist(tbuds1)));
  % t1: the *track slot* of the merged-into track, an element of tlist
  % t1_inx: the index of the merged-into track within tbuds1
  t1 = tracker_merge_tracks(hooks, tlist(tbuds1));
  t1_inx = find(t1 == tlist(tbuds1));
  % all of tbuds1, except the merged-into element
  t_merged_away = [t_merged_away tbuds1([1:end] ~= t1_inx)];
  t_merged_into = [t_merged_into t1];
  % add the overlap back into the row of rho we will retain
  rho(tbuds1(t1_inx),:) = sum(rho(tbuds1,:));
end;
% update local variables to reflect the merges
Nt0 = ROI_rgn.Nt0; % number of active tracks
% remove the merged rows of rho; these rows were summed above
rho(t_merged_away,:) = []; 
tlist = setdiff(tlist, tlist(t_merged_away));
clear t1 t1_inx tbuds1 % don't need this, but retain tbuds

%
% match tracks to rois
%
% note:  pi1 is of length Nt0, maps to 1..Nr,  0 for unmatched
%        pi2 is of length Nr,  maps to 1..Nt0, 0 for unmatched
% and, tlist is of length Nt0, maps to 1..length(ROI_t)
[pi1,pi2,score] = assignment(ceil(rho)); % fn takes only ints
% combine regions if they diverge from a common track
for r = 1:Nr,
  if pi2(r) > 0, continue; end; % rgn r already matched
  if isempty(rho(:,r)) || max(rho(:,r)) == 0, continue; end; % rgn can't match
  tinx = find(rho(:,r) == max(rho(:,r)));
  tinx = tinx(1); % ensure a unique match
  r2 = pi1(tinx); % the matched region
  if r2 == 0, disp('r2=0, this cannot happen'); continue; end;
  % merge r into r2
  roimap(roimap == r) = r2;
  pi2(r) = -1; % code to ignore this region
end;

% 
% create new tracks as necessary
% 
for r = 1:Nr,
  % check for new ROIs
  if pi2(r) == 0,
    % allocate the new track
    t = tracker_new_track(ROI_rgn.Ntot); % t=new track index, in 1..length(ROI_t)
    % record in pi2 and tlist for subsequent loops; key condition is:
    % t = track index satisfies: t = tlist(pi2(r))
    pi2(r) = length(tlist) + 1;
    tlist(pi2(r)) = t;
    % born at: current time + 1 sec per new rgn, for uniqueness of birth times
    ROI_t(t).birth = s_time + sum(pi2 == 0) / (24*3600);
    % this (time,map) value plays well with the data assimilation step below
    %   (all the values below will be reset in the loop below)
    ROI_t(t).time    = -Inf; % last time prev. seen
    ROI_t(t).map     =    0; % map prev. seen (like: zeros(size(im),tkp_maptype))
    ROI_t(t).maptime = -Inf; % time when map above was observed (unused)
    ROI_t(t).mapgeom = [0 0 0 0 0]; % map geometry (unused)
  end;
end;

% 
% data assimilation
%    - combine existing ROI_t.map with roimap to make updated map
%    - plus, simple bookkeeping updates to snapshot of ROI_t(t)
% 
% Also finds the all-track "posterior most probable" map -- the most likely 
% track at each image pixel.
roimap_p  = zeros(size(roimap)); % most likely track (p = posterior)
roiweight = zeros(size(roimap), tkp_maptype); % and its weight
for r = 1:Nr,
  if pi2(r) < 0,
    continue; % r was combined above: skip
  end;
  % existing track (possibly new)
  t = tlist(pi2(r)); % ROI r matches track t
  % compute the "data-assimilated" map (in 0..tkp_mapmax)
  %   ROI_t(t).map was transformed to current image geometry above
  %   it is usually an image-size integer, or scalar 0 for new track
  %   using ROI_t.time, not .maptime, ensures dt increases if track not seen
  dt = s_time - ROI_t(t).time; % in days, > 0, Inf for new track
  newmap = tkp_dbl2map(tracker_weight_map(ROI_t(t).map, roimap==r, dt, ...
                                          tkp_maprate, tkp_mapmax));
  % save map, set up .mapgeom and .maptime in case of new track
  ROI_t(t).map     = newmap;
  ROI_t(t).maptime = s_time;
  ROI_t(t).mapgeom = s_geom;
  % update the posterior-most-probable map
  % if the new map has higher weight, the region "r" wins
  in_rgn = (newmap > roiweight); % pixel is in most likely region
  roimap_p (in_rgn) = r; % plug in the region number
  roiweight(in_rgn) = newmap(in_rgn); % plug in the higher weights
  % update the track snapshot ROI_t(t) with the matched ROI
  % (.birth does not change, and we updated .map above)
  ROI_t(t).state = 0;      % r matched t: make t a first-class active track
  ROI_t(t).time  = s_time; % most recent time track was found
  ROI_t(t).last  = ROI_rgn.imagenum; % most recent image# matched
end;
clear dt in_rgn newmap roiweight;  % do not need, but can keep as a diagnostic

% save roimap and other side info
if ischar(hooks.loop_hook),
  eval(hooks.loop_hook);
end;

%
% gather roi statistics, one row per roi
% 
% remake the bb's for roimap_p
b = region_bb(roimap_p, [1 1]); % 1-based coordinates
% compute statistics using roimap_p
s = roi_stats_mag(roimap_p, ...
                  double(im==hooks.tracker_params.active), ...
                  mag, s_geom, -1, 'sesw');

% 
% matched regions: update ROI_t and ROI_s
% 
for r = 1:Nr,
  if pi2(r) < 0,
    continue; % r was combined above: skip
  end;
  % existing or new track
  t = tlist(pi2(r)); % ROI r matches track slot t
  % add region r to track t: update the segmentation scroll, ROI_s{t}
  tracker_cons_ROI(t, tracker_mask_map_encode(im, roimap_p, b(r,:), r), ...
                   b(r,1:2), s_time, s_geom, s(r,:), fn, r, ...
                   r > Nr_orig, any(t == t_merged_into));
end;

% 
% Track maintenance: adjust current state, and write and reset if ROI vanishes
%
% note, this loop does not cover any track just created above
for t = tlist,
  if ((s_time-ROI_t(t).time)>tkp_final_time) && (ROI_t(t).state<tkp_final_num),
    % AR has vanished: close it out
    tracker_close_track(t, hooks, 'final'); % finalized: write summary file
    continue;
  end;
  if ROI_t(t).state < 0,
    % AR was unmatched by an ROI in current frame
    roimap1 = (ROI_t(t).map > 0);
    [i,j] = find(roimap1); % indexes within roi1
    if isempty(i),
      % all pixels on AR have vanished: close it out
      tracker_close_track(t, hooks, 'final'); % finalized: write summary file
      continue; 
    end;
    % insert a placeholder to ROI_s{t}
    % roimap1 was extrapolated from last-seen map, it is a logical
    b1 = [max(min(i)-tkp_pixpad,1)          max(min(j)-tkp_pixpad,1) ...
          min(max(i)+tkp_pixpad,size(im,1)) min(max(j)+tkp_pixpad,size(im,2))];
    % mask_map_encode last argument is 1 because roimap1 == 1 is on-roi
    % note, the stats argument will never be output, but it may be merged-in
    stats1 = roi_stats_mag(0, 0, [], s_geom, 1, 'sesw');
    tracker_cons_ROI(t, tracker_mask_map_encode(im, roimap1, b1, 1), ...
                     b1(1:2), s_time, s_geom, stats1, fn, 0, false, false);
    % restore old map, so it's not repeatedly rotated
    if ~isempty(map_save(t).map),
      ROI_t(t).map     = map_save(t).map;
      ROI_t(t).mapgeom = map_save(t).mapgeom;
      ROI_t(t).maptime = map_save(t).maptime;
    end;
  end;
end;

% output something useful
tinx = [ROI_t(:).tid].' > 0;
Ntot = nnz(tinx); % total # tracks
res.nums = [Nr_orig Nr nnz(pi2==0) ROI_rgn.Nt0 ROI_rgn.Nt1]'; 
% Commented out: currently unused.
if 0,
  res.rho = rho;
  res.pi1 = pi1;
  res.pi2 = pi2;
  res.score = score;
  res.nacr = nnz(im == hooks.tracker_params.active); % #active pixels
  res.tracks = [ROI_t(tinx).tid].';
  res.tlist = tlist;
  res.t_merged_away = t_merged_away; % track slots: merged 
  res.t_harder = t_harder; % track slots: tried harder to match and succeeded
  % geom: Ntot by whatever
  % res.geom = [[ROI_t(tinx).time].' reshape([ROI_t(tinx).mapgeom], [], Ntot).'];
end;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialization: set up data structures, establish cat-file headers
% ROI_rgn: set most fields to empty or zero
% ROI_t, ROI_s: clear
% Start with some empty slots in the track array
%
function tracker_loop_begin(hooks)
global ROI_t ROI_s ROI_rgn % state information

% 1: ROI_rgn (data pertaining to all the tracks)
%    note - ROI_rgn is not null at outset - check it
if ~exist('ROI_rgn'),
  error('bad initialization: ROI_rgn must exist as a global');
end;
% Must either load existing state, or init to empty
if hooks.first_track == 0,
  % attempt to load existing state
  % TODO: abstract this into a hook, it's very hmi-specific
  fn = sprintf(hooks.filename.template, 'jsoc', '-prior', '', 'mat');
  if exist(fn, 'file'),
    % prior state file exists; load it
    prior = load(fn);
    ROI_rgn = prior.ROI_rgn;
    ROI_t   = prior.ROI_t;
    ROI_s   = hmi_rois_truncate(prior.ROI_s, hooks.retain_history);
    ROI_rgn.curr_run = ROI_rgn.curr_run + 1;
    ROI_rgn.firstfile = ''; % first file examined in this run
    fprintf('Restart: track ID = %d, pending tracks = %d, history_retain = %d\n', ...
            ROI_rgn.Ntot, ROI_rgn.Nt0, hooks.retain_history);
  else,
    % announce we could not restart
    error('Could not find prior-state file "%s", fatal.', fn);
  end;
else,
  fprintf('Fresh start at track %d (no prior-state file)\n', hooks.first_track);
  ROI_rgn.Ntot = hooks.first_track;     % index of next track we will create
  ROI_rgn.Nt0 = 0;      % #tracks currently active
  ROI_rgn.Nt1 = 0;      % #tracks all together
  ROI_rgn.filenum = 0;  % #files processed (unused)
  ROI_rgn.imagenum = 0; % #image frames processed
  ROI_rgn.currfile = '';% current filename is TBD
  ROI_rgn.currtime = 0; % current file's time is TBD
  ROI_rgn.curr_run = 1; % current run number
  ROI_rgn.firstfile = ''; % first file examined in this run
  % 2: ROI_t, ROI_s: tracks and scrolls are empty
  ROI_t = struct([]);
  ROI_s = {};
  % set up a few empty slots (not really necessary)
  for t = 1:32,
    tracker_new_track(-1); % -1 => unallocated track
  end;
end;
% use supplied hook to make headers
hooks.init_run(-Inf, hooks, ROI_rgn);
% Now ROI_rgn, ROI_s, ROI_t are OK
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Termination: write out partial results
%
function tracker_loop_end(hooks)
global ROI_t ROI_s ROI_rgn % state information

% save state in all cases
% (fast to write the file, allows later restart, do not overwrite prior state)
fn = sprintf(hooks.filename.template, 'jsoc', '-post', '', 'mat');
[junk1,junk2] = mkdir(fileparts(fn)); % ensure the dir exists
% truncate ROI_s
ROI_s = hmi_rois_truncate(ROI_s, hooks.retain_history);
% save as -v7.3 because ROI_s can be larger than 2GB after unpacking
% (v7.0 is the default)
save(fn, 'ROI_t', 'ROI_s', 'ROI_rgn', '-v7.3');
% loop over all allocated tracks
for t = find([ROI_t(:).tid] >= 0),
  tracker_close_track(t, hooks, 'pending'); % still pending
end;
% write file tails, if desired
hooks.init_run(Inf, hooks, ROI_rgn);
% no need to keep these around
clear ROI_t ROI_s ROI_rgn
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Find/allocate new, empty track
%
% input: new track id (int); output: new track index in heap (>=1)
% if new track id < 0, a completely empty track-mosaic is allocated,
% with an empty scroll.  it would then be filled in later by a
% subsequent call; this might be useful for initialization or for
% flushing of a track
function t = tracker_new_track(tid)
global ROI_t ROI_s ROI_rgn % state information

% 0: find a slot
if tid >= 0 & ~isempty(ROI_t),
  t = find([ROI_t(:).tid] < 0, 1); % take first one (if none, t = [])
else
  t = []; % only at program initialization, when ROI_t = []
end;
if isempty(t),
  % make a new slot (includes case when ROI_t = ROI_s = [])
  t = length(ROI_t) + 1; % one more track
end;
% 1: update global information (ROI_rgn)
if tid >= 0, 
  ROI_rgn.Nt0 = ROI_rgn.Nt0 + 1;  % new track will be active
  ROI_rgn.Nt1 = ROI_rgn.Nt1 + 1;  % any new track ups this
end;
if tid == ROI_rgn.Ntot, ROI_rgn.Ntot = tid+1; end; % up total #tracks
% 2: clear scroll part (ROI_s)
ROI_s{t} = struct([]);
% 3: clear track part (ROI_t)
tracker_clear_track(t);
ROI_t(t).tid = tid; % this must be set correctly; it indicates used/unused
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Close out track
%
% inputs:
%   t: track slot number in heap (>=1)
%   hooks: struct of hooks
%   status: 'merged' 'pending' 'final'
% this routine: 
%   computes and writes the scroll connected with track indexed t
%   updates the global information to indicate the demise of track t
%   clears the side information connected with the track (ROI_t(t))
function tracker_close_track(t, hooks, status)
global ROI_t ROI_s ROI_rgn % state information

% write its scroll (and metadata), if desired
if any(strcmp(status, {'pending', 'final'})),
  % use hook for file/db saving
  hooks.save_track(ROI_s{t}, ROI_t(t).tid, ROI_t(t).birth, ...
                   status, hooks.run_name, hooks.filename.template);
end;
ROI_s{t} = struct([]); % clear scroll part
% update global track information
if ROI_t(t).state <= 0,
  ROI_rgn.Nt0 = ROI_rgn.Nt0 - 1; % one fewer active track (Nt0)
end;
ROI_rgn.Nt1 = ROI_rgn.Nt1 - 1; % one fewer track (Nt1)
% clear track part
tracker_clear_track(t);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Wipe a track data structure clean
%
% input: track slot number in heap (>=1)
% clears all side information connected with the track, ROI_t(t)
function tracker_clear_track(t)
global ROI_t ROI_s ROI_rgn % state information

ROI_t(t).tid = -1;    % this must be set correctly; it indicates used/unused
% below here, we do not attempt to set up "correct" values
ROI_t(t).state = 0;   % just-seen or dying
ROI_t(t).map = [];    % current map indicator
ROI_t(t).maptime = 0; % matlab time corresponding to map
ROI_t(t).mapgeom = [0 0 0 0 0]; % geometry for map
ROI_t(t).birth = 0;   % matlab time of the first frame the track was found
ROI_t(t).time = 0;    % time when track was last seen
ROI_t(t).last = 0;    % image number when track was last seen
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Merge tracks
%
% input: 
%   hooks, side information for closing merged tracks
%   tlist, a nonempty list of track slot numbers (1..Nt) to merge
% output:
%   t1, the track slot which was merged-into, a member of tlist
%
% The merged track slots, tlist \ {t1}, are closed and written out
%
% Currently t1 is chosen as the track with the longest history
% among tlist
% 
function t1 = tracker_merge_tracks(hooks, tlist)
global ROI_t ROI_s % state information

% 1: Find the track number to preserve
% t1_inx is index of the longest scroll among tlist
[junk,t1_inx] = max(cellfun('length', ROI_s(tlist)));
t1 = tlist(t1_inx); % will merge (tlist \ {t1}) into t1
% 2: Merge the scrolls, ROI_s(tlist)
rois = [ROI_s{tlist}]; % all rois in all tracks
sp1 = rois(1); % prototype to ensure the field names agree
sp1.chip = sp1.chip(1); % economize the prototype, preserving type
[flist,i,j] = unique([rois.frame]); % all unique image frame numbers
nf = length(flist); % number of frames in sp
% merge the scrolls into a new scroll `sp(f)', one frame at a time
for f = 1:nf,
  sp(f) = sp1; % set the names up correctly upon creation
  % set up the fields of sp(f)
  sp(f).frame = flist(f);
  rois1 = rois(j == f); % just the rois for this frame
  nr1 = length(rois1); % how many of these
  sp(f).time = rois1(1).time; % all times in these roi's will be the same
  sp(f).geom = rois1(1).geom; % all geoms are the same
  sp(f).fn   = rois1(1).fn;   % all filenames are the same
  sp(f).tag  = max([rois1.tag]); % 2 beats 1, 1 beats 0, 0 beats 1
  sp(f).rgn  = sort(cat(2, rois1.rgn)); % cat region IDs together
  % merge statistics
  sp(f).stats = roi_stats_combine(reshape([rois1.stats], ...
                                          size(rois1(1).stats, 2), [])');
  % merge coords
  coords1 = min(cat(1, rois1.coords), [], 1); % min(array of all corners) ...
  sp(f).coords = coords1; % ... = "top-left" corner = new corner
  % loop to find size of roi containing each rois1(:).chip
  % (do so algorithmically by finding its "bottom-right" corner coords2)
  coords2 = zeros(1,2);
  for r = [1:nr1], 
    coords2 = max(coords2, rois1(r).coords + size(rois1(r).chip) - 1);
  end;
  % loop to plug all the chips in
  chip1 = zeros(coords2 - coords1 + 1, 'uint8'); % empty field
  for r = [1:nr1], 
    off1 = rois1(r).coords - coords1; % offset: old origin - new origin (>=0)
    box1 = size(rois1(r).chip); % just the chip size
    % define the box to update
    inx1 = [off1(1)+1:off1(1)+box1(1)];
    inx2 = [off1(2)+1:off1(2)+box1(2)];
    % FIXME: need smarter update to not overwrite good info with zeros
    chip_now = chip1(inx1,inx2);
    % chip1(inx1,inx2) = max(chip_now, full(rois1(r).chip)); % (new)
    chip1(inx1,inx2) = full(rois1(r).chip); % (original - bug)
  end;
  sp(f).chip = chip1; % plug it in
end;
% (merged-to status is indicated in .tag of present frame, not yet in ROI_s)
ROI_s{t1} = sp; % plug the new entry in
% 3: Merge the track state, ROI_t(tlist), into ROI_t(t1)
ROI_t(t1).state = max([ROI_t(tlist).state]); % highest state
ROI_t(t1).time  = max([ROI_t(tlist).time]);  % latest time seen
ROI_t(t1).birth = min([ROI_t(tlist).birth]); % earliest appearance
% do not alter: geom, last, maptime
% get rid of redundant tracks
for t = tlist,
  if t == t1, continue; end; % do not get rid of the merged-into track
  % maps are already commensurate: same coord system (.mapgeom) and .maptime
  ROI_t(t1).map = ROI_t(t1).map + ROI_t(t).map; % combine uint8's
  tracker_close_track(t, hooks, 'merged'); % note, this dumps ROI_s{t} also!
end;
% (returns t1)
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% add the ROI chip to track slot t, with useful side info
%
function tracker_cons_ROI(t, chip, coords, s_time, s_geom, stats, fn, r, harder, merge)
global ROI_t ROI_s ROI_rgn % state information

nr = length(ROI_s{t}) + 1; % the new region number
ROI_s{t}(nr).chip   = uint8(chip);
ROI_s{t}(nr).coords = coords; % (i,j) coords, not (x,y)
ROI_s{t}(nr).frame  = ROI_rgn.imagenum;
ROI_s{t}(nr).time   = s_time; % numeric time in matlab format
ROI_s{t}(nr).geom   = s_geom;
ROI_s{t}(nr).stats  = stats;
ROI_s{t}(nr).fn     = fn;     % current filename as a string
% set up roi.rgn, a list.  note, .rgn's are concatenated after merge!
% (do not use 0 for placeholder because 0 is not a valid region number)
if r > 0,
  ROI_s{t}(nr).rgn = r;    % region number (after merges, can be a list)
else
  ROI_s{t}(nr).rgn = [];   % is a placeholder region
end;  
% set up roi.tag, an integer.
if ROI_t(t).state == 0,
  if merge,
    ROI_s{t}(nr).tag = 2;  % regular ROI, merged-into on this frame
  elseif harder,
    ROI_s{t}(nr).tag = 0;  % regular ROI, had to try harder
  else,
    ROI_s{t}(nr).tag = 1;  % regular ROI
  end;
else,
  ROI_s{t}(nr).tag   = -1; % just a placeholder ROI
end;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Encode both the mask (0/1/2 for offdisk, quiet, active) and the
% roi map into one 8-bit number.  This follows the convention for
% HMI patches.  The encoded value once was:
%  code = mask + 32 * [within_other_roi] + 64 * [within_this_roi]
% but is now:
%  code = mask + 16 * [within_other_roi] + 32 * [within_this_roi]
% To extract just the present patch, use: (code > 32).
% Arguments: mask (region type), map (roi map), b (bounding box),
% region number to compute within_this_roi
function coded = tracker_mask_map_encode(mask, map, b, r)

mask1 = mask(b(1):b(3), b(2):b(4));
map1  = map (b(1):b(3), b(2):b(4));
within_roi = (map1 == r);
coded  = mask1 + 16*((map1 > 0) & ~within_roi) + 32*within_roi;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% tracker_weight_map: combine an old map (of integer type) 
% and a current map (0/1), at a time difference dt (days) apart.
% They are already in the same image geometry (center, scale, angle)
% A standard exponential formula based on a rate (in inverse days)
% is used to compute weights, summing to one, to linearly combine 
% the two maps.
%
% NOTE: Currently, the weights do not sum to one.  The 0/1 ROI,
% instead, saturates the result to mapmax, if it is 1.
%
% The result is returned as a `single' type in 0..mapmax, which can be
% converted to the map type (e.g., uint8) by the caller
function newmap = tracker_weight_map(oldmap, roi_0_1, dt, maprate, mapmax)

% lower the weight on the prior map as dt rises
% (dt = Inf is OK, and results in zero weight on the old map)
oldwt = exp(-dt * maprate);
% this expression is of the form:
%   floor(scalar_double * single + scalar_double * single)
% first, it evaluates faster in single than in uint8
% second, using singles allows us to round down (floor), preventing errors
% in which fading oldmap entries of small uint8 values cannot ever round 
% down to zero: e.g., newmap = uint8(0.4 * 0 + 0.6 * 1) = uint8(0.6) = 1
newmap = floor(mapmax*single(roi_0_1) + oldwt*single(oldmap));
return

