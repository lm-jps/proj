%track_hmi_production	script to find tracks in JSOC data series
%
% * This script is the top-level Matlab driver for HMI production.
% If calling from outside Matlab, you can use a companion Bourne
% shell script that invokes this routine.
% * All four of these variables must be defined:
%   mask_series:  
%       the JSOC data series for the mask images
%   dest_dir:
%       the destination directory to preserve the tracking
%       system state
%   first_track:
%       if 0, indicates we are restarting the tracker where we
%       left off, using the state in dest_dir.
%       if greater than 0, gives the number of the first track.
%       This implies DELETING THE OLD STATE, so it is to be used
%       with caution.
%   retain_history:
%       if inf, indicates we retain all history in saved checkpoint
%       files.  if a positive integer, only this many past frames
%       of each track is retained.
%       it truncates history both in loading and saving checkpoints.
% The following variables are optional:
%   mag_series: defaults to `hmi.M_720s'
%       for NRT processing, set to: `hmi.M_720s_nrt'
%   make_movie: defaults to false
%       if true or nonzero, compiles a large .avi movie of the
%       time history of the tracker; it's wise to transcode this
%       movie into a more compact file format
%   setup: defaults to empty
%       a function taking and returning three arguments, allowing modification
%       of tracker_params, hooks, and movie_params, in that order, for 
%       instrument or data product customization.
%   params: defaults to empty
%       contains a semicolon-separated list of modified tracker parameters.
%       if both setup and params are given, params is done after setup.
%   gap_fill: integer, defaults to nan (= no gap filling)
%       if non-nan, new track ID numbers will be chosen in agreement
%       with future values, or recycled from previously-unused numbers.
%       The gap_fill integer, typically zero, is added to ROI_rgn.Ntot,
%       the track ID counter, at the end of the run, to provide another
%       way to synch the track IDs of the gap-fill series with existing 
%       IDs.
%  
% Inputs:
%   string mask_series
%   string dest_dir
%   int first_track
%   int retain_history (or +inf)
%
% Outputs:
%   (to disk and to global variables)
% 

% tracker state
global ROI_s ROI_t ROI_rgn

fprintf('\n');
fprintf('%s: Path info: now running: %s\n', mfilename, mfilename('fullpath'));
fprintf('%s: Path info: track_hmi_loop is at: %s\n', mfilename, which('track_hmi_loop'));
fprintf('\n');

% require mask_series
if ~exist('mask_series', 'var'),
  error('Require "mask_series" to be defined');
end;
% require dest_dir
if ~exist('dest_dir', 'var') || length(dest_dir) == 0,
  error('Require "dest_dir" to be defined and nonempty');
end;
% require first_track
if ~exist('first_track', 'var'),
  error('Require "first_track" to be defined (= 0 to restart)');
end;
% require retain_history
if ~exist('retain_history', 'var'),
  error('Require "retain_history" to be defined (= inf to retain all)');
end;

% magnetogram data series: optional, but default to M_720s
if ~exist('mag_series', 'var'),
  mag_series = 'hmi.M_720s{magnetogram}';
end;
% make_movie: optional, default to false
if ~exist('make_movie', 'var'),
  make_movie = false;
end;
% gap_fill: optional, default to false
if ~exist('gap_fill', 'var'),
  gap_fill = nan;
end;
% setup function: optional, default to none
if exist('setup', 'var') && ~isempty(setup),
  if ~exist(setup, 'file'), error('Supplied setup function (%s) not present', setup); end;
  param_setup = str2func(setup); % turn into a function handle
else,
  setup = '';
  param_setup = @(tp1,h1,mp1)deal(tp1,h1,mp1); % dummy that copies three args thru
end;
% params: optional, default to none
if ~exist('params', 'var'),
  params = '';
end;

% set up NRT mode if mag_series contains _nrt
%   (used internally for WCS and QUALITY)
if ~isempty(strfind(mag_series, '_nrt')),
    hmi_property('set', 'nrt_mode', 1);
end;

% Echo what we're about to do
if first_track == 0,
  first_run_string = 'Appending to ';
else,
  first_run_string = 'Removing';
end;
fprintf('%s: Tracking from "%s" and "%s".\n', mfilename, mask_series, mag_series);
fprintf('%s: %s existing files in "%s".\n', mfilename, first_run_string, dest_dir);
clear first_run_string

% delete existing files if the tracker is being restarted
if first_track > 0,
  system(sprintf('rm -rf %s/Tracks', dest_dir));
end;

%% establish a text tag for this run
% (tag is a metacharacter-free encoding of the run name)
run_tag = mask_series;
% re-encode the tag to remove any weird chars (/ is the main issue)
% RE is of form: [^xyz] where xyz are good
% note - has to be the last char because it's a [] metacharacter
% note - we use system() on this string, so it has to have NO metacharacters
not_goodchar = '[^\w@+=%.,-]';
run_tag = regexprep(run_tag, not_goodchar, '_');
% kill multiple _'s
run_tag = regexprep(run_tag, '_+', '_');
% kill _ at end, if any
run_tag = regexprep(run_tag, '_$', '');


%
% Tracker setup
%

% tracker internal parameters
% turmon dec 2010 -- these are rescaled from the standard MDI 
% choices, and are the same as those used in various HMI test runs.
% turmon mar 2011 -- added in maprate and final* parameters
% changed tau2 from 0.5 to 0.8
tracker_params1.kwid = 1;
tracker_params1.klat = 4;     % stretch factor along latitude
tracker_params1.tau  = 0.018; % threshold
tracker_params1.tau2 = 0.8;   % tau reduction factor for "try harder"
tracker_params1.active = 2;   % mask code for "activity"
tracker_params1.mode = 0;     % plain (0), or laplacian (1) [unsupported]
tracker_params1.final_num  = -2;  % delay in frames before declaring final (<0)
tracker_params1.final_time = 2; % delay in days before declaring final
tracker_params1.maprate = 15*log(2); % gives wt, weight for old map
                                     % wt = exp(-maprate * dt [days])
                                     % 30 log(2): wt = 1/2 @ dt = 4*720sec

% handlers for utilities, abstracted from the tracker itself
hooks1.load_meta      = @(t)(hmidisk(t, 'geom,vector,quiet'));
hooks1.load_quality   = @hmiquality;
hooks1.save_track     = @hmi_save_track_production;
hooks1.init_run       = @hmi_init_run_production;
hooks1.loop_hook      = 'hmi_loop_hook_local'; % a script, not a function
hooks1.find_roi       = @hmi_find_roi; % uses tracker_params
hooks1.first_track    = first_track;
hooks1.retain_history = retain_history;
hooks1.run_name       = mask_series; % text tag, just for identification
hooks1.gap_fill       = gap_fill; % gap-fill mode, if non-nan (usually nan)
% filename templates
hooks1.filename.template = fullfile(dest_dir, '/Tracks/%s/track%s%s.%s');
hooks1.filename.tracknum = '%06d';
hooks1.filename.imagenum = '%06d';
% depending on gap_fill, set up hook to install new track IDs
if ~isnan(hooks1.gap_fill),
  hooks1.update_tid = @hmi_gapfill_update_tid;
else,
  hooks1.update_tid = @(varargin)(disp('Unexpected negative TID'));
end;

% movie-maker parameter, loaded with standard HMI arguments
% movie_params1.load_meta = hooks1.load_meta;
movie_params1.string_title = 'SDO/HMI Tracker Diagnostic';
movie_params1.string_arname = 'HARPs';
% track metadata file template
movie_params1.file_track_mat = sprintf(hooks1.filename.template, 'mask', '-mask.%s', '', 'mat');
% movie output file template
movie_params1.file_movie     = sprintf(hooks1.filename.template, 'movie', '-movie.', [run_tag '_%d'], 'avi');
movie_params1.frame_dims = 1024; % target movie dimens
movie_params1.frames_per_sec = 20;

% allow the params to be transformed
if ~isempty(setup), fprintf('%s: Running setup: %s\n', mfilename, setup); end;
[tracker_params,hooks,movie_params] = param_setup(tracker_params1, hooks1, movie_params1);

% allow modification of tracker_params through passed-in `params'
%   find statements like: key1=val1;key2=val2;key3=val3
paramsep = unique([0 strfind(params, ';') length(params)+1]);
for i = 1:length(paramsep)-1,
  param1 = params(paramsep(i)+1:paramsep(i+1)-1);
  if isempty(param1), continue; end;
  fprintf('%s: \tparam = <%s>\n', mfilename, param1);
  inx1 = strfind(param1, '=');
  if length(inx1) ~= 1,
    error('Did not find key=val in param %s from %s', param1, params);
  end;
  key = param1(1:inx1-1);
  val = param1(inx1+1:end);
  if isempty(val),
    error('Did not find key=val in param %s from %s', param1, params);
  end;
  if ~isfield(tracker_params, key),
    error('Did not find known key %s in param %s from params %s', key, param1, params);
  end;
  value = str2double(val);
  if isnan(value),
    % currently, no need to pass in NaNs, so can use the error code
    error('Could not convert "%s" in param "%s" to double', val, param1)
  end;
  % at long last
  tracker_params.(key) = value;
end;
clear paramsep i param1 inx1 key val

if length(params) > 0,
  fprintf('%s: Modified tracker_params according to params:\n', mfilename);
  disp(tracker_params);
end;

% NB: the tracker params are passed in via "hooks"
hooks.tracker_params = tracker_params;

%% (end of parameter setup)

jsoc_ds = {mask_series, mag_series};

% display tracker parameters for capture in the log file
fprintf('%s: Tracker parameters used:\n', mfilename);
disp(hooks.tracker_params)

% switches from loading images by http versus filesystem
% (see jsoc_trec_loop for more)
if ~isempty(regexpi(getenv('HOST'), 'stanford')),
  smpl = 1; % filesystem
else,
  smpl = -1; % http
end;

%% Run the tracker
fprintf('%s: Beginning tracking.\n', mfilename);
[fn,res,fnfail,fnlog]=jsoc_trec_loop(jsoc_ds, [smpl 1 1], 'track_hmi_loop', hooks);

% preserve log file and failures, if any
track_hmi_postrun_files(fullfile(dest_dir, 'Tracks', 'jsoc'), ...
                        hooks.run_name, fnlog, fnfail);

fprintf('%s: Completed tracking.\n', mfilename);

%% Make a movie

if make_movie && length(fn) > 0,
  % turn off annoying warnings about NOAA metadata
  warning('off','hmi_base:no_noaa_info');
  % interval between frames
  dt = 1; 
  % master cat file
  fnIC = sprintf(hooks.filename.template, '', '-fd-instant', '', 'cat');
  % loop which makes the movie
  fprintf('%s: Beginning movie.\n', mfilename);
  [fn,res] = file_loop_cat(fnIC, [0 dt 1], @track_hmi_movie_loop, movie_params);
  fprintf('%s: Completed movie.\n', mfilename);
end;

fprintf('%s: Done.\n', mfilename);
return
