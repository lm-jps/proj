%track_hmi_harp_movie_script -- script for harp quicklook frames or movie
%
% * This script is the top-level Matlab driver for making quick-look images or
% movies using harps already ingested into JSOC.  If calling from outside Matlab, 
% you should probably use the companion shell-script driver (see below).
% * Options are passed by variables in the workspace.
% * These variables must be defined:
%   mask_series:  
%       the JSOC data series for the masks, with trec range and {mask}
%   harp_series:  
%       the JSOC data series for the HARPs, typically hmi.Mharp_720s
%   dest_dir:
%       the destination directory where images are to be placed
%       (it will be created if it does not exist)
% * These variables are optional:
%   mode:
%       a comma-separated list -- default = "frame"
%       if it contains "frame", output frames as png's for each T_REC
%       if it contains "movie", output a movie for each T_REC
%   setup:
%       A function, default null, which can reset some default values.
%       This is used to implement per-instrument movie modes.
%   harp_selector:
%       A SQL qualifier string selecting a subset of HARPs at each T_REC, e.g.
%       'HARPNUM > 1000' to choose only these HARPs.  Will be surrounded with
%       [?...?] and used in the query finding HARPs to show.
% * The times of mask_series control which frames are produced, e.g.:
%   mask_series = 'hmi.Marmask_720s[2012.09.21_00:00_TAI/2d]{mask}';
%   mask_series = 'hmi.Marmask_720s[2012.08.27_22:48_TAI/1d@4h]{mask}';
% 
% Inputs:
%   string mask_series
%   string harp_series
%   string dest_dir
%   opt string mode = 'frame'
%   opt string setup = ''
%   opt string harp_selector = ''
%
% Outputs:
%   (to disk)
%
% See also: track_hmi_harp_movie_driver.sh

fprintf('\n');
fprintf('%s: Path info: now running: %s\n', mfilename, mfilename('fullpath'));
fprintf('%s: Path info: track_hmi_harp_movie_loop is at: %s\n', mfilename, which('track_hmi_harp_movie_loop'));
fprintf('\n');

% require mask_series
if ~exist('mask_series', 'var'),
  error('Require "mask_series" to be defined');
end;
% require harp_series
if ~exist('harp_series', 'var'),
  error('Require "harp_series" to be defined');
end;
% require dest_dir
if ~exist('dest_dir', 'var') || length(dest_dir) == 0,
  error('Require "dest_dir" to be defined and nonempty');
end;
% make frames by default
if ~exist('mode', 'var'),
  mode = 'frame';
end;
% allow for sub-selection of harps at each T_REC
if ~exist('harp_selector', 'var'),
  harp_selector = '';
end;
% setup function: optional, default to none
if exist('setup', 'var') && ~isempty(setup),
  if ~exist(setup, 'file'), error('Supplied setup function (%s) not present', setup); end;
  param_setup = str2func(setup); % turn into a function handle
else,
  setup = '';
  param_setup = @(hooks1)hooks1; % dummy that copies arg thru
end;

% set up NRT mode if mask_series contains _nrt
%   (used internally for WCS and QUALITY)
if ~isempty(strfind(mask_series, '_nrt')),
    hmi_property('set', 'nrt_mode', 1);
end;

% set up hooks -- misc. parameters
hooks1 = struct();
hooks1.mode = mode;
hooks1.tag = track_hmi_clean_runname(mask_series);
hooks1.mp4 = 'clean';
% standard hmi geometry-getter
hooks1.load_meta = @(t)(hmidisk(t, 'geom,vector,quiet'));
% sub-selection of HARPs at each T_REC
hooks1.harp_selector = harp_selector;
% 12 is good for hourly cadence, 20 for 12-min cadence
hooks1.fps = 20;
% 4x spatial subsampling
hooks1.sub_sample = 4;
% label strings
hooks1.labels.title = 'SDO/HMI Tracked AR (HARP)';
hooks1.labels.arname = 'HARPs';
hooks1.labels.archar = 'H';

% allow the hooks to be transformed
if ~isempty(setup), fprintf('%s: Running setup: %s\n', mfilename, setup); end;
hooks = param_setup(hooks1);

% output frames
fn_pat = fullfile(dest_dir, 'harp.%s.png');

% ensure output dir exists
if ~exist(dest_dir, 'dir'), mkdir(dest_dir); end; 

% turn off annoying warnings about NOAA metadata
warning('off', 'hmi_base:no_noaa_info');

% use query engine?
if ~isempty(strfind(mode, 'engine')),
  hmi_property('set', 'jsoc_host', 'server');
end;

% switches from loading images by http versus filesystem
% (see jsoc_trec_loop for more)
if ~isempty(regexpi(getenv('HOST'), 'stanford')),
  smpl = 1; % filesystem
else,
  smpl = -1; % http
end;

dt = 1; % temporal subsampling

% loop which makes the frames
fprintf('%s: Beginning frames.\n', mfilename);
[fn,res] = jsoc_trec_loop(mask_series, [smpl dt 1], ...
                          @track_hmi_harp_movie_loop, ...
                          harp_series, fn_pat, hooks);
fprintf('%s: Completed frames.\n', mfilename);

% shut down query engine
if ~isempty(strfind(mode, 'engine')),
  server_jsoc({'op', 'exit'});
end;

fprintf('%s: Done.\n', mfilename);
return

