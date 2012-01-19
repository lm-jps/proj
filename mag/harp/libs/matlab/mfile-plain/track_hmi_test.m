% script - run the tracker over pieces of HMI data
%
% * This is an ad hoc script, which is frequently modified, to invoke
% the hmi tracker using certain input arguments.  Depending on the 
% setting of a variable `for_real', different input datasets, output
% directories, or parameters may be used.
%
% See also: track_hmi_production
%

% tracker state
global ROI_s ROI_t ROI_rgn

% require for_real variable
if ~exist('for_real', 'var'),
  error('Require "for_real" to be defined');
end;

% the magnetogram data series is always the same
mag_series = 'hmi.M_720s{magnetogram}';
use_mags = true; % can override this

% tracker internal parameters
tracker_params1.kwid    = 1;
tracker_params1.klat    = 4;     % stretch factor along latitude
tracker_params1.tau     = 0.018; % threshold
tracker_params1.tau2    = 0.8;   % tau reduction factor for "try harder"
tracker_params1.active  = 2;     % mask code for "activity"
tracker_params1.mode    = 0;     % plain (0), or laplacian (1) [unsupported]
tracker_params1.final_num  = -2;  % delay in frames before declaring final (<0)
tracker_params1.final_time = 1/8; % delay in days before declaring final
tracker_params1.maprate = 30*log(2); % gives wt, weight for old map
                                     % wt = exp(-maprate * dt [days])
                                     % 30 log(2): wt = 1/2 @ dt = 4*720sec

% handlers for utilities, abstracted from the tracker itself
hooks.load_meta      = @hmidisk;
hooks.save_track     = @hmi_save_track_local;
hooks.init_run       = @hmi_init_run_local;
hooks.loop_hook      = 'hmi_loop_hook_local'; % a script, not a function
hooks.tracker_params = tracker_params1;
hooks.find_roi       = @hmi_find_roi; % uses tracker_params
hooks.first_track    = 1; % assume this is not a restart
hooks.run_name       = 'HMI experiment'; % text tag, to be replaced below

% ensure the dependence checker depfun sees we depend on hooks.loop_hook
if rand(1) > 10,
  hmi_loop_hook_local;
end;

% filename hooks 
hooks.filename.template = strcat(sprintf('/tmp/fr%+03.0f', for_real*10),...
                                 '/Tracks/%s/track%s%s.%s');
hooks.filename.tracknum = '%06d';
hooks.filename.imagenum = '%06d';

if for_real == -1, 
  % trivial test (5 masks, no mags)
  use_mags = false;
  mask_series = 'su_turmon.Marmaskv1_720s[2011.03.31_12:00_TAI/1h]{mask}';

elseif for_real == 0,
  % trivial test (5 masks + mags)
  % mask_series = 'hmi.Marmask_720s[2011.01.15_12:00_TAI/1h]{mask}';
  mask_series = 'su_turmon.Marmaskv1_720s[2011.03.31_12:00_TAI/1h]{mask}';

elseif for_real == 1,
  % longer simple test (125 masks)
  use_mags = false;
  mask_series = 'hmi.Marmask_720s[2011.01.15_12:00_TAI/25h]{mask}';

elseif for_real == 2,
  % even longer test (500 masks)
  use_mags = false;
  mask_series = 'hmi.Marmask_720s[2011.01.15/96h]{mask}';

elseif for_real == 3,
  % long test (1500 masks)
  use_mags = false;
  mask_series = 'hmi.Marmask_720s[2011.01.14/13d]{mask}';

elseif for_real == 3.1,
  % long test (1500 masks)
  use_mags = false;
  mask_series = 'hmi.Marmask_720s[2011.01.14/13d]{mask}';
  hooks.tracker_params.maprate = Inf; % old wt = 0 for any dt > 0

elseif for_real == 3.2,
  % long test (1500 masks)
  use_mags = false;
  mask_series = 'hmi.Marmask_720s[2011.01.14/13d]{mask}';
  hooks.tracker_params.maprate = 60*log(2); % wt = 1/2 @ 2*720s

elseif for_real == 3.3,
  % long test (1500 masks)
  use_mags = false;
  mask_series = 'hmi.Marmask_720s[2011.01.14/13d]{mask}';
  hooks.tracker_params.maprate = 30*log(2); % wt = 1/2 @ 4*720s

elseif for_real == 4,
  % jsoc "production" (3 masks + 3 mags)
  % mask_series = 'hmi.Marmask_720s[2010.09.10_12:00_TAI/40m]{mask}';
  mask_series = 'su_turmon.Marmaskv1_720s[2010.09.10_12:00_TAI/36m]{mask}';
  hooks.save_track = @hmi_save_track_production;
  hooks.init_run   = @hmi_init_run_production;
elseif for_real == 4.1,
  % jsoc "production" (3 masks + 3 mags, add-on to "4.0" above)
  mask_series = 'hmi.Marmask_720s[2010.09.10_12:36_TAI/36m]{mask}';
  hooks.save_track = @hmi_save_track_production;
  hooks.init_run   = @hmi_init_run_production;
  hooks.first_track = 0; % restart
elseif for_real == 4.2,
  % jsoc "production" (3 masks + 3 mags, add-on to "4.1" above)
  mask_series = 'hmi.Marmask_720s[2010.09.10_13:12_TAI/36m]{mask}';
  hooks.save_track = @hmi_save_track_production;
  hooks.init_run   = @hmi_init_run_production;
  hooks.first_track = 0; % restart

elseif for_real == 5,
  % jsoc "production" (125 masks + 125 mags)
  mask_series = 'hmi_test.Marmask_720s[2010.08.09_06:00:00_TAI/25h]{mask}';
  hooks.filename.template = '/tmp/hmi-jsoc-2/Tracks/%s/track%s%s.%s';
  hooks.save_track = @hmi_save_track_production;
  hooks.init_run   = @hmi_init_run_production;

elseif for_real == 6,
  % hard feb 2011 test
  use_mags = false;
  mask_series = 'su_turmon.test_armask_720s[2011.02.01/12d]{mask}';
  hooks.tracker_params.maprate = 30*log(2); % wt = 1/2 @ 4*720s

elseif for_real == 7,
  % hard feb 2011 test
  use_mags = false;
  mask_series = 'su_turmon.test_armask_720s[2011.02.01/28d]{mask}';
  hooks.tracker_params.maprate = 30*log(2); % wt = 1/2 @ 4*720s

elseif for_real == 8,
  %% "feb 14 ar" -- with mags
  % mask_series = 'su_turmon.test_armask_720s[2011.02.10_12:00/12d]{mask}';
  % hooks.filename.template = '/tmp/hmi-jsoc-test/Tracks/%s/track%s%s.%s';

  % new inputs for frame info test, june 2011
  mask_series = 'su_turmon.Marmaskv1_720s[2011.03.07_12:00_TAI/2d]{mask}';

  hooks.tracker_params.tau2 = 0.8; % don't try as hard
  hooks.tracker_params.maprate = 30*log(2);
  hooks.save_track = @hmi_save_track_production;
  hooks.init_run   = @hmi_init_run_production;

elseif for_real == 9,
  % test of the improved "try harder" mechanism
  use_mags = false;
  mask_series = 'su_turmon.test_armask_720s[2011.02.01/3d]{mask}';
  hooks.tracker_params.tau2 = 0.8; % don't try as hard
  hooks.tracker_params.maprate = 30*log(2); % wt = 1/2 @ 4*720s

else,
  disp('Need for_real to be set up, see the script');
  return;
end;

if use_mags,
  jsoc_ds = {mask_series, mag_series};
else,
  jsoc_ds = mask_series;
end;

% text tag for the run
hooks.run_name = mask_series;

% switches from loading images by http versus filesystem
% (see jsoc_trec_loop for more)
if ~isempty(regexpi(getenv('HOST'), 'stanford')),
  smpl = 1; % filesystem
else,
  smpl = -1; % http
end;

% say what we're about to do
if use_mags, use_mags_string = 'with mags'; else use_mags_string = 'no mags'; end;
fprintf('Tracking from "%s" (%s).\nPutting files in "%s".\n', ...
        mask_series, use_mags_string, hooks.filename.template);

[fn,res,fnfail,fnlog]=jsoc_trec_loop(jsoc_ds, [smpl 1 1], 'track_hmi_loop', hooks);

% preserve log file and failures, if any, to /tmp/.../Tracks/jsoc/...
dest_dir = fileparts(sprintf(hooks.filename.template, 'jsoc', 'x', 'x', 'ext'));
track_hmi_postrun_files(dest_dir, hooks.run_name, fnlog, fnfail);

return
