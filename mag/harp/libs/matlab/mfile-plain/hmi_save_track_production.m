function err=hmi_save_track_production(rois, tid, birth, status, run_name, tmpl)
%hmi_save_track_production	save a track for later JSOC ingestion
% 
% err=hmi_save_track_production(rois, tid, birth, tmpl)
% * Saves a set of chips representing a track.
% 
% Inputs:
%   struct rois(nt)
%   int tid
%   real birth
%   string status in { 'merged' 'pending' 'final' }
%   string run_name
%   string tmpl
% 
% Outputs:
%   int err
% 
% See Also:

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 31 Aug 2010
% Copyright (c) 2010.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [6]), error ('Bad input arg number'); end;
%if all(nargout ~= [0 1 2 3 5]), error ('Bad output arg number'); end;

%
% Computation
% 

save_movie_information = true;

% edit the ROIs to remove times where the ROI was not found in the image
roi_present = [rois.tag] >= 0;
rois = rois(roi_present);

fprintf('closing track %d of %d chips [%s]\n', tid, length(rois), status);

% FIXME: should allow setting this another way
harpDir = sprintf('-%06d/', tid);
Nf = length(rois);

% save chip-list to a mat-file
% (jsoc ingestion does not use this, but leaving it in allows movies 
% to be made)
if save_movie_information,
  fn = sprintf(tmpl, 'mask', '-mask.', num2str(tid, '%06d'), 'mat');
  [junk1,junk2] = mkdir(fileparts(fn)); % suppress warning if exists
  save(fn, 'rois', 'tid');
end;

%% save track status
% 1: global status
if strcmp(status, 'final'),
  % cumulative finalized-tracks list
  fn = sprintf(tmpl, 'jsoc', '-final', '', 'txt');
  [junk1,junk2] = mkdir(fileparts(fn)); % suppress warning if exists
  fp = fopen(fn, 'a'); % note, append
  fprintf(fp, '%d\n', tid);
  fclose(fp);
  % newly-finalized-tracks list
  fn = sprintf(tmpl, 'jsoc', '-new', '', 'txt');
  [junk1,junk2] = mkdir(fileparts(fn)); % suppress warning if exists
  fp = fopen(fn, 'a');
  fprintf(fp, '%d\n', tid);
  fclose(fp);
end;
if strcmp(status, 'pending'),
  % pending-tracks list
  fn = sprintf(tmpl, 'jsoc', '-pending', '', 'txt');
  [junk1,junk2] = mkdir(fileparts(fn)); % suppress warning if exists
  fp = fopen(fn, 'a'); % note, append
  fprintf(fp, '%d\n', tid);
  fclose(fp);
end;
% 2: per-directory status
fn = sprintf(tmpl, 'jsoc', harpDir, 'track-status', 'txt');
[junk1,junk2] = mkdir(fileparts(fn)); % suppress warning if exists
fp = fopen(fn, 'a'); % note, append
if fp < 0,
  warn('TID %d: failed to write %s', tid, fn);
  err = 1;
  return;
end;
fprintf(fp, '%s\t%s\t%s\n', status, run_name, datestr(now, 30));
fclose(fp);


%% save chips as fits
for f = 1:Nf,
  roi1 = rois(f);
  % save bitmap
  fn = sprintf(tmpl, 'jsoc', harpDir, sprintf('patch-%06d', f), 'fits.gz');
  [junk1,junk2] = mkdir(fileparts(fn)); % suppress warning if exists
  savefitsa(double(full(roi1.chip)), fn, 8);
end;

%% save stats as ascii
fn = sprintf(tmpl, 'jsoc', harpDir, 'track-stats', 'txt');
[junk1,junk2] = mkdir(fileparts(fn)); % suppress warning if exists
statsAll = reshape([rois(:).stats], [], Nf)';
save(fn, 'statsAll', '-ascii', '-double');

%% save the lat/lon boundary as a separate file
% convert to length-Nf structure array for clarity in field references
statSt = roi_stats2struct(statsAll);
% rotational constants
obs_v = 0.9865; % observer's angular velocity in deg/day
coefs = [14.6,-2.2,0]; % solar sidereal rotation, deg/day, powers of sin2(lat)
% rotation rate for this HARP
% OLD: lat_avg = median([statSt.ar_fwt_lat]); % median of flux-weighted means
lat_avg = median(([statSt.rgn_min_lat] + [statSt.rgn_max_lat])/2); % median of bbox centers
sin2lat = sind(lat_avg).^2; % sin in degrees
deg_day = -obs_v + coefs(1) + coefs(2)*sin2lat + coefs(3)*sin2lat.^2; % for now
% rotation angle
times = [rois.time]; % this is in days
times = times - times(1); % relative to start 
rotation = deg_day * times; % expected rotation in degrees
% centered longitudes
lontop = [statSt.rgn_max_lon] - rotation;
lonbot = [statSt.rgn_min_lon] - rotation;
% origin and width of lon bounding box come from centered longitudes 
lontop_max = max(lontop);
lonbot_min = min(lonbot);
lonOrigin = mean([lontop_max lonbot_min]);
lonWidth = lontop_max - lonbot_min;
% compose the lon bounding boxes, and un-center them
lontop_new = lonOrigin + (lonWidth*0.5) + rotation;
lonbot_new = lonOrigin - (lonWidth*0.5) + rotation;
% find bracketing lat, which is unchanging over time
lat_max = max([statSt.rgn_max_lat]);
lat_min = min([statSt.rgn_min_lat]);
% write out the lat/lon bounding boxes
%   (lats are all the same)
LatLonBB = [repmat(lat_min, [Nf 1]) lonbot_new(:) ...
            repmat(lat_max, [Nf 1]) lontop_new(:)];

%% save pixel and lat/lon BBs as ascii
fn = sprintf(tmpl, 'jsoc', harpDir, 'track-frame', 'txt');
[junk1,junk2] = mkdir(fileparts(fn)); % suppress warning if exists
fp = fopen(fn, 'w');
if fp < 0,
  warn('TID %d: failed to write %s', tid, fn);
  err = 1;
  return;
end;
for f = 1:Nf,
  % size of outline chip
  PixSzs    = sprintf('%d ', size(rois(f).chip));
  % upper-left corner of outline chip (.fits above)
  PixBBs    = sprintf('%d ',   rois(f).coords);
  % lat/lon bounding box vs. time
  LatLonBBs = sprintf('%.4f ', LatLonBB(f,:));
  % add a line to the file
  fprintf(fp, '%d %s %s %s %s %.4f\n', tid, rois(f).fn, ...
          PixSzs, PixBBs, LatLonBBs, deg_day);
end
fclose(fp);

% keyboard

err = 0;
return;
end
