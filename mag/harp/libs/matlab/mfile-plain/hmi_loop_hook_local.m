%hmi_loop_hook_local	script, within-loop dump of results
% 
% * Runs within the tracker_loop as a script, currently just dumps the
% full-disk ROI image.
% 
% See Also:

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 31 Aug 2010
% Copyright (c) 2010.  All rights reserved.

%
% Computation
% 

% we use fp
if exist('fp', 'var'),
  warning('%s: found "fp" variable in workspace, will overwrite', mfilename);
end;

%% segmentation file
fn_s = sprintf(hooks.filename.template, 's', '-fd-segment.', s_imagenumS, 'fits.gz');
[junk1,junk2] = mkdir(fileparts(fn_s)); % suppress warning if exists
savefitsa(im, fn_s, 8, 1, 0, 255);

%% boundary file
if exist('roimap_p', 'var'),
  roimap_disk = roimap_p; % use the posterior one if it exists
else,
  roimap_disk = roimap;
end;
roimap_disk(im == 0) = NaN;

fn_b = sprintf(hooks.filename.template, 'rb', '-fd-bound.', s_imagenumS, 'fits.gz');
[junk1,junk2] = mkdir(fileparts(fn_b)); % suppress warning if exists
savefitsa(roimap_disk, fn_b, 8, 1, 0, 255);

%% also save original boundary file, if there were two
if exist('roimap_p', 'var'),
  roimap_disk = roimap; % now this is the ORIGINAL roimap
  roimap_disk(im == 0) = NaN;

  fn_b0 = sprintf(hooks.filename.template, 'rb', '-fd-bound0.', s_imagenumS, 'fits.gz');
  savefitsa(roimap_disk, fn_b0, 8, 1, 0, 255);
end;

clear roimap_disk disk_mask junk1 junk2

%% geometry, tag, time
fn_m = sprintf(hooks.filename.template, 'xfd', '-fd-meta.', s_imagenumS, 'mat');
save(fn_m, 's_geom', 'fn', 's_time', 'Nr');

%% instant file
fn_i = sprintf(hooks.filename.template, 'ifd', '-fd-instant.', s_imagenumS, 'instant');
fp = fopen(fn_i, 'w');
fprintf(fp, '#\n# .instant metafile\n');
fprintf(fp, '# linking full-disk images with their metadata\n#\n');
% fprintf(fp, 'region %s\n', fnr);
fprintf(fp, 'segment %s\n',  fn_s);
fprintf(fp, 'boundary %s\n', fn_b);
if exist('roimap_p', 'var'),
  fprintf(fp, 'boundary0 %s\n', fn_b0); % the original boundary
end;
fprintf(fp, 'metadata %s\n', fn_m);
fclose(fp);

%% update instant cat-file
% fulldisk instant cat file (.../Tracks/track-fd-instant.cat)
% NB: must match above
fn_ic = sprintf(hooks.filename.template, '', '-fd-instant', '', 'cat');
fp = fopen(fn_ic, 'a');
fprintf(fp, '%s %s\n', fn_i, s_timeS);
fclose(fp);

clear fn_b fn_b0 fn_m fn_s fn_i fn_ic
clear fp

%% -- could save other things, do plots, etc. --

return;