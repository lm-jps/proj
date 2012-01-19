%
% Script, make movie showing track boundaries and segmentations
%
% This script takes in files produced by track_hmi_test.  It uses
% file_loop_lst to iterate over these files, and for each one, the
% inner-loop function track_hmi_movie_loop is called.
%
% See also: track_hmi_test (produces files this routine needs)
%

% mode defaults to 1
if ~exist('mode', 'var'),
  mode = 1;
  fprintf('Defaulting to mode = %d (display resolution)\n', mode);
end;

% colormap: white, black, gray, then a series of easily distiguishable colors
% the text rendering in the moviemaker relies upon the first 3 colors being
% this way
cm = prism2(40); 
cm = [1,1,1; 0,0,0; 0.4,0.4,0.4; 0.7,0.7,0.7; cm];

% inner loop command
cmd = @track_hmi_movie_loop;

rootd = sprintf('/tmp/fr%+03.0f', for_real*10);
fnIC = strcat(rootd, '/Tracks/track-fd-instant.cat');
fnTs = strcat(rootd, '/Tracks/mask/track-mask.%s.mat');
    
%% source list
if for_real == 0,
  key = '1h';
elseif for_real == -1,
  key = '1h-nomag';
elseif for_real == 1,
  key = '25h';
elseif for_real == 2,
  key = '96h';
elseif for_real == 3,
  key = '13d';
elseif for_real == 3.1,
  key = '13d-inf';
elseif for_real == 3.2,
  key = '13d-60';
elseif for_real == 3.3,
  key = '13d-30';
elseif for_real == 4.0,
  key = 'jsoc-30m-a';
elseif for_real == 4.1,
  key = 'jsoc-30m-b';
elseif for_real == 4.2,
  key = 'jsoc-30m-c';
elseif for_real == 6,
  key = 'feb-2011-12d';
elseif for_real == 7,
  key = 'feb-2011-28d';
elseif for_real == 8,
  % rootd = sprintf('/tmp/hmi-jsoc-feb14');
  % rootd = sprintf('/tmp/hmi-jsoc-test');
  % fnIC = strcat(rootd, '/Tracks/track-fd-instant.cat');
  % fnTs = strcat(rootd, '/Tracks/mask/track-mask.%s.mat');
  key = 'mar-2011-2d';
elseif for_real == 9,
  key = 'feb-2011-3d';
else,
  disp('Did not get a good for_real');
  return;
end;

% temporal subsampling
dt = 1; % usual case
% dt = 10; % for example

%% make movie
if mode == 1,
  % set up for 1024^2 frames
  clear res;
  Mname = sprintf('/tmp/hmi-patch-%s-all_%%d.avi', key);
  [fn,res]=file_loop_cat(fnIC, [0 dt 1], cmd, 4, cm, fnTs, Mname);

elseif mode == 2,
  % 512^2 frames
  clear res;
  Mname = sprintf('/tmp/hmi-patchS-%s-all_%%d.avi', key);
  [fn,res]=file_loop_cat(fnIC,[0 dt 1], cmd, 8, cm, fnTs, Mname);

elseif mode == 3,
  % 256^2 frames, not stored in a disk file as-we-go, but saved later
  clear M;
  [fn,M]=file_loop_cat(fnIC,[0 dt 1], cmd, 16, cm, fnTs);
  movie2avi(M, '/tmp/hmi-patch-SS.avi', 'fps', 8);

end;

