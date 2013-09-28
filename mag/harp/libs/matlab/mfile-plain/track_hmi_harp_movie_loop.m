function res=track_hmi_harp_movie_loop(trec, s_time, im, ds_harp, fn_pat, hooks)
%track_hmi_harp_movie_loop	make quick-look image from mask and HARP
%
% res=track_hmi_harp_movie_loop(trec, s_time, im, ds_harp, fn_pat, hooks)
% * Make a single frame, and/or accumulate a movie, for a single T_REC and
% image, probably within a larger series.
% * jsoc_trec_loop should control this function.  The first three
% arguments are supplied automatically when invoked by jsoc_trec_loop.
% * Input trec is the T_REC corresponding to a mask image -- it is
% used to get image metadata, and query for HARPs.
% * Usage: Call the image-map function
%  jsoc_trec_loop('hmi.Marmask_720s[2010.07.01/1d]', [-1 1 1],...
%                 'track_hmi_harp_movie_loop', ds_harp, fn_pat, hooks);
% Using skip=[-1 1 1] here indicates we want the masks loaded
% fully, don't want to temporally subsample masks, and want begin/end hook
% calls made to this function by jsoc_trec_loop (this is non-optional).
% * If skip(1) is negative, the image-loader will use http to get the
% images; if positive, it will look in the filesystem.  See 
% jsoc_trec_loop and jsoc_imload.
% * If res is returned as empty, this indicates the file was skipped.
% 
% Inputs:
%   string trec;   -- tag naming a mask
%   real s_time;   -- matlab time number
%   real im(m,n);  -- the contents of fn
%   string ds_harp
%   string fn_pat
%   struct hooks
% 
% Outputs:
%   int res
% 
% See Also:  track_hmi_harp_movie_script

% Written by Michael Turmon (turmon@jpl.nasa.gov) in 2012
% Copyright (c) 2010.  All rights reserved.  

% videowriter object
persistent Movie

% 
% Error checking
% 

if all(nargin  ~= [6]), error ('Bad input arg number'); end;
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end;

% write frame?
write_frame = ~isempty(strfind(hooks.mode, 'frame'));
% write movie?
write_movie = ~isempty(strfind(hooks.mode, 'movie'));
% tag for filename
if isfield(hooks, 'tag'), tag = hooks.tag; else tag = 'no_tag'; end;
% FPS?
if isfield(hooks, 'fps'), fps = hooks.fps; else, fps = 12; end;
% sub_sample?
if ~isfield(hooks, 'sub_sample'), hooks.sub_sample = 1; end;
% harp qualifier?
if ~isfield(hooks, 'harp_select'), hooks.harp_select = ''; end;

% 
% Setup/teardown
% 
if isnumeric(trec),
  % begin/end sentinel
  res = 0; % set it now, for consistency (non-empty)
  if trec < 0,
    % Loop-begin setup
    if write_movie,
      % MovieType = 'Uncompressed AVI';
      MovieType = 'Motion JPEG AVI';
      fn = fullfile(fileparts(fn_pat), ['harp.' tag '.avi']);
      Movie = VideoWriter(fn, MovieType);
      Movie.FrameRate = fps; 
      Movie.Quality = 100; % only for motion jpeg
      open(Movie);
    end;
    return;
  else,
    % Loop-end teardown
    if write_movie,
      % preserve filename
      fn = fullfile(Movie.Path, Movie.Filename);
      close(Movie);
      % without the clear, this segfaults on the stanford system
      clear Movie;
      % optionally run ffmpeg to transcode to mp4
      if isfield(hooks, 'mp4'),
        opts = struct(); % by default, nothing
        % optionally clean the avi
        if ischar(hooks.mp4) && strcmp(hooks.mp4, 'clean'),
          opts.clean = 1; % value is unimportant
        end;
        status = avi2mp4(fn, opts, 'qscale', 4, 'v', '0');
        if status ~= 0,
          fprintf('Error: Failed to convert avi to mp4!\n');
        else,
          fprintf('Converted avi to mp4.\n');
        end;
      end;
    end;
    return;
  end;
  % it's not a hook we know
  error('Bad input args: trec is numeric');
end;

if write_movie && isempty(Movie),
  error('Movie is empty, did you initialize (set jsoc_trec_loop skip(3) nonzero)?');
end;

%
% Computation
% 

% output frame name
fn_out = sprintf(fn_pat, trec);

% categorical color map for movies
% cm = [1,1,1; 0,0,0; 0.4,0.4,0.4; 0.7,0.7,0.7; prism2(40)]; % old -- not as good
white = [1 1 1];
cm = [1.0*white; ...
      1,1,0; ...
      0.0*white; ...
      0.4*white; ...
      0.7*white; ...
      colormap_protovis];
cm(end-3,:) = []; % this is too close to yellow (cm(2,:) overlay)

% parameters to get from each HARP
params = {'key', ...
          'HARPNUM,CRPIX1,CRPIX2,CRSIZE1,CRSIZE2,H_MERGE,H_FAINT,NPIX,T_FRST1,T_LAST1'};

% image metadata
[geom,msg] = hooks.load_meta(trec);
if ~isempty(msg),
  res = msg;
  return;   % this will skip it
end;

% harp selector, if provided
if isempty(hooks.harp_selector),
  selector = '';
else,
  selector = sprintf('[?%s?]', hooks.harp_selector);
end;
% get HARP numbers and HARP metadata
[res,msg] = rs_list(sprintf('%s[][%s]%s', ds_harp, trec, selector), params);
if ~isempty(msg),
  res = msg; % skip frame, and indicate error
  return;
end;
if res.count == 0,
  % will be taken care of elsewhere
  harps = [];
  Nharp = 0;
else,
  harps = jsoc_cell2struct_keys(res.keywords);
  Nharp = length(harps.harpnum);
end;

% get bitmaps for each HARP
bitmaps = cell(1,Nharp);
for inx = 1:Nharp,
    h = harps.harpnum(inx);
    [im1,err] = jsoc_imload([ ds_harp '[' int2str(h) ']' ], trec, 'bitmap');
    if ~isempty(err),
      res = err; % skip frame, and indicate error
      return;
    end;
    bitmaps{inx} = im1;
end;

frame = track_hmi_movie_frame(s_time, geom, im, harps, bitmaps, ...
                              hooks.sub_sample, hooks.labels);

% write frame if desired
if write_frame,
  imwrite(frame, cm, fn_out, 'png');
end;

% write movie if desired
if write_movie,
  writeVideo(Movie, im2frame(frame, cm));
end;

% fprintf('Finished %s\n', fn_out);

res = 0;
return;
end


