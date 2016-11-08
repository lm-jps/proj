function frame=track_hmi_movie_frame(t, geom, mask, harps, bitmaps, smpl, labels)
%track_hmi_movie_frame	make one frame of a tracker movie
% 
% frame=track_hmi_movie_frame(t, geom, mask, harps, bitmaps, smpl, labels)
% * Render one frame of a HARP diagnostic movie.  The output is a single
% indexed image.  This uses code similar to track_hmi_movie, but is 
% designed to retain only the frame-making parts, leaving the I/O to 
% a separate driver routine.
% * This routine was originally intended to be independently callable, 
% but too much extra information is needed to render the frame.  We need
% HARP information (harps structure) and data segments (bitmaps cell).
% Thus, breaking it off helped modularity, but this routine is locked
% in to its place in the pyramid.
% * You can, however, call the function with harps and bitmaps empty.
% * The integer smpl gives spatial downsampling from mask to frame. 
% For example, smpl=4 implies downsampling by a factor of 4 in each 
% direction, yielding roughly 1Kx1K images for HMI.
% * The labels structure contains strings used for text labels in the frame.
% 
% Inputs:
%   real t;       -- datenum time
%   real geom(5)  -- a geom
%   real mask(m,n)
%   struct harps(1) or []
%   cell bitmaps(nr) of real
%   int smpl;
%   struct labels
% 
% Outputs:
%   int frame(mp,np)
% 
% See Also: track_hmi_harp_movie_loop

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 03 Oct 2012
% Copyright (c) 2012.  All rights reserved.  

% 
% Error checking
% 
% if all(nargin  ~= [7]), error ('Bad input arg number'); end;
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end;
if length(smpl) ~= 1, error('Need one number in smpl'); end;
if isempty(harps),
    harps = struct('harpnum', [], ...
                   'crpix1',  [], ...
                   'crpix2',  [], ...
                   'crsize1', [], ...
                   'crsize2', [], ...
                   't_frst1', [], ...
                   't_last1', [], ...
                   'h_merge', [], ...
                   'h_faint', [], ...
                   'npix',    []);
end;

nR = length(harps.harpnum);

rgn = [harps.crpix1' ...
       harps.crpix2' ...
       harps.crpix1' + harps.crsize1' - 1 ...
       harps.crpix2' + harps.crsize2' - 1];
% down-sampled region boundary
rgnS = floor((rgn-1) / smpl) + 1;

rgnAge = (t - harps.t_frst1); % 1 if new
rgnPad = (t > harps.t_last1) - (t < harps.t_frst1); % -1/0/+1
rgnMrg = (harps.h_merge >  0); % 1 if merged
rgnTag = double(harps.h_faint == 0); % 1 = normal, 0 = faint...
rgnTag((harps.h_faint == 1) & isnan(harps.npix)) = -1; % ...-1 = placeholder
rgnTid = (harps.harpnum)';

SkipEmptyBox = false; % show box around ROI for placeholders?

frame_num = NaN; % unknown
prior_t = NaN; % unknown

mask = mask(1:smpl:end,1:smpl:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NOAA Regions

% load NOAA regions -- allow extrapolation
NOAAar = hmi_noaa_info_interp(t, 'extrap');
% project into our disk geometry
[NOAAx, NOAAy, NOAAz] = hmi_latlon2image([NOAAar.latitudehg ]', ...
                                         [NOAAar.longitudecm]', geom);
NOAAx = round(NOAAx/abs(smpl));
NOAAy = round(NOAAy/abs(smpl));
NOAAz = round(NOAAz/abs(smpl));
NOAAv = (NOAAz >= 0); % on-disk filter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Disk Geometry

% disk latlon-center in image coordinates -- ensure it is integral.
%   due to b-angle, latlon-center != X0,Y0.  if the two are much
%   different, the "plot numbers along the equator" scheme wouldn't work.
[CTRx, CTRy, CTRz] = hmi_latlon2image(0, 0, geom);
CTRx = round(CTRx / abs(smpl));
CTRy = round(CTRy / abs(smpl));
% poles in image coordinates
[NPx, NPy, NPz] = hmi_latlon2image([90;-90], [0;0], geom);
NPx = round(NPx/abs(smpl));
NPy = round(NPy/abs(smpl));
NPv = (NPz >= 0); % visible on-disk

% p-angle
p0 = geom(5);

% discretized p0 -- 0, 90, 180, etc. -90 maps to 270.
p0d = mod(round(mod(p0, 360) / 90) * 90, 360);
% rotation by p0
P0_mat  = [cosd(p0)  sind(p0); -sind(p0) cosd(p0)];
% pole-pole unit vector, and its complement (see below for interpretation)
%   used to transform from image coords to sun-centric coordinates
%   we need sun-centric coordinates in two places:
%   (1) To put NOAA-AR labels on the disk, because they need to be squashed
%       toward the equator, and tagged with arrows as northern or southern;
%   (2) To decide northern/southern hemisphere for the NOAA AR and HARP 
%       "legends" at the right edge of the frame
U_pole = P0_mat(2,:)';
U_equa = P0_mat(1,:)';
% image rotation depends on p0d, not p0 -- just flip, don't interpolate
SPIN = round(p0d / 90);
SPIN = 1-SPIN; % chosen so that rot90(im,SPIN) orients im correctly
if nR > 0,
  % find bounding box centroid in image coordinates
  %   note: box centroid can be off-disk, so lat/lon may not work
  centroidIMG = [mean(rgnS(:,[1 3]),2) mean(rgnS(:,[2 4]), 2)];
  % transform to sun-centric coordinates by rotating backwards by p0
  %   E = along equator; P = pole-to-pole
  %   subtract sun center, rotate by -p0
  centroidEP = bsxfun(@minus, centroidIMG, [CTRx CTRy]) * P0_mat';
  % is HARP above or below equator (above if polar distance > 0)
  UPPER = centroidEP(:,2) > 0;
else
  UPPER = [];
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Colors, margins

% margins
MRGN = 0.01; % [0,1] units
MRGNPIX = floor(MRGN*size(mask, 1))-1; % pixels: extra (-1) is empirical

% original region map:
%   offdisk = NaN, ondisk = 0, regions = 1..Nr
% remapped into new integers made to align with the color map
HILITE = 1;  % color index for text, overlays, etc. (intended as white)
CONTRA = 2;  % contrast for text, overlays (yellow)
LOLITE = 3;  % complementary to HILITE (black)
NEUTA  = 4;  % boring (dark gray)
NEUTB  = 5;  % boring (brighter gray)
FRSTROI= 6;  % first ROI color
TCOLOR  = [LOLITE NEUTA NEUTB HILITE]; % text, dark-to-bright
TCOLOR2 = [LOLITE NEUTA NEUTB CONTRA]; % text/alt, dark-to-bright

% load the fonts for text
[M,N] = size(mask); % it has been downsampled
if M < 400, 
  DFONT = text_on_image_font('Menlo.ttf', [], [16 0]);
  SFONT = text_on_image_font('ProggyClean_font.m', [], [10 0 0 1]);
else, 
  DFONT = text_on_image_font('Menlo.ttf', [], [18 0 0 1]);
  SFONT = text_on_image_font('ProggyClean_font.m', [], [10 0 0 1]);
end;

% Nt_show = 40; % old (2010) contrasty prism2 colors
Nt_show = 17; % new (2012) protovis colors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make a full-disk image with segmentation, boxes, box-labels superimposed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 0: Put the blobs into the mask-size image
im1 = zeros(size(mask));
im1(mask == 0) = LOLITE; % off-disk
im1(mask >  0) = NEUTA;  % on-disk
for inx = 1:nR,
  % value for on-blob region
  imval = FRSTROI + rem(rgnTid(inx)-1, Nt_show);
  % indicator for on-blob region, plug active in later
  new = bitmaps{inx}(1:smpl:end,1:smpl:end) > 32; 
  % box to index into the image
  c1 = floor(harps.crpix1(inx)/smpl);
  c2 = floor(harps.crpix2(inx)/smpl);
  inx1 = c1 + [0:size(new,1)-1];
  inx2 = c2 + [0:size(new,2)-1];
  % image chunk
  now = im1(inx1,inx2); % may have overlays already
  now(new) = imval;
  % plug back in
  im1(inx1,inx2) = now;
end;

% 1: (we do not have have both posterior and prior boundary)

% 2: active are LOLITE
im1(mask == 2) = LOLITE;

% 3: box around each region
% Note: region corners are in the ims() coordinates, unflipped/untransposed
rgnS = floor((rgn-1) / smpl) + 1;
if false,
  % (double lines disabled)
  % L-R and T-B hairlines offset 0 and 1
  tk1 = [1 1 3 3]; % four hairlines...
  tk2 = tk1 + 1;
  del = [0 1 0 1]; % ...offset by 0, 1, 0, 1
  flag = 7; % size of merged-region panel
else,
  % hairlines not offset
  tk1 = [1 3]; % two hairlines
  tk2 = tk1 + 1;
  del = [0 0]; % not offset
  flag = 3;
end;
[M,N] = size(im1);
for r=1:nR,
  if SkipEmptyBox && rgnTag(r) < 0, continue; end; % no enclosing box for placeholder regions
  % throughout, have to clip to range of indexes to prevent trouble
  im1(rgnS(r,1):rgnS(r,3), range(rgnS(r,tk2) + del, 1, N)) = HILITE;
  im1(range(rgnS(r,tk1) + del, 1, M), rgnS(r,2):rgnS(r,4)) = HILITE;
  if rgnTag(r) < 0 || rgnPad(r) ~= 0,
    % note: this case is sometimes skipped entirely (above)
    % dotted box for a placeholder or padding region
    im1(rgnS(r,1):2:rgnS(r,3), range(rgnS(r,tk2) + del, 1, N)) = NEUTA;
    im1(range(rgnS(r,tk1) + del, 1, M), rgnS(r,2):2:rgnS(r,4)) = NEUTA;
  elseif rgnTag(r) == 0,
    % dotted horizontal lines for a tried-harder region
    im1(range(rgnS(r,tk1) + del, 1, M), rgnS(r,2):2:rgnS(r,4)) = NEUTA;
  end;
  % flag a merged region (clip to image dims)
  im1(range(rgnS(r,1) + [1:flag          ], 1, M), ...
      range(rgnS(r,2) + [1:flag*rgnMrg(r)], 1, N)) = HILITE;
  % flag a new region (clip to image dims)
  im1(range(rgnS(r,1) - [1:flag                     ], 1, M), ...
      range(rgnS(r,2) - [1:flag*double(rgnAge(r)==0)], 1, N)) = HILITE;
end;

% NOAA AR text labels: along equator
% ARROW = up-arrow char in SFONT; down-arrow is ARROW+1
ARROW = 143; 
% NOAAdist: increases going to north pole
% NOAAequa: increases along direction of sunspot motion
NOAAdist = bsxfun(@minus, [NOAAx NOAAy], [CTRx CTRy]) * U_pole;
NOAAequa = bsxfun(@minus, [NOAAx NOAAy], [CTRx CTRy]) * U_equa;
NOAAsquash = @(x)(sign(x).*sqrt(abs(x))); % squash AR location toward equator
NOAAtxt = cellstr([int2str([NOAAar(NOAAv).regionnumber]') ...
                   char(ARROW + double(NOAAdist(NOAAv) < 0))]);

% rotate (equator, polar) distance by p0 to put onto image (which is rotated by p0)
%   The complex construction below rotates by P0 and offsets onto the disk center,
%   then pads by "0-SPIN" to control rotation of the numerals themselves.
%   After numerals are put on here, the rot90 by +SPIN below leaves them vertical.
% C_OFF: pixel offset given to the 5-number AR label (empirical, for centering)
%   SPIN = 1 (p0=0) needs it; SPIN = -1 (p0=+/-180) and SPIN=-2 (p0=-90) do not.
%   I don't know about SPIN = 0 (p0 = 90).
if SPIN == 1, C_OFF = 30; else C_OFF = 0; end;
if nnz(NOAAv) > 0,
  % "*" below can fail if nothing is visible
  im1 = text_on_image(im1, NOAAtxt, [], CONTRA, ...
                      [bsxfun(@plus, ...
                              [NOAAequa(NOAAv)-C_OFF 3*NOAAsquash(NOAAdist(NOAAv))] * P0_mat, ...
                              [CTRx CTRy]), ...
                      zeros(nnz(NOAAv),1)-SPIN], ...
                      SFONT);
end;

% place NOAA AR markers on the ARs
im1 = marker_on_image(im1, [NOAAx(NOAAv) NOAAy(NOAAv)], ...
                      '+', [15 3], [], CONTRA);

% HARP-number labels: scattered within the image
%   (these are put on after NOAA ARs, to over-write them)
%   Once again, using the same logic as above, the numerals themselves are
%   rotated by -SPIN by text_on_image to compensate for later rotation by 
%   +SPIN by the global rotation below.
if nR > 0,
  harp_lbl_locn = [rgnS(:,1:2)+1 zeros(nR,1)-SPIN];
  % play: harp_lbl_locn = [bsxfun(@minus, rgnS(:,3:4), [9 5*7]) zeros(nR,1)-SPIN];
  % allow to skip labeling the placeholder boxes
  if SkipEmptyBox, inx = (rgnTag >= 0); else, inx = true(size(rgnTid)); end;
  im1 = text_on_image(im1, cellstr(int2str(rgnTid(inx))), ...
                      [], HILITE, harp_lbl_locn(inx,:), SFONT);
  clear inx
end;

% Label the poles, if they're visible
NP_txt = {'N'; 'S'};
NPspc = [-20;10]; % offset text from marker
if any(NPv),
  im1 = marker_on_image(im1, [NPx(NPv) NPy(NPv)], 'x', [15 3], [], HILITE);
  im1 = text_on_image(im1, NP_txt(NPv), [], HILITE, ...
                      [NPx(NPv) NPy(NPv)+NPspc(NPv) zeros(nnz(NPv),1)-SPIN], SFONT);
end;

% keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make the peripheral text overlays
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We have to rotate the image into final form as-presented (North up)
% before we put the text labels on.
% By contrast, the bounding boxes were done in original
% Matlab coordinates, since they were originally defined that way.

% rotate the image up
im1 = rot90(im1, SPIN);
M = size(im1, 1);

% make space in im1 for overlay
if     M < 400, pad = [32 96]; 
elseif M < 800, pad = [ 0 96]; 
else            pad = [ 0 96];
end;
% imT is the image-with-text
imT = [zeros(M, pad(1))+LOLITE im1 zeros(M, pad(2))+LOLITE];

% Time overlay -- upper-left corner
t1code = datestr(t, 'yyyy/mm/dd');
t2code = datestr(t, 'HH:MM');
if prior_t > 0,
  t3code = sprintf('dt = %s', datestr(t - prior_t, 'dd+HH:MM'));
elseif ~isnan(prior_t),
  t3code = 'first frame';
else,
  t3code = '';
end;
% frame number overlay
if ~isnan(frame_num),
  fnumT = num2str(frame_num, '%06d');
else
  fnumT = '';
end;
% initial space makes it play nicer with qt player on mac
txt = {{labels.title, t1code, t2code, ' ', t3code, ' ', fnumT}};
imT = text_on_image(imT, txt, [], TCOLOR, [MRGN MRGN], DFONT);

% if p-angle is not near a multiple of 90, the image will be tilted
%   note: rem, not mod
if abs(rem(p0 - p0d, 360)) > 1,
  txt = {{sprintf('P0 = %.1f ', p0)}};
  imT = text_on_image(imT, txt, [], TCOLOR, [0.10 MRGN], DFONT);
end;

% Region-count overlay -- lower-left corner
% nR is the number of tracks, accounting for merges
txt = {{...
    sprintf('%s = %d', labels.archar, nR), ...
    ' ', ...
    sprintf('! = %d (new)',         nnz(rgnAge == 0)), ...
    sprintf('+ = %d (merge)',       nnz(rgnMrg >  0)), ...
    sprintf('( = %d (pad before)',  nnz(rgnPad <  0)), ...
    sprintf(') = %d (pad after)',   nnz(rgnPad >  0)), ...
    sprintf('~ = %d (use past)',    nnz(rgnTag == 0)), ...
    sprintf('? = %d (placeholder)', nnz(rgnTag <  0)), ...
      }};
% put region-count overlay on
imT = text_on_image(imT, txt, [], TCOLOR, [1-5*MRGN MRGN], DFONT);
% explanatory note
txt = {{...
    sprintf('%s: numbered boxes; active region colored', labels.arname), ...
    'NOAA ARs: crosses; numerical label shifted to near equator'}};
imT = text_on_image(imT, txt, [], TCOLOR, [1-MRGN MRGN], SFONT);

% Insert HARP ID numbers into frame legend
NC_max = floor(size(imT, 1) / (2*19)); % max track ID's in one column
if nR > 0,
  % label regions
  rgnID_t = num2str(rgnTid, '%04.0f'); % text format
  rgnID_x = blanks(nR)';
  rgnID_x(rgnTag == 0) = '~'; % indicate try-harder tracks
  rgnID_x(rgnTag <  0) = '?'; % indicate placeholder tracks
  rgnID_x(rgnMrg >  0) = '+'; % indicate merges
  rgnID_x(rgnAge == 0) = '!'; % indicate new tracks
  rgnID_x(rgnPad <  0) = '('; % indicate pad-before tracks
  rgnID_x(rgnPad >  0) = ')'; % indicate pad-after tracks
  rgnID_t = [rgnID_x rgnID_t];
  % define text blocks
  TBu = columnize(rgnID_t( UPPER,:), NC_max, 1);
  TBd = columnize(rgnID_t(~UPPER,:), NC_max, 0);
  imT = text_on_image(imT, strvcat({labels.arname; TBu}), [], TCOLOR, [MRGN   0.97], DFONT);
  imT = text_on_image(imT, TBd, [], TCOLOR, [1-MRGN 0.97], DFONT);
  % put on the per-track color key
  pix_per_line = 19; % for DFONT
  WID = 20; % width of the color key
  % upper key
  nline = nnz(UPPER);
  track_colors = FRSTROI + rem(rgnTid(UPPER)-1, Nt_show);
  track_block = kron(track_colors, ones(pix_per_line, WID));
  imT(MRGNPIX+pix_per_line+[1:(pix_per_line*nline)],(end-1)-[1:WID]) = track_block;
  % lower key
  nline = nnz(~UPPER);
  track_colors = FRSTROI + rem(rgnTid(~UPPER)-1, Nt_show);
  track_block = kron(track_colors, ones(pix_per_line, WID));
  imT([end-(pix_per_line*nline)+1:end]-MRGNPIX,(end-1)-[1:WID]) = track_block;
else,
  imT = text_on_image(imT, strvcat({labels.arname; '(none)'}), [], TCOLOR, [MRGN   0.97], DFONT);
end;

% Insert NOAA AR numbers into frame legend
if nnz(NOAAv) > 0,
  NOAAlat = [NOAAar.latitudehg]; % all lat's
  NOAAtid = [NOAAar.regionnumber]';
  NOAAspt = [NOAAar.spotcount]';
  % put on northern tags
  NOAAinx = NOAAv &  (NOAAlat(:) >= 0); % valid and north
  NOAAtid1 = NOAAtid (NOAAinx); % track id's
  NOAAspt1 = NOAAspt (NOAAinx); % spot counts
  NOAAeq1  = NOAAequa(NOAAinx); % along-equator values
  [junk,NOAAup] = sort(NOAAeq1, 1, 'ascend'); % sort (valid&north) by equa
  NBu = num2str(NOAAtid1(NOAAup), '%05.0f');
  % NBs = num2str(NOAAspt1(NOAAup), ':%2d'); % remove these
  imT = text_on_image(imT, strvcat({'NOAA ARs'; [NBu]}), ...
                      [], TCOLOR2, [MRGN 0.88], DFONT);
  % put on southern tags
  NOAAinx = NOAAv &  (NOAAlat(:) < 0); % valid and south
  NOAAtid1 = NOAAtid (NOAAinx); % track id's
  NOAAspt1 = NOAAspt (NOAAinx); % spot counts
  NOAAeq1  = NOAAequa(NOAAinx); % along-equator values
  [junk,NOAAup] = sort(NOAAeq1, 1, 'ascend'); % sort (valid&south) by equa
  NBd = num2str(NOAAtid1(NOAAup), '%05.0f');
  % NBs = num2str(NOAAspt1(NOAAup), ':%2d'); % remove these
  imT = text_on_image(imT, [NBd], [], TCOLOR2, [1-MRGN 0.88], DFONT);
else,
  imT = text_on_image(imT, strvcat({'NOAA ARs'; '(none)'}), ...
                      [], TCOLOR2, [MRGN 0.88], DFONT);
end;

frame = imT;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% columnize: break text block "t" into columns having at most nlin rows
%
function tp = columnize(t, nlin, up)

npad = 0; % #blanks between columns

[m1,n1] = size(t);
% fits on 1 column
if m1 <= nlin, 
  tp = t; 
  return; 
end;
% put into 2 columns
mf = 2*nlin - m1; % number of fill lines needed
% fail gracefully if too big for 2 columns
if mf < 0,
  mf = 0;
end;
fill = repmat(blanks(n1), [mf 1]); % fill block
pad = repmat(blanks(npad), [nlin 1]); % blank padding
mx = nlin - mf;
% incomplete column goes high or low
if up,
  col1 = strvcat(t(1:mx,:), fill);
else,
  col1 = strvcat(fill, t(1:mx,:));
end;

tp = [col1 pad t(mx+1:end,:)];
return;


