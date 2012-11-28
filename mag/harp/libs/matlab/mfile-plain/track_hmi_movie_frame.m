function frame=track_hmi_movie_frame(t, geom, mask, harps, bitmaps, smpl)
%track_hmi_movie_frame	make one frame of a tracker movie
% 
% frame=track_hmi_movie_frame(t, mask, harps, bitmaps, smpl)
% * Make one frame of a harp diagnostic movie.  The output is a single
% indexed image.  This uses code similar to track_hmi_movie, but is intended
% to retain only the frame-making parts, leaving the I/O to a separate
% routine.
% * Usage: 
%   = xxx('fulldisk-instant.cat', 
% 
% where 4 gives the spatial reduction (factor of 4 in each direction 
% yields 1024^2 images for HMI), cm is a colormap (try prism2), and 
% fnTs is a pattern for *track* summary mat-files.
% 
% Inputs:
%   real t;       -- datenum time
%   real mask(m,n)
%   struct harps
%   cell bitmaps(nr) of real
%   int smpl;
% 
% Outputs:
%   int frame(mp,np)
% 
% See Also: track_hmi_test

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 03 Oct 2012
% Copyright (c) 2009.  All rights reserved.  

% 
% Error checking
% 
% if all(nargin  ~= [5]), error ('Bad input arg number'); end;
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

rgnAge = (t - harps.t_frst1); % 1 if new
rgnPad = (t > harps.t_last1) - (t < harps.t_frst1); % -1/0/+1
rgnMrg = (harps.h_merge >  0); % 1 if merged
rgnTag = double(harps.h_faint == 0); % 1 = normal, 0 = faint...
rgnTag((harps.h_faint == 1) & isnan(harps.npix)) = -1; % ...-1 = placeholder
rgnTid = (harps.harpnum)';

frame_num = NaN; % unknown
prior_t = NaN; % unknown

mask = mask(1:smpl:end,1:smpl:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NOAA Regions

% load NOAA regions
NOAAar = hmi_noaa_info_interp(t);
% project into our disk geometry
[NOAAx, NOAAy, NOAAz] = hmi_latlon2image([NOAAar.latitudehg ]', ...
                                         [NOAAar.longitudecm]', geom);
NOAAx = round(NOAAx/abs(smpl));
NOAAy = round(NOAAy/abs(smpl));
NOAAz = round(NOAAz/abs(smpl));
NOAAv = (NOAAz >= 0); % on-disk filter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Colors, margins

% set FLIP
% +1 for no flip (P=0), -1 for flip (P=180)
p0 = geom(5); 
FLIP = 2*(cosd(p0) > 0) - 1;

% margins
MRGN = 0.01; % [0,1] units
MRGNPIX = floor(MRGN*size(mask, 1)); % pixels

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

% inputs to r2t remapper come from znan(...) construction, encoded as:
%    offdisk starts as NaN, mapped to 1
%    ondisk starts as 0, mapped to 2
%    blobs start as >=1, mapped to >= FRSTROI
% outputs of r2t mapper: 
%    offdisk => LOLITE, ondisk => NEUTA, blobs => FRSTROI + 0, 1, ...

% keep the initialization simple
im1 = zeros(size(mask));
im1(mask == 0) = LOLITE; % off-disk
im1(mask >  0) = NEUTA;  % on-disk
for inx = 1:nR,
    % value for on-blob region
    imval = FRSTROI + rem(rgnTid(inx)-1, Nt_show);
    % indicator for on-blob region
    new = bitmaps{inx}(1:smpl:end,1:smpl:end) > 32; % plug active later
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make a full-disk image with segmentation + boundary info superimposed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1: (we do not have have both posterior and prior boundary)

% 2: active are LOLITE
im1(mask == 2) = LOLITE;

% 3: box around each region
% Note: region corners are in the ims() coordinates, unflipped/untransposed
rgnS = floor((rgn-1) / smpl) + 1;
if smpl == 1,
  % L-R and T-B hairlines offset 0 and 1
  tk1 = [1 1 3 3]; % four hairlines...
  tk2 = tk1 + 1;
  del = [0 1 0 1]; % ...offset by 0, 1, 0, 1
  flag = 6; % size of merged-region panel
else,
  % hairlines not offset
  tk1 = [1 3]; % two hairlines
  tk2 = tk1 + 1;
  del = [0 0]; % not offset
  flag = 3;
end;
[M,N] = size(im1);
for r=1:nR,
  if rgnTag(r) < 0, continue; end; % no enclosing box for placeholder regions
  % throughout, have to clip to range of indexes to prevent trouble
  im1(rgnS(r,1):rgnS(r,3), range(rgnS(r,tk2) + del, 1, N)) = HILITE;
  im1(range(rgnS(r,tk1) + del, 1, M), rgnS(r,2):rgnS(r,4)) = HILITE;
  if rgnTag(r) < 0,
    % NB: this case is now skipped directly (above)
    % dotted box for a placeholder region
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

% put on NOAA AR markers 
im1 = marker_on_image(im1, [NOAAx(NOAAv) NOAAy(NOAAv)], ...
                      '+', [15 3], [], CONTRA);
im1 = text_on_image(im1, cellstr(int2str([NOAAar(NOAAv).regionnumber]')), [], CONTRA, ...
                    [NOAAx(NOAAv)-15 ...
                    2*sign(NOAAy(NOAAv)-M/2).*sqrt(abs(NOAAy(NOAAv)-M/2))+M/2 ...
                    zeros(nnz(NOAAv),1)-FLIP], SFONT);

% put on the poles, if they're there
[NPx, NPy, NPz] = hmi_latlon2image([90;-90], [0;0], geom);
NPx = round(NPx/abs(smpl));
NPy = round(NPy/abs(smpl));
NPv = (NPz >= 0); % visible on-disk
NP_txt = {'N'; 'S'};
NPspc = [-20;10]; % offset text from marker
im1 = marker_on_image(im1, [NPx(NPv) NPy(NPv)], 'x', [15 3], [], HILITE);
im1 = text_on_image(im1, NP_txt(NPv), [], HILITE, ...
                    [NPx(NPv) NPy(NPv)+NPspc(NPv) zeros(nnz(NPv),1)-FLIP], SFONT);

% HARP-number labels -- only the non-placeholder ones
if nR > 0,
  harp_lbl_locn = [rgnS(:,1:2)+1 zeros(nR,1)-FLIP];
  im1 = text_on_image(im1, cellstr(int2str(rgnTid(rgnTag >= 0))), ...
                      [], HILITE, harp_lbl_locn(rgnTag >= 0,:), SFONT);
end;

% keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make the text overlays
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We have to rotate the image into final form as-presented (North up)
% before we put the text labels on.
% By contrast, the surrounding boxes have to be done in original
% Matlab coordinates, since they were originally located that way.

% rotate the image up
% FIXME: not sure about this!
im1 = rot90(im1, FLIP);
M = size(im1, 1); % scale info

% make space in im1 for overlay
if     M < 400, pad = [32 96]; 
elseif M < 800, pad = [ 0 96]; 
else            pad = [ 0 96];
end;
% imT is the image-with-text
imT = [zeros(M, pad(1))+LOLITE im1 zeros(M, pad(2))+LOLITE];

% time overlay
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
txt = {{'SDO/HMI Tracked AR (HARP)', t1code, t2code, ' ', t3code, ' ', fnumT}};
imT = text_on_image(imT, txt, [], TCOLOR, [MRGN MRGN], DFONT);

% nR is the number of tracks, accounting for merges
txt = {{...
    sprintf('T = %d', nR), ...
    ' ', ...
    sprintf('! = %d (new)',         nnz(rgnAge == 0)), ...
    sprintf('( = %d (pad before)',  nnz(rgnPad <  0)), ...
    sprintf(') = %d (pad after)',   nnz(rgnPad >  0)), ...
    sprintf('+ = %d (merge)',       nnz(rgnMrg >  0)), ...
    sprintf('~ = %d (use past)',    nnz(rgnTag == 0)), ...
    sprintf('? = %d (placeholder)', nnz(rgnTag <  0)), ...
      }};
% put overlay on
imT = text_on_image(imT, txt, [], TCOLOR, [1-MRGN MRGN], DFONT);

% track number overlay
NC_max = floor(size(imT, 1) / (2*19)); % max track ID's in one column
if nR > 0,
  % compile sizes
  rgnSize = (rgn(:,3)-rgn(:,1)) .* (rgn(:,4) - rgn(:,2)); % NB: orig. units
  [junk,botSize] = sort(rgnSize);
  botSize(rgnMrg > 0) = []; % don't double-indicate
  % compile locations
  UPPER = (rgnS(:,2)+rgnS(:,4)) > M; % center in Northern hemi
  if FLIP < 0, UPPER = ~UPPER; end;
  % now label regions
  rgnID_t = num2str(rgnTid, '%04.0f'); % text format
  rgnID_x = blanks(nR)';
  if 0 == 1 && length(botSize) > 0,
    % disabled
    rgnID_x(botSize(1))   = 'v'; % smallest
    rgnID_x(botSize(end)) = '^'; % biggest
  end;
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
  imT = text_on_image(imT, strvcat({'HARPs'; TBu}), [], TCOLOR, [MRGN   0.97], DFONT);
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
end;

% put on NOAA track numbers
if nnz(NOAAv) > 0,
  NOAAlat = [NOAAar.latitudehg]; % all lat's
  NOAAtid = [NOAAar.regionnumber]';
  NOAAspt = [NOAAar.spotcount]';
  % put on northern tags
  NOAAinx = NOAAv & (NOAAlat(:) >= 0); % valid and north
  NOAAtid1 = NOAAtid(NOAAinx); % track id's
  NOAAspt1 = NOAAspt(NOAAinx); % spot counts
  NOAAx1   = NOAAx  (NOAAinx); % x-values
  [junk,NOAAxup] = sort(NOAAx1, 1, 'descend'); % sort valid + north by x
  NBu = num2str(NOAAtid1(NOAAxup), '%05.0f');
  NBs = num2str(NOAAspt1(NOAAxup), ':%2d');
  imT = text_on_image(imT, strvcat({'NOAA ARs'; [NBu NBs]}), ...
                      [], TCOLOR2, [MRGN 0.88], DFONT);
  % put on southern tags
  NOAAinx = NOAAv & (NOAAlat(:) < 0); % valid and south
  NOAAtid1 = NOAAtid(NOAAinx); % track id's
  NOAAspt1 = NOAAspt(NOAAinx); % spot counts
  NOAAx1   = NOAAx  (NOAAinx); % x-values
  [junk,NOAAxup] = sort(NOAAx1, 1, 'descend'); % sort valid + south by x
  NBd = num2str(NOAAtid1(NOAAxup), '%05.0f'); % text format
  NBs = num2str(NOAAspt1(NOAAxup), ':%2d');
  imT = text_on_image(imT, [NBd NBs], [], TCOLOR2, [1-MRGN 0.88], DFONT);
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


