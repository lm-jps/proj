function out = text_on_image(img, txts, icolor, tcolor, pos, font)
%text_on_image	overlay text onto image
% 
% out = text_on_image(img, txts, icolor, tcolor, pos, font)
% * Overlay text (txts) on image (img).  If img is empty, just the
% overlay is returned.  Each overlay pixel is a scalar value,
% so it can serve as an index or a weight on an RGB color.
% * Both uint8 truecolor images and indexed images are supported, with 
% somewhat different icolor/tcolor interpretations.
% * Non-printable ascii characters are turned into spaces.  Character 
% matrices and cell arrays of strings are stacked.
% * The image color is controlled by icolor; the text color by tcolor.
% For truecolor images, scalar colors are expanded to RGB triples of the 
% given color, and everything must be between 0 and 1.
% For indexed images, icolor and tcolor are the (integer) colormap
% indexes.  (Real values like NaN and Inf are "correctly" treated
% by Matlab, also.)  
% * Additionally, for scalar images, vector tcolor is OK.  It is
% interpreted as defining a mapping from the overlay range [0,1] to
% a color index tcolor.  This can give smoother outlines, even if 
% tcolor is only of length 3.
% * For truecolor images, the color of `out' is the *sum* of the 
% background (image), and the foreground (letter-sequence bitmaps), where:
%   bg = cat(3, icolor(1)*img(R), icolor(2)*img(G), icolor(3)*img(B))
%   fg = cat(3, tcolor(1)*bitmap, tcolor(2)*bitmap, tcolor(3)*bitmap)
% For blue text, use tcolor = [0 0 1].  For a black background, use 
% icolor = 0.  To keep the background intact, use icolor = 1.  Off-scale
% colors (outside of [0,255] after linear combination) are clipped to
% stay in-range.
% * For indexed images, and simple scalar icolor/tcolor, the color
% index of `out' is icolor if the overlay is < 128, or tcolor otherwise.
% For vector tcolor, we map the overlay range thru tcolor.  In all
% cases, if the overlay is "0", icolor is placed in the output.  If
% icolor is empty, *the output receives the image value* where the
% overlay equals 0.  This is the default for indexed images.
% * The position of the letter bitmaps within the image is given by 
% pos(1:2,:).  Scaled coordinates (i.e., (1,1) for lower right, 
% (0,0) for upper left, and (0.5,0.5) for centered) are assumed, unless
% some position exceeds 1, in which case pixel coordinates are assumed.
% * The third component of pos, if given, encodes the rotation.  The 
% most useful value is 1, which gives standard vertical text.
% * The font is a font structure which is output from 
% text_on_image_font, or text_on_image_font_old.  See those routines 
% for more.
% 
% Inputs:
%   uint8 img(m,n,3) or double img(m,n)
%   char txts(mt,nt) or cell txt{mt} of char()
%   opt real icolor(3) or (1)
%      (truecolor: icolor = [0 0 0]. indexed: icolor = [].)
%   opt real tcolor(3) or (1) or (nc)
%      (truecolor: tcolor = [0.7 0.7 0.7]. indexed: icolor = 1.)
%   opt real pos(3) = [1 1 0]  -- in [0,1] x [0,1] x {-2,-1,0,1,2}
%   opt struct font = text_on_image_font_old('dejavu24');
% 
% Outputs:
%   uint8 out(m,n,3) or double out(m,n)
%    
% See Also: text_on_image_font for font construction

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 25 Sep 09.
% Copyright (c) 2009.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [2 3 4 5 6]), error ('Bad input arg number'); end
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end  

% truecolor or indexed?
P = size(img, 3); % 1 or 3
% set up args
if nargin < 6, 
  font = text_on_image_font_old('dejavu24');
end;
if nargin < 5 || isempty(pos), pos = [1 1]; end;
if size(pos, 2) == 2, 
  pos = [pos zeros(size(pos,1),1)]; % plug in final 0
end; 
if nargin < 4, 
  % (text) truecolor: gray. indexed: "color #2"
  if P == 3, tcolor = 0.7; else, tcolor = 2; end;
end;
if nargin < 3, 
  % (image) truecolor: black.  indexed: "color #1"
  if P == 3, icolor = 0.0; else, icolor = []; end;
end;
% convert char txts to one-entry cell
if ischar(txts), 
  txt1 = {txts}; txts = txt1; clear txt1; 
end;
% convert pos to have one entry per txt
if length(txts) > 1 && size(pos,1) == 1,
  pos = repmat(pos, length(txts), 1); % preserves uint/double status
end;
% expand scalar colors to length-P
if P == 3,
  if length(icolor) == 1,
    icolor = icolor + zeros(1,P);
  end;
  if length(tcolor) == 1,
    tcolor = tcolor + zeros(1,P);
  end;
end;
% some more checking
if ~isstruct(font), error('Font must be a struct'); end;
if length(txts) >= 1 && size(pos, 1) ~= length(txts), 
  error('pos must be a row, or have one row per txt'); 
end;
if length(tcolor) ~= P && P == 3, error('Bad tcolor'); end; %P=1:vector OK
if length(icolor) ~= P && ~isempty(icolor), error('Bad icolor'); end;
if ~iscell(txts), error('txts must be a string matrix, or a cell'); end;
if all(P ~= [1 3]), error('Need indexed or RGB image'); end;
% if all coords are <= 1, assume they are scaled (allow small negatives too)
scaled_coords = all(all(abs(pos(:,1:2)) <= 1));

%
% Computation
% 

[M,N,P] = size(img);

% text alignment: implemented, but no interface
align1 = 1;
align2 = 1;

%% Compute dictionary info
bitmaps = font.bitmaps;
dict = font.dict;
% bitmap letters are Mb by Nb
[Mb,Nb,junk] = size(bitmaps);
% find index within dict of the char we will call blank
blankchar = find(dict == ' ', 1); % first space
if isempty(blankchar), blankchar = 1; end;

%% Loop over text blocks, placing them in the overlay
overlay = zeros(M, N, 'uint8');
Ntxt = length(txts);
for i = 1:Ntxt,
  % char() converts a cell array of strings to a string matrix
  txt = char(txts{i});
  [Mt,Nt] = size(txt);

  %% Find dictionary index for each char in txt
  % distance from txt to dict (Ntxt x Ndict)
  dist = abs(repmat(dict, [Mt*Nt 1]) - repmat(txt(:), [1 length(dict)]));
  [exact,inx] = min(dist, [], 2); % min along dict axis
  % create indexes into dict
  indexes = zeros(size(txt)) + blankchar; % all blank for now
  indexes(exact == 0) = inx(exact == 0); % only replace exact matches

  % bitmaps for each letter, stacked along third dim
  letters = bitmaps(:,:,indexes(:));
  % reshape so a multi-line txt can be permuted
  block = reshape(letters, [Mb Nb Mt Nt]);
  % permute and reshape into the text as displayed; rotate as desired
  block = rot90(reshape(permute(block, [1 3 2 4]), [Mb*Mt Nb*Nt]), double(pos(i,3)));

  % determine where to place the block
  [m,n] = size(block);
  % FIXME: are we handling fenceposts correctly?
  if scaled_coords,
    origin = [(M-m) (N-n)] .* pos(i,1:2);
  else
    origin = pos(i,1:2) - 1;
  end;
  if align1 == 0,
    offset1 = m/2; % centered
  elseif align1 == -1,
    offset1 = -m;
  else,
    offset1 = 0;
  end;
  if align2 == 0,
    offset2 = n/2; % centered
  elseif align2 == -1,
    offset2 = -n;
  else,
    offset2 = 0;
  end;
  inx1 = round([1:m] + origin(1) + offset1);
  inx2 = round([1:n] + origin(2) + offset2);
  % allow for clipping of the block
  inrange1 = (inx1 > 0) & (inx1 <= M);
  inrange2 = (inx2 > 0) & (inx2 <= N);
  % place block in overlay
  overlay(inx1(inrange1),inx2(inrange2)) = block(inrange1,inrange2);
end;

%% Set up the output
if P == 3,
  % truecolor case
  % blend icolor with tcolor
  out = uint8(bsxfun(@times, double(img),     reshape(icolor, [1 1 3])) + ...
              bsxfun(@times, double(overlay), reshape(tcolor, [1 1 3])));
else,
  % indexed case
  % set up output color map
  if length(tcolor) == 1,
    % text in two colors: icolor and tcolor
    if isempty(icolor), 
      outcolor = [tcolor tcolor]; % just have 2 points
    else,
      outcolor = [icolor tcolor];
    end;
  else,
    % text in (vector) tcolor
    outcolor = tcolor;
  end;
  % map range of the entries in overlay(:) thru outcolor
  inx = linspace(0, 255, length(outcolor));
  out = reshape(interp1(inx, outcolor, double(overlay(:)), 'nearest'), ...
                size(overlay));
  % ensure output is unchanged where bitmap == 0
  if ~isempty(icolor),
    % case where icolor was given explicitly
    out(overlay == 0) = icolor;
  else,
    out(overlay == 0) = img(overlay == 0);
  end;
end;
return;
