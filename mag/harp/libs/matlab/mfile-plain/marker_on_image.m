function out = marker_on_image(img, pos, m, mprop, icolor, mcolor)
%marker_on_image	overlay markers onto image
% 
% out = marker_on_image(img, pos, m, mprop, icolor, mcolor)
% * Overlay marker shape `m' at positions `pos' within image `img'.  
% * To give pixel coordinates, for pos, use values according to size(img).
% For scaled coordinates, use (1,1) for lower right, (0,0) for upper
% left, and (0.5,0.5) for centered.  We switch automatically between
% the two coordinate types by looking at max(pos(:))
% * Properties of the marker are given in `mprop'.  mprop(1) is an 
% integer size in pixels, and mprop(2) is an integer stroke thickness
% in pixels.
% * Both uint8 truecolor images and indexed images are supported, with 
% somewhat different icolor/mcolor interpretations.  
% Colors are a bit rococo.  
% * The image color is controlled by icolor; the marker color by mcolor.
% For truecolor images, scalar colors are expanded to RGB triples of the 
% given color, and everything must be between 0 and 1.
% For indexed images, icolor and mcolor are the (integer) colormap
% indexes.  (Real values like NaN and Inf are "correctly" treated
% by Matlab, also.)  
% * Additionally, for scalar images, vector mcolor is OK.  It is
% interpreted as defining a mapping from the overlay range [0,1] to
% a color index mcolor.  With smoother marker shapes, this can give 
% smoother outlines, even if mcolor is only of length 3.  But, our 
% markers are currently 0/1 so only scalar mcolor makes sense.
% * For truecolor images, the color of `out' is the *sum* of the 
% background (image), and the foreground (marker bitmaps), where:
%   bg = cat(3, icolor(1)*img(R), icolor(2)*img(G), icolor(3)*img(B))
%   fg = cat(3, mcolor(1)*bitmap, mcolor(2)*bitmap, mcolor(3)*bitmap)
% For blue markers, use mcolor = [0 0 1].  For a black background, use 
% icolor = 0.  To keep the background intact, use icolor = 1.  Off-scale
% colors (outside of [0,255] after linear combination) are clipped to
% stay in-range.
% * For indexed images, and simple scalar icolor/mcolor, the color
% index of `out' is icolor if the overlay is 1, or mcolor otherwise.
% For vector mcolor, we map the overlay range thru mcolor.  In all
% cases, if the overlay is "0", icolor is placed in the output.  If
% icolor is empty, *the output receives the image value* where the
% overlay equals 0.  This is the default for indexed images.
% * The marker is output from marker_on_image_marker. See that routine 
% for more on marker shapes.
% 
% Inputs:
%   uint8 img(m,n,3) or double img(m,n)
%   real pos(np,2)
%   string m          -- '+', 'x', '*', etc.
%   opt int mprop(2) = [15 1]
%   opt real icolor(3) or (1)
%      (truecolor: icolor = [0 0 0]. indexed: icolor = [].)
%   opt real mcolor(3) or (1) or (nc)
%      (truecolor: mcolor = [0.7 0.7 0.7]. indexed: icolor = 1.)
% 
% Outputs:
%   uint8 out(m,n,3) or double out(m,n)
%    
% See Also: marker_on_image_marker

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 25 Sep 11.
% Copyright (c) 2009.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [3:6]), error ('Bad input arg number'); end
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end  

% truecolor or indexed?
P = size(img, 3); % 1 or 3
% set up args
if nargin < 6, 
  % (marker) truecolor: gray. indexed: "color #2"
  if P == 3, mcolor = 0.7; else, mcolor = 2; end;
end;
if nargin < 5, 
  % (image) truecolor: black.  indexed: "color #1"
  if P == 3, icolor = 0.0; else, icolor = 1; end;
end;
if nargin < 4, mprop = [15 1]; end;
% expand scalar colors to length-P
if P == 3,
  if length(icolor) == 1,
    icolor = icolor + zeros(1,P);
  end;
  if length(mcolor) == 1,
    mcolor = mcolor + zeros(1,P);
  end;
end;
% trivial case
if isempty(pos), out = img; return; end;
% some more checking
if size(pos,2) ~= 2, error('pos must have (x,y) tuples'); end;
if ~ischar(m), error('Marker must be a string'); end;
if length(mprop) ~= 2, error('mprop has two entries'); end;
if length(mcolor) ~= P && P == 3, error('Bad mcolor'); end; %P=1:vector OK
if length(icolor) ~= P && ~isempty(icolor), error('Bad icolor'); end;
if all(P ~= [1 3]), error('Need empty, indexed, or RGB image'); end;

%
% Computation
% 
Nmark = size(pos, 1);
[M,N,P] = size(img);
mscale = mprop(1);
mwidth = mprop(2);
mhwid  = (mscale - 1)/2;  % half-width of scale X scale box

% get scaling on "pos"
if any(pos(:) > 1)
  % pixel coords given
  scalex = 1; 
  scaley = 1;
  offset = 0; 
else,
  scalex = N;
  scaley = M;
  offset = 1; % off-by-one matlab coordinates
end;

% get marker stencil -- currently just 0/1
[mx, my] = marker_on_image_marker(m, mscale, mwidth);

% translate markers to image places
m_inx = [];
for i = 1:Nmark,
  % don't forget to move the origin down by half a marker-width
  mx1 = round(mx + scalex * pos(i,1) + offset - mhwid);
  my1 = round(my + scaley * pos(i,2) + offset - mhwid);
  % crop to be on-image
  offscale = (mx1 < 1) | (mx1 > N) | (my1 < 1) | (my1 > M);
  mx1(offscale) = [];
  my1(offscale) = [];
  % accumulate index sets
  inx1 = sub2ind([M N], mx1, my1);
  m_inx = cat(2, m_inx, inx1);
end;
block = zeros(M,N);
block(m_inx) = 255;

% keyboard

% set up out1, the text block in the image
if P == 3,
  % blend icolor with mcolor
  out = uint8(double(img)                    .* repmat(reshape(icolor, [1 1 3]), [M N 1]) + ...
              double(repmat(block, [1 1 3])) .* repmat(reshape(mcolor, [1 1 3]), [M N 1]));
else,
  % set up output colors
  if length(mcolor) == 1,
    % marker in two colors: icolor and mcolor
    if isempty(icolor), 
      outcolor = [mcolor mcolor]; % just have 2 points
    else,
      outcolor = [icolor mcolor];
    end;
  else,
    % text in (vector) mcolor
    outcolor = mcolor;
  end;
  % map range of the entries in block(:) thru outcolor
  inx = linspace(0, 255, length(outcolor));
  out = reshape(interp1(inx, outcolor, double(block(:)), 'nearest'), ...
                [M N]);
  % if icolor was given, ensure output=icolor where the bitmap == 0
  if ~isempty(icolor),
    out(block == 0) = icolor;
  else,
    out(block == 0) = img(block == 0);
  end;
end;
return;
