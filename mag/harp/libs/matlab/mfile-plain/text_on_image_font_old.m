function font=text_on_image_font_old(fontfile)
%text_on_image_font_old	load font (old version) for text_on_image
% 
% font=text_on_image_font_old(fontfile)
% * Load a stored font `fontfile' containing letters from a given dictionary
% rendered at a given size.  Return a structure encapsulating the font.
% * The fontfile is transformed by a filename convention into a `load'able 
% file.  Supported fontfiles are dejavu{18,24,36} (Deja Vu Monospace).  
% * If supplied as empty, fontfile is taken as `dejavu24'.
% * The character bitmaps are stored as a 3d array containing one entry 
% for ascii SPC thru ascii ~.  The characters are thus all the same size.  
% * These bitmaps were found by screen grabs of the fonts in TextEdit, and 
% chopping a long string into individual characters.  Their bitmaps are read 
% from a mat-file.
% * To skip repeated file reads, the output structure is cached across 
% function invocations (the "persistent" mechanism).  To clear the cache
% and return the cached structure, supply no input arguments.
% 
% Inputs:
%  string fontfile = 'dejavu24'
% 
% Outputs:
%  struct font
% 
% See Also:  text_on_image; ./Related/text_on_image for bitmap information
% text_on_image_font for other font construction

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 25 Sep 09.
% Copyright (c) 2009.  All rights reserved.

% This is the string we used, taken from `man ascii'.
% ' !"#$%&''()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

persistent text_on_image_bitmaps;
if isempty(text_on_image_bitmaps),
  text_on_image_bitmaps = struct([]);
end;
% hook to return the bitmaps and/or clear the cache
if nargin == 0, 
  font = text_on_image_bitmaps;
  clear text_on_image_bitmaps;
  return;
end;

% Error checking
% 
% if all(nargin  ~= [1 2 3]), error ('Bad input arg number'); end
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end  
if isempty(fontfile), fontfile = 'dejavu24'; end;

%
% Computation
% 

% check the cache
if isfield(text_on_image_bitmaps, fontfile),
  % it is cached; return it
  font = text_on_image_bitmaps(1).(fontfile);
  return;
end;

% otherwise, we must load the font's bitmaps
fn = sprintf('text_on_image_%s.mat', fontfile);
if ~exist(fn, 'file'),
  error(sprintf('Could not find font <%s> from file <%s>', fontfile, fn));
end;
S = load(fn);
bm = S.bitmaps;

% this is the only dictionary we use for these fonts
dict = [' !"#$%&''()*+,-./0123456789:;<=>?@' ...
        'ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'];

% dump the information into a structure
font.bitmaps = bm;
font.dict = dict;
font.fontfile = fn;
font.name = sprintf('%s_%s', fontfile(1:end-2), fontfile(end-1:end)); % crude
font.pos = zeros(1, size(bm, 3));

% cache it
text_on_image_bitmaps(1).(fontfile) = font;

return;
end

