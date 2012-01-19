function font=text_on_image_font(fontfile, dict, height, width)
%text_on_image_font	load font for text_on_image
% 
% font=text_on_image_font(fontfile, dict, height, width)
% * Render letters from dict, using a given TrueType fontfile, into
% bitmaps and other information that is packed into a structure 
% `font' which can be handed to text_on_image.
% * The default dict is a large set of printable ascii characters.
% If dict is given as empty, this default is chosen as well.
% * See text_on_image_font_helper for full interpretation of height 
% and width.  Height is the more important, it is in the form:
%  height = [char_height allow_space_for_descenders top_pad bottom_pad]
% where char_height is in pixels, the next is 0/1, and the 
% top and bottom pad are in pixels and can be negative.
% * To skip repeated file reads, the output structure is cached across 
% function invocations (the "persistent" mechanism).  To clear the cache
% and return the cached structure, supply no input arguments.
% 
% Inputs:
%  string font
%  opt string dict(n) = ' !"#$%&''()*+,-./0123456789:;<=>?@' + ...
%        'ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'
%  opt int height(0, 1, 2, 3, or 4) = [32 1 0 0]
%  opt int width(0, 1, 2, or 3) = [double('0') 0 0]
% 
% Outputs:
%  struct font
% 
% See Also:  text_on_image

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

% 
% Error checking
% 
% if all(nargin  ~= [1 2 3]), error ('Bad input arg number'); end
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end  
if nargin < 2 || isempty(dict), 
  dict = [' !"#$%&''()*+,-./0123456789:;<=>?@' ...
          'ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'];
end;
if nargin < 3, height = [32]; end;
if nargin < 4, width  = [];   end;

% check the cache
[inx,sign1] = args_in_cache(text_on_image_bitmaps, fontfile, dict, height, width);
if inx > 0,
  % it is cached; return it
  font = text_on_image_bitmaps(inx).font;
  return;
end;

% not in cache now: load, render, and cache

% get the bitmaps and char widths
if exist(fontfile, 'file'),
  fontpath = which(fontfile);
else,
  error('Could not find fontfile %s in path', fontfile);
end;
[bm,cpos] = text_on_image_font_helper(fontpath, dict, height, width);

% just the basis string for the font name
[junk1,fontbase,junk2] = fileparts(fontfile);
clear junk1 junk2

% dump it into a structure
font.bitmaps = permute(bm, [2 1 3]); % flip x and y
font.dict = dict;
% font.sign = sign1; % the signature
font.fontfile = fontfile;
font.name = sprintf('%s_%d', fontbase, height(1));
font.pos = cpos;

% cache it
text_on_image_bitmaps(end+1).sign = sign1;
text_on_image_bitmaps(end  ).font = font;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% return index number of a cache match, or -1

function [m,sign1] = args_in_cache(info, fontfile, dict, height, width)

sign1 = args2sign(fontfile, dict, height, width);
for k = 1:length(info),
  if all(size(info(k).sign) == size(sign1)),
    if all(info(k).sign == sign1),
      m = k; % match
      return;
    end;
  end;
end;
% not found
m = -1;
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% convert argument list to signature vector

function sign = args2sign(fontfile, dict, height, width)

% something illegal in any argument, to eliminate accidental alignments
% of the arguments, which are all variable length
% (all negative fractions are illegal)
GAP = [-1e10; -1e-20];
% concatenate the arguments
sign = [double(fontfile(:)); GAP; double(dict(:)); GAP; height(:); GAP; width(:)];
return;
end
