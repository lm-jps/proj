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

if exist(fontfile, 'file'),
  fontpath = which(fontfile);
else,
  error('Could not find fontfile %s in path', fontfile);
end;
% the basis string for the font name, and the extension
[junk1,fontbase,font_ext] = fileparts(fontfile);
clear junk1 junk2


% get the bitmaps and char widths
if strcmp(font_ext, '.ttf'),
  [bm,cpos] = text_on_image_font_helper(fontpath, dict, height, width);
  bm = permute(bm, [2 1 3]); % flip x and y
elseif strcmp(font_ext, '.m'),
  info = feval(fontbase); % loads into "info"
  if ~isstruct(info),
    error('M-file font must load into info structure');
  end
  % TODO: probably want the font's own height here, from info.x
  [bm,cpos] = text_on_image_font_mfile(info, dict, height, width);
else,
  error('Unrecognized font type.');
end;

% dump it into a structure
font.bitmaps = bm;
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
% return bitmaps (bm) and character positions (cpos)
%
% We have to allow for padding, which causes a lot of special cases
% 

function [bm,cpos] = text_on_image_font_mfile(info, dict, height, width)

    %% defaults and constants
    Nmap = length(info.bitmaps);
    Nchar = length(dict);
    % height defaults
    if length(height) < 1, height(1) = 32; end; % unused for bitmaps
    if length(height) < 2, height(2) = 1; end;
    if length(height) < 3, height(3) = 0; end;
    if length(height) < 4, height(4) = 0; end;
    % width defaults
    if length(width)  < 1, width(1) = double('0'); end;
    if length(width)  < 2, width(2) = 0; end;
    if length(width)  < 3, width(3) = 0; end;
    % set top/bottom padding
    if height(2) == 0,
        padB = (info.ascent - info.height) + height(4);
    else
        padB = height(4);
    end;
    padT = height(3);
    % set left/right padding
    padL = width(2); 
    padR = width(3);
    %% load bitmaps
    % set character size -- default uses size of '0'
    [Isiz1,Isiz2] = size(info.bitmaps{width(1)});
    Osiz1 = Isiz1 + padB + padT;
    Osiz2 = Isiz2 + padL + padR;
    % bitmaps wil be truncated this much
    min_bm_width = max(0, -padL) + max(0, -padR);
    % resize 0x0 bitmaps to this height
    canonical_bm_height = Isiz1;
    %% Loop to fill in bm and cpos
    bm = zeros(Osiz1, Osiz2, Nchar, 'uint8');
    cpos = zeros(2, Nchar);
    for i = 1:Nchar,
        inx = double(dict(i))+1;
        if inx > Nmap, continue; end; % should never happen
        % mfile is in raster order
        bm1 = uint8(255*info.bitmaps{inx});
        % ensure correct height
        if isempty(bm1),
          bm1 = zeros(canonical_bm_height, 0, 'uint8');
        end;
        [siz1,siz2] = size(bm1);
        % adjust narrow chars so they may be truncated by padL/R
        if siz2 < min_bm_width,
          bm1 = [bm1 zeros(canonical_bm_height, min_bm_width-siz2, 'uint8')];
        end;
        %% pad on four sides -- wordy but functional
        % right pad (do before left)
        if padR < 0, 
            bm1(:,end+padR+1:end) = []; % remove
        elseif padR > 0,
            bm1(:,end+1:end+padR) = 0; % add
        end;
        % left pad
        if padL < 0, 
            bm1(:,1:-padL) = []; % remove
        elseif padL > 0,
            bm1 = [zeros(size(bm1,1),padL) bm1]; % add columns
        end;
        % bottom pad (do before top)
        if padB < 0, 
            bm1(end+padB+1:end,:) = []; % remove
        elseif padB > 0,
            bm1(end+1:end+padB,:) = 0; % add
        end;
        % top pad
        if padT < 0, 
            bm1(1:-padT,:) = []; % remove
        elseif padT > 0,
            bm1 = [zeros(padT,size(bm1,2)); bm1]; % add rows
        end;
        % force into given size if necessary
        bmplug = bm1;
        bm_width = size(bm1, 2);
        if bm_width < Osiz2,
          bmplug(:,end+1:Osiz2) = 0; % augment
        elseif bm_width > Osiz2,
          bmplug(:,Osiz2+1:end) = []; % truncate
        end;
        % plug it in
        bm(:,:,i) = bmplug;
        % these char widths seem to be unused
        cpos(1,i) = 1+padL;
        cpos(2,i) = bm_width;
    end;
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
