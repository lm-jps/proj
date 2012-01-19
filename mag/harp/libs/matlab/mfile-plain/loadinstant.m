function [im,rgn,rgnC,rgnT,cats,rgnR,fns]=loadinstant(fn, keys, skip)
%loadinstant	load components of a StarTool instant file
% 
% [im,rgn,rgnC,rgnT,cats,rgnR,fns]=loadinstant(fn, keys, skip)
% * loads the named components of an instant file.  Loosely, such a file
% contains filenames of image files that are associated with each other,
% together with region-of-interest information, in a 
% 	key<whitespace>filename
% format (# signals comments).
% * If only certain keys are wanted, this can be signalled by
% passing them in via the keys cell array.  If keys is numeric, all
% images are read.
% * The integer skip gives the spatial subsampling.  Its absolute
% value is given as the second argument to loadfits/loadfitsa.
% This is as in file_loop_lst.
% * If skip is negative, "loadfitsa" is used instead of the default
% "loadfits".  The latter transposes while loading.
% * The "region" key is special.  When given, it causes the
% corresponding file to be opened, and all region information within
% is read into the rgn,rgnC,rgnT,rgnR outputs.  Old-format ".rgn" files 
% contain entries like:
%   # 1996/05/19 19:08:35
%   # 2/5+4/5@(241:288,465:512)
%   +box(40.5,40.5,48.0,48.0,0) # yellow
% The first comment line is a date and time.  The second gives the
% region number (or numbers, if more than one, separated by a + sign)
% within the file at that time, the total number of
% regions in the file, and the coordinates of the region within the
% original image (not the mosaic image).  The "box" command
% gives the region location in the mosaic image in saoTNG coordinates.
% * New-format ".reg" files contain entries like:
% box(930.0,544.5,97.0,96.0,0) # color=magenta 
%   text={1996/11/17 23:58:35; 2/2@(937:1033,17:112)}
% (without the line break present above).  The elements are just formatted
% slightly differently.
% * Images are returned in im(), region boundaries in rgn and rgnC, 
% associated times in rgnT, and the cat-files for accessing the 
% original images in cats.  The latter are extracted from the region
% file header.  The region number info (2/5+4/5) is returned in rgnR.
% * To extract the kth region from im, use
% >> im(rgn(k,1):rgn(k,3),rgn(k,2):rgn(k,4),nI)
% * To draw boxes around the regions on a plotted image, use
% >> line(rgn(:,[2 4 4 2 2])', rgn(:,[1 1 3 3 1])')
% * Note, if 'region' or 'reg' is not within keys, the region info
% will not be read.
% 
% Inputs:
%  string fn;  -- a .instant filename
%  opt string keys(nK) = 1;  -- read all images
%  opt int skip = 1;  -- loadfits
% 
% Outputs:
%  real im(m,n,nD);
%  real rgn(nR,4);
%  real rgnC(nR,4);
%  real rgnT(nR);
%  string cats(nI);
%  string rgnR(nR);
%  
% 
% See Also:  

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 04 Mar 99.
% Copyright (c) 1999.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [1 2 3]), error ('Bad input arg number'); end;
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end;
if (nargin < 2), keys = 1; end; % anything numeric
if isempty(keys), keys = 'NO_SUCH_KEY'; end; % nothing will be read
if (nargin < 3), skip = 1; end; % 1 => loadfits, no subsampling

%
% Computation
% 
% must initialize outputs in case no match occurs
im = []; rgn = []; rgnC = zeros(0,4); rgnT = []; cats = []; rgnR = {};
keys_found = {}; % image types found so far
fns = {}; % file names
fp = fopen(fn, 'r');
if (fp < 0), error('could not read <%s>', fn); end;
rgn_fn = '';
% load image components
nI = 0; % number of images so far
lin = fgetl_nocomment(fp);
while ischar(lin),
  [key1,lin] = strtok(lin); % first component
  key1fn = strtok(lin); % second component
  % enable the read if this key is wanted
  % formerly was: strncmp(keys, key1, 3)
  if isnumeric(keys) | (length(key1) >= 3 & any(strcmp(keys, key1))),
    if strcmp(key1, 'region'),
      % this is the regions file: only take the first one
      if isempty(rgn_fn),
	rgn_fn = key1fn;
      end;
    elseif ~isempty(strfind(key1, '.name')),
      % just report the name in this case
      fns{end+1} = key1fn;
    else
      % it is an image file: load it
      keys_found{nI+1} = key1; % save the image type (string)
      if length(regexp(key1fn, '\.fits[^/]*$')) == 1,
        if skip > 0
	  [temp,err] = loadfits(key1fn, abs(skip)); % get FITS image
        else,
	  [temp,err] = loadfitsa(key1fn, abs(skip)); % get FITS image
        end;
	if (err == 0),
	  nI = nI + 1; % #images
	  im = cat(3,im,temp);
	end;
      else,
	% try imread
	try,
	  temp = imread(key1fn); % get image, PPM, tiff, etc
	  nI = nI + 1; % #images
	  % take the easy way out when promoting classes
	  if length(im) > 0 && strcmp(class(temp),class(im)) == 0, 
	    im = double(im); temp = double(temp);
	  end;
	  im = cat(3,im,temp);
	  err = 0; % is OK
	catch,
	  err = 1; % not OK: just give a consistent message
	end;
      end;
      % report errors
      if err,
	error('failed to load file %s tag %s', key1fn, key1);
        fclose(fp);
      end;
    end; 
  end;
  lin = fgetl_nocomment(fp);
end;
fclose(fp);
% 
% load regions (boxes and catfiles)
% 
if isempty(rgn_fn),
  % nothing there; just return
  rgn = []; % no regions
  rgnC = []; % no region coords
  rgnT = []; % no region times
  cats = []; % no catfiles
  return;
end;  
% decide on old vs new region file
if     regexp(rgn_fn, '\.rgn$') > 0, rgn_type = 1; % old SAOtng .rgn
elseif regexp(rgn_fn, '\.reg$') > 0, rgn_type = 2; % new ds9 .reg
else   error('unrecognized region file type <%s>', rgn_fn); end;
% open file
fp = fopen(rgn_fn, 'r');
if (fp < 0), error('could not read <%s>', rgn_fn); end;
% 1: scan header for cat-files
% this works for either rgn_type
% set cat-file-list to null
for k = [1:nI], cats{k} = ''; end;
lin = fgetl(fp);
while ischar(lin),
  % break if we leave header
  if length(lin) == 0, break; end;
  if lin(1) ~= '#', break; end;
  lin = lin(2:end); % kill #
  % separate into key: value if possible
  [key1,lin] = strtok(lin); % find first word (keyword?)
  if ~isempty(key1) & any(key1 == ':') & length(key1) > 3, % if yes: keyword
    match = find(strncmp(keys_found, key1, 3)); % look for match
    if length(match) == 1, % found it
      cats{match} = strtok(lin); % set the corresponding catfile up
    end;
  end;
  lin = fgetl(fp);
end;
% 2: loop to get times and boxes
% we just reset the fp because there are several ways to exit the loop above
frewind(fp); 
if rgn_type == 1,
  lin = fgetl_nocomment(fp); % skip beyond first box, which is a dummy
elseif rgn_type == 2,
  lin = fgetl_nocomment(fp); % skip beyond global properties
  % check that the line really is global props; if not, back up
  if isempty(regexp(lin, '^global')), frewind(fp); end;
end;
nR = 1; % number of regions
while 1,
  % strategy is to put the date_time, coords, and box into lin1, lin2, lin3
  % for both file formats, then pick up from there
  if rgn_type == 1, 
    % Old .rgn region format
    %   # 1996/05/19 19:08:35
    %   # 2/5+4/5@(241:288,465:512)
    %   +box(40.5,40.5,48.0,48.0,0) # yellow
    lin1 = fgetl(fp); % date, time
    lin2 = fgetl(fp); % coordinates
    lin3 = fgetl(fp); % box spec
    if ~ischar(lin1) | ~ischar(lin2), break; end; % no more lines
    % clean up the pieces
    lin1 = regexprep(lin1, '^# *', ''); % strip comment mark
    lin2 = regexprep(lin2, '^# *', ''); % strip comment mark
    lin3 = regexprep(lin2, '^[+-]', ''); % strip leading + or -
  else,
    % new .reg format
    % box(930.0,544.5,97.0,96.0,0) # color=magenta 
    %   text={1996/11/17 23:58:35; 2/3@(937:1033,17:112)}
    lin = fgetl(fp); % box(...) # text={date time; coords}
    if ~ischar(lin), break; end; % no more lines
    if lin(1) == '#', continue; end; % skip comments
    % grab the contents of the {} in ...#...text={FOOBAR}...
    linT = regexprep(lin, '.*#.*text=\{([^\}]*)\}.*', '$1');
    lin1 = regexprep(linT, ';.*', ''); % date time: kill from first ; back
    lin2 = regexprep(linT, '.*; *', ''); % n/N@(...): kill up to last ;
    lin3 = regexprep(lin, '\s*#.*', ''); % strip from # back to obtain box
  end;
  % time is in lin1
  tnums = sscanf(lin1, '%d/%d/%d %d:%d:%f');
  if (length(tnums) ~= 6), 
    error('bad time in %s box %d (%s)', rgn_fn, nR, lin1);
  end;
  rgnT(nR) = datenum(tnums(1),tnums(2),tnums(3),tnums(4),tnums(5),tnums(6));
  % coords are in lin2
  atsign = strfind(lin2, '@');
  if length(atsign) ~= 1, 
    error('bad format in %s box %d (%s)', rgn_fn, nR, lin2);
  end;
  cnums = sscanf(lin2(atsign:end), '@(%d:%d,%d:%d)');
  if (length(cnums) ~= 4), 
    error('bad coordinates in %s box %d (%s)', rgn_fn, nR, lin2);
  end;
  rgnR{nR,1} = lin2(1:atsign-1); % n/N or n/N+m/M, etc. part
  rgnC(nR,:) = cnums(1:4)'; % put in as-is
  % region is in lin3
  rgnnums = sscanf(lin3, 'box(%f,%f,%f,%f'); % the fifth number is rotation
  if (length(rgnnums) ~= 4), 
    error('bad box in %s box %d (%s)', rgn_fn, nR, lin3);
  end;
  % put rgn in matlab coords: x -> n; y -> m
  rgn(nR,:) = [ceil(rgnnums(2:-1:1)-rgnnums(4:-1:3)/2) ; ...
	      floor(rgnnums(2:-1:1)+rgnnums(4:-1:3)/2) ]';
  nR = nR + 1;
end;
fclose(fp);
% change to vector orientation
cats = cats';
rgnT = rgnT';
return;


% get lines, discarding comments
function lin=fgetl_nocomment(fp)

lin = '#';
while ischar(lin) && ((length(lin) == 0) || (lin(1) == '#')),
  lin = fgetl(fp);
end;
return;
