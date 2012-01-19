function [image, err, params] = loadfitsa_jsoc(file, block, type)
%loadfitsa_jsoc	load array from file in FITS format, cached
%
% [image, err, params] = loadfitsa_jsoc(file, block, type)
% * This is a thin wrapper for loadfitsa which implements a
% fits-file cache for remote files.  That is, if given a 
% url-format file beginning with `http://jsoc', it first checks
% a local cache for the file and returns that if present.
% Otherwise, it loads the url and then caches it the resulting
% file.
% * File is the name of the fits file to be read, which is an
% arbitrary-dimensional fits file.  It is loaded into Matlab in
% the "natural" ordering -- the native fits ordering, Fortran
% progession of coordinates.  The standard fitsio naming shortcuts
% are accepted, including 'file.fits[1]' for the first image
% extension, useful for loading compressed images in bintables.
% * Block is the blocking factor; only elements 1, block+1, ...
% are loaded (in every dimension).
% * type is the destination Matlab type, e.g. 'double' or 'int32'.
% When using the integer types, be aware that the input FITS file
% must fit within the range of the given type, or there will be
% an error.
% * The resulting array is returned; this is 0x0 if an error occurs.
% * Error is 1 if there was an error, 0 otherwise.  An error
% synopsis is printed if error = 1. 
% * To allow robust scripting, if error is requested, a matlab 
% error is not raised even if a file-related error occurs.  The
% error flag is set, and a message is printed, but mexErrMsgTxt is
% not called.  (Exception: if the arguments are malformed, as 
% opposed to a file-reading problem, a mex error *is* raised.)
% * loadfitsa and savefitsa are inverse functions.
%
% Inputs:
%   string filename
%   opt int block = 1
%   opt string type = 'double'
%
% Outputs:
%   real image(...)
%   int err
%   real params(4)
%
% See Also: loadfitsa

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 02 Sep 2010
% Copyright (c) 2010.  All rights reserved.

% 
% Error checking
% 
% if all(nargin  ~= []), error ('Bad input arg number'); end;
% if all(nargout ~= [0 1 2]), error ('Bad output arg number'); end;
if nargin < 3, type = 'double'; end;
if nargin < 2, block = 1; end;

%
% Computation
% 

sentinel = 'http://jsoc';
NS = length(sentinel);

% bail out of regular-file case ASAP
if length(file) < NS || strcmp(file(1:NS), sentinel) == false,
  [image, err, params] = loadfitsa(file, block, type);
  return;
end;

% from here down, file begins with sentinel, http://...
% (i.e., we are now in the http case)

% cache setup
cacheroot = '/data/hmi/cache';
cache_enabled = exist(cacheroot); % at stanford, cache is not enabled

% make the cache filename
slashes = find(file == '/');
file_root = file(slashes(end)+1:end); % file.fits
file_path = file(slashes(3)+1:slashes(end)-1); % /path/to/dir
file_full = file(slashes(3)+1:end); % /path/to/dir/file.fits

cachefile = fullfile(cacheroot, file_full);
cachedir  = fileparts(cachefile); % just the enclosing dir

% load the file from the cache if possible
if cache_enabled && exist(cachefile, 'file'),
  [image,err,params] = loadfitsa(cachefile, block, type);
  return;
  % (if it was uncached, just fall through)
end;

% fetch a new copy, unblocked
[image,err,params] = loadfitsa(file, 1, type);
% if this load fails, no caching can happen, so just return
if err,
  return;
end;

% cache it, unblocked, if desired
if cache_enabled,
  % cachedir and cachefile are from above
  if ~exist(cachedir)
    mkdir(cachedir);
  end;
  error2 = savefitsa(image, cachefile, params(1), params(2), params(3), params(4));
  if error2,
    % not fatal
    warning('%s: cache operation to %s failed', mfilename, cachefile);
  end;
end;

% block the return value
if block > 1,
  image = image(1:block:end, 1:block:end);
end;

return;
