function [trec,res,fnfail,logname] = jsoc_trec_loop(series,skip,cmd,varargin)
%jsoc_trec_loop	loop over data series, finding cmd(varargin)
% 
% [trec,res,fnfail,logname] = jsoc_trec_loop(series,skip,cmd,varargin)
% * Generic loop over the trec indexes in a data series, evaluating cmd on
% each, returning results in res.
% * Like the "map" construct in functional programming.
% * Specifically, for each trec F in the series,
% >> res = cmd(F, T, loadfitsa(F,skip(1)), varargin{:});
% is evaluated.  The trec's used are returned in trec, and the
% list of all ans's is returned in res.  As remarked above,
% the time is omitted if file_trex_loop is used.
% * If res is returned as a char, it is assumed that a non-fatal error
% caused that input to be skipped.  That res is taken to be the reason.
% A 1x2 cell of the offending T_REC and the reason is appended to
% fnfail.  The rationale is that char's should seldom be returned
% as a legitimate answer, because they don't "cat" well anyway.  To
% return a non-error string res, make it a 1x1 cell.
% * If a cell array of data series is given, the "map" is vectorized
% by passing a cell array of trec's and images in the command
% above, one cell per data series.  This is similar to the python 
% construct:  for (x,y) in zip(xlist,ylist): cmd((x,y)).
% * If skip(1) < 0, we use loadfitsa_jsoc() rather than loadfitsa(),
% with the absolute value of the skip specified.  If skip(1) > 0,
% we use loadfitsa on the SUMS filename.
% * If skip(1) = 0, it is assumed that the image is not needed, so
% >> ans = cmd(F, T, varargin{:});
% is evaluated, with results and filenames treated as usual.
% * If skip(2) is supplied, it gives the temporal subsampling, so
% that files 1, skip(2)+1, 2*skip(2)+1, is processed.  Default is 1.
% * If skip(3) is nonzero, the cmd is called at the beginning (-Inf)
% and end (Inf) of the iteration, with the image empty, and the
% filename not a string, but a real number, -Inf or +Inf.  This
% is a hook to allow initialization and cleanup.  No return value
% is expected or stored from these calls to the routine.
% 
% Inputs:
%   string series   -- data series
%   int skip(3);    -- skip(2)=1, skip(3)=0 if not supplied
%   function cmd;   -- string, or inline
%   cell varargin;
% 
% Outputs:
%   cell trec{nfile} of string; -- cell array of trec's
%   cell res{nfile};            -- cell array of results
%   cell fnfail{nfail} of cell; -- {trec,reason} for fails only
%   char logname;               -- log file name in /tmp
% 
% See Also: the matlab function_handle concept introduced by @

% 
% Error checking
% 
if nargin < 3, error ('Bad input arg number'); end
if all(nargout ~= [2 3 4]), error ('Bad output arg number'); end  
if isempty(skip), skip = 1; end;
if length(skip) == 1, skip = [skip 1]; end;
if length(skip) == 2, skip = [skip 0]; end;
% normalization of inputs
if ischar(series),
  lsts = {series}; % a one-entry cell array
else,
  lsts = series; 
end;
nF = length(lsts);
% assume cat files (with times) or just lst files?
cat_mode = ~isempty(strfind(mfilename, 'trec'));

%
% Computation
% 

% prepare for toc within loop below
tic;
% Extract from lsts:
%    T_REC specifier, data series names (hmi.mag_720s), segment names (mask)
[ok, trec_spec, dsnames, segnames] = jsoc_parse_lstnames(lsts);
if ~ok,
  % trec_spec contains error message in this case
  error('Error in a data series (%s)', trec_spec)
end;

%
% 1: Read reference series to get times
%
% retrieve T_REC's and data segments
a = rs_list(sprintf('%s key=T_REC', lsts{1}), 'web_access');
if a.status > 0,
  error('Could not get T_RECs for %s (rs_list failed)', lsts{1});
end;
% exhaustively test the response
if ~isfield(a, 'keywords') || ~iscell(a.keywords) || length(a.keywords) ~= 1,
  error('Could not get T_RECs for %s (empty response from rs_list)', lsts{1});
end;
if ~isfield(a.keywords{1}, 'values'),
  error('Could not get T_RECs for %s (no values from rs_list)', lsts{1});
end;
% this is the T_REC indexes as a cell array of strings
TRECs = a.keywords{1}.values';
nTREC = length(TRECs);   % number of TREC's
% number of images passed in per function call
if skip(1) ~= 0,
  nImg = length(dsnames); 
else,
  nImg = 0;
end;

%
% 2: Log file (attempt to be unique)
%
logname = sprintf('/tmp/trec_loop.%s.%06d-%03d.log', ...
                  getenv('LOGNAME'), floor(1e6*rem(now,1)), fix(1000*rand));
efp = fopen(logname, 'w');
if efp < 0
  error('log file <%s> not writable', logname);
else
  if isa(cmd, 'function_handle'), cmdAsString = func2str(cmd); 
  else, cmdAsString = cmd; end;
  fprintf(1, 'logging to <%s>\n', logname);
  fprintf(efp, '#\n# Log file created %s\n', datestr(now));
  fprintf(efp, '# Command: %s\n', cmdAsString);
  fprintf(efp, '# Series: %s\n', lsts{1});
  fprintf(efp, '#\n\n');
end;

%
% 3: Call initial hook if needed
%
if skip(3) ~= 0,
  % optional time, image inputs handled via conversion of cell arrays to lists
  if cat_mode, TimeArg = {-Inf}; else, TimeArg = {}; end; % 1x1 or 0x0 cell
  ImageArg = cell(sign(nImg)); % {}, or {[]}
  % call the command (don't record the result)
  res1 = feval(cmd, -Inf, TimeArg{:}, ImageArg{:}, varargin{:});
end;  

%
% 4: Loop over files
% 

% initialize to empty
trec   = {}; % successful trec's
fnfail = {}; % unsuccessful trec's
res  = []; res_accum = [];
dots  = 0; % dots printed on this line

for inxTREC = [1:skip(2):nTREC],
  trec1 = TRECs{inxTREC}; % trec as a string
  time1 = hmiobstime(trec1); % matlab datenum for the trec
  % display status update
  fprintf(1,'.'); dots = dots+1;
  if dots > 79; fprintf(1,'\n'); dots = 0; end;
  % note the call to toc below (first tic is above)
  fprintf(efp, '# clock time = %s; dt = %.2f\n%s\n', datestr(now), toc, trec1);
  tic; % (toc is above)
  %
  % Prepare the arguments
  % optional args handled via conversion of cell arrays to arg lists
  %
  % `filename' input
  FileArg = trec1; % a string
  % time input (optional)
  if cat_mode, TimeArg = {time1}; else, TimeArg = {}; end;
  % image input (optional)
  % ImageArg is a cell:
  %   nImg = 0 --> {}              --> ImageArg{:} unpacks to nothing
  %   nImg = 1 --> {x}             --> ImageArg{:} unpacks to x
  %   nImg > 1 --> {{x1, ..., xN}} --> unpacks to {x1, ..., xN}
  if nImg == 0,
    ImageArg = {};
  else,
    [ImageArg, err1] = jsoc_load_images(trec1, dsnames, segnames, skip(1));
    if ~isempty(err1),
      fprintf(efp, '# %s\n',   err1);
      fprintf(1,   '\n%s\n', err1);
      fnfail{end+1} = {FileArg, err1}; % 1x2 cell
      continue;
    end;
  end;
  % 
  % Do command 
  % 
  % TBD: a try/catch here, and augment fnfail ?
  %   (at least under some circumstances)
  res1 = feval(cmd, FileArg, TimeArg{:}, ImageArg{:}, varargin{:});
  % record offending indexes
  if ischar(res1),
    err1 = sprintf('%s: %s', trec1, res1);
    fprintf(efp, '# %s\n',   err1);
    fprintf(1,   '\n%s\n', err1);
    fnfail{end+1} = {FileArg, res1}; % 1x2 cell
    continue; % do not append this
  end;
  %
  % Append the result
  %
  trec{end+1} = FileArg;
  % do what we can to make the results commensurate
  if isstruct(res1),
    res1 = orderfields(res1);
  end;
  % put individual results into res_accum, transfer over to the (larger) 
  % res every so often to avoid O(n^2) behavior.
  if isempty(res_accum),
    res_accum = res1;
  else,
    res_accum = cat(2, res_accum, res1);
  end;
  % append res_accum to res every so often
  if length(res_accum) > max(10, sqrt(nTREC/skip(2))),
    if isempty(res),
      res = res_accum; res_accum = [];
    else,
      res = cat(2, res, res_accum); res_accum = [];
    end;
  end;
end
% take last bit from res_accum
if isempty(res),
  res = res_accum;
elseif ~isempty(res_accum),
  res = cat(2, res, res_accum); 
end;
% call final hook if needed
if skip(3) ~= 0,
  % optional time, image inputs handled via conversion of cell arrays to lists
  if cat_mode, TimeArg = {-Inf}; else, TimeArg = {}; end; % 1x1 or 0x0 cell
  ImageArg = cell(sign(nImg)); % {}, or {[]}
  % call the command (don't record the result)
  res1 = feval(cmd, +Inf, TimeArg{:}, ImageArg{:}, varargin{:});
end;  
% clean up
fprintf(1,'\n');
fprintf(efp, '\n#\n# Log file closed %s\n#\n', datestr(now));
fclose(efp);
if nargout < 4,
  % caller did not ask for the log file
  delete(logname);
  fprintf(1, 'Removed log file.\n');
end;
return;
end

%
% load image(s) from a list of data series names
%
% return images as cell array, and an error code
%

function [arg,err] = jsoc_load_images(trec1, dsnames, segnames, skip1)

arg = {}; % so early returns under error conditions are easy
Nimg = length(segnames);
x = cell(1,Nimg);
for i = 1:Nimg, 
  ds1 = dsnames{i};
  seg1 = segnames{i};
  obj1 = sprintf('%s[%s]', ds1, trec1);
  % query JSOC to locate the segment
  a = rs_list(sprintf('%s&seg=%s', obj1, seg1), 'web_access');
  if a.status > 0,
    err = sprintf('Failed to load %s, not found by rs_list\n', obj1);
    return;
  end;
  % exhaustively test the response
  % TBD: abstract this into jsoc_imload
  if ~isfield(a, 'segments') || ~iscell(a.segments) || length(a.segments) ~= 1,
    err = sprintf('Failed to load %s, no segments\n', obj1);
    return;
  end;
  s = a.segments{1};
  if ~isfield(s, 'values') || ~iscell(s.values) || length(s.values) ~= 1,
    err = sprintf('Failed to load %s, no segment values\n', obj1);
    return;
  end;
  v = s.values{1};
  if ~ischar(v) || isempty(v),
    err = sprintf('Failed to load %s, segment values bad\n', obj1);
    return;
  end;
  if v(1) ~= '/',
    % 2/2011, offline gives v = 'NoDataDirectory'
    err = sprintf('Failed to load %s, segment seems offline\n', obj1);
    return;
  end;
  % this is the SUMS filename
  file1 = v;
  % load the file
  if skip1 > 0,
    % direct read from filesystem, not http
    [x1,err1] = loadfitsa(file1, skip1);
  else,
    % load thru http, with optional cache (loadfitsa understands URLs)
    file1 = sprintf('http://jsoc.stanford.edu/%s', file1);
    [x1,err1] = loadfitsa_jsoc(file1, -skip1);
  end;
  % Check for an error
  if err1,
    % summarize
    err = sprintf('Failed to load %s from %s\n', obj1, file1);
    return;
  end;
  x{i} = x1;
end;
% if multiple images, wrap up again, so that they will unpack into
% a single cell array
if Nimg > 1,
  arg = {x}; % {{x1, x2, ..., xN}}
else
  arg = x; % {x1}
end;
err = '';
return
end
  
%%
% From each data series name in lsts, given as: 
%   series[TREC]{segname}
% where {segname} is optional, extract these three parts.
% We take pains to validate each series name, and to ensure
% the segname is in that series.
%
function [ok, trec_spec, dsnames, segnames] = jsoc_parse_lstnames(lsts)

ref = 1; % data series in `lsts' that serves as the t_rec reference
% default return values
ok = false;
trec_spec = '';
dsnames  = [];
segnames = [];

% extract trec_spec from one of lsts
trec_spec = regexprep(lsts{ref}, '^.*\[(.*)\].*$', '$1');
if strcmp(trec_spec, lsts{ref}),
  % no replacement => error
  trec_spec = sprintf('Could not find TREC in data series %s', lsts{ref});
  return;
end;

% extract dsnames en masse (allow ., _, and alphanumerics)
% note, a lsts entry could consist of a series name and nothing else.
dsnames = regexprep(lsts, '^([.\w]+).*$', '$1');
% validate each dsname
for i = 1:length(dsnames),
  ds1 = dsnames{i};
  % ensure it's made of the right chars
  if isempty(regexp(ds1, '^[.\w]+$')),
    trec_spec = sprintf('Could not find series name in %s', lsts{i});
    return;
  end;
  % check that it's present
  s = rs_summary(ds1, 'web_access');
  if s.status > 0,
    trec_spec = sprintf('Series %s appears invalid to series_info', ds1);
    return;
  end;
end;

% extract segnames, within {...}, or **ALL** if not
segnames = regexprep(lsts, '.*\{(.*)\}$', '$1');
% segname not given -> no replacement -> segname == lst -> change to **ALL**
segnames(strcmp(lsts, segnames) == true) = {deal('**ALL**')}; 
% validate each segname
for i = 1:length(segnames),
  ds1  = dsnames{i};
  seg1 = segnames{i};
  % get the series info (detailed version)
  s = series_struct(ds1, 'web_access');
  if s.status > 0,
    trec_spec = sprintf('Series %s appears invalid to series_struct', ds1);
    return;
  end;
  if ~isfield(s, 'segments') || length(s.segments) == 0,
    trec_spec = sprintf('Series %s has no segments', ds1);
    return;
  end;
  % check that seg1 is in the series
  if strcmp(seg1, '**ALL**'),
    if length(s.segments) ~= 1,
      trec_spec = sprintf('Series %s has != 1 segment, specify one', ds1);
      return;
    end;
  else,
    if nnz(cellfun(@(sg) strcmp(sg.name, seg1), s.segments)) ~= 1,
      trec_spec = sprintf('Series %s does not have segment %s', ds1, seg1);
      return;
    end;
  end;
end;

ok = true;
return
end

