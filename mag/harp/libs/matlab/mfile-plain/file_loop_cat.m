function [fn,res,fnfail,logname] = file_loop_lst(lstname,skip,cmd,varargin)
%file_loop_lst	loop over files, finding cmd(varargin)
% 
% function [fn,res,fnfail,logname] = file_loop_lst(lstname,skip,cmd,varargin)
% * Generic loop over the files in a .lst file, evaluating cmd on
% each, returning results in res.  Also works for .cat files.
% If invoked as file_loop_lst, only the filename of cat/lst
% files is accessible.  
% * If invoked as file_loop_cat, the time is also read and given 
% unto the command 'cmd' after the filename argument.  At least one 
% cat-file must be present, and cat-file times must agree to
% within one second.
% * Comment and blank lines are treated properly.  Like the "map"
% construct in functional programming.
% * Specifically, for each filename F and time T in the .lst/.cat file, 
% >> ans = cmd(F, T, loadfits(F,skip(1)), varargin{:});
% is evaluated.  The filenames used are returned in fn, and the
% list of all ans's is returned in res.  As remarked above,
% the time is omitted if file_loop_lst is used.
% * If a cell array of lst files is given, the "map" is vectorized
% by passing a cell array of filenames and images in the command
% above, one cell per lst file.  This is similar to the python 
% construct:  for (x,y) in zip(xlist,ylist): cmd((x,y)).
% * If skip(1) < 0, we use loadfitsa() rather than loadfits(),
% with the absolute value of the skip specified.
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
%   string lstname  -- filename
%   int skip(3);    -- skip(2)=1, skip(3)=0 if not supplied
%   function cmd;   -- string, or inline
%   cell varargin;
% 
% Outputs:
%   cell string fn{nfile}; -- cell array of strings
%   cell res{nfile};       -- cell array of results
%   cell fnfail{nfail};    -- files that failed (typically empty)
%   char logname;          -- log file name in /tmp
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
if ischar(lstname),
  lsts = {lstname}; % a one-entry cell array
else,
  lsts = lstname; 
end;
nF = length(lsts);
% assume cat files (with times) or just lst files?
cat_mode = ~isempty(strfind(mfilename, 'cat'));

%
% Computation
% 

tic; % tic/toc is used later
% Read files
ftype = zeros(1,nF); % 0=unread, 1=lst, 2=cat
for f = [1:nF],
  lstname = lsts{f};
  ifp = fopen(lstname, 'r');
  if ifp < 0
    error('image .lst file <%s> not found', lstname);
  end;
  % try opening as a cat-file (with times)
  % "%n" seems to include floating point numbers (bad), but converts
  % into a double (good).  The datenum conversion needs doubles
  % to capture milliseconds, so we use %n for month numbers, etc.
  Cr = textscan(ifp, '%s%n/%n/%n%n:%n:%f', 'CommentStyle', '#');
  if feof(ifp),
    % it *was* a cat-file; get the times
    CrT = datenum([Cr{2:end}]); % concatenating into N-by-6 matrix
    ftype(f) = 2; % it was a .cat
  else,
    frewind(ifp); % not at eof => conversion failed.  Back to start.
  end;
  % if not read already, try for lst file
  if ftype(f) == 0,
    Cr = textscan(ifp, '%s', 'CommentStyle', '#');  % lst-file format
    if ~feof(ifp), 
      error('file %s not .lst/.cat file (.lst error at %d bytes)',...
                    lstname, ftell(ifp));
    end;
    % converted OK as .lst
    CrT = nan(length(Cr{1}), 1); % times are NaNs
    ftype(f) = 1; % it was a .lst
  end;
  fclose(ifp); % file read OK
  % compile file names and times
  if f == 1, 
    nI = length(Cr{1}); % first sets nI, the number of images
    ImN = cell(nF, nI); % image names
    ImT = zeros(nF, nI); % image times
  end; 
  % provide diagnostic for length mismatches
  if length(Cr{1}) ~= nI,
    error('lst length mismatch.\n%s: %d files\n%s: %d files\n',...
                  lsts{1}, nI, lstname, length(Cr{1}));
  end;
  % these are used later
  ImN(f,:) = Cr{1}'; % filenames
  ImT(f,:) = CrT'; % file times or nan
end;
% validate file times; loop below allows good error diagnostics
ref = find(ftype == 2, 1); % first file with times
for f = find(ftype == 2),
  if f == ref, continue; end;
  % times are in units of days, 1/24*3600 is 1 sec
  crit = abs(ImT(f,:) - ImT(ref,:)) > 2.001/(24*3600);
  if any(crit),
    first = find(crit, 1);
    error([sprintf('file time mismatch %s vs %s\n', lsts{ref}, lsts{f}) ...
           sprintf('%d bads, first is %s vs %s', nnz(crit), ...
                   datestr(ImT(ref,first), 'yyyy/mm/dd HH:MM:SS.FFF'), ...
                   datestr(ImT(  f,first), 'yyyy/mm/dd HH:MM:SS.FFF'))]);
  end;
end;
ImT = ImT(ref,:); % just take one

% log file (attempt to be unique)
logname = sprintf('/tmp/fileloop.%06d-%03d.log', floor(1e6*rem(now,1)), fix(1000*rand));
efp = fopen(logname, 'w');
if efp < 0
  error('log file <%s> not writable', logname);
else
  if isa(cmd, 'function_handle'), cmdAsString = func2str(cmd); 
  else, cmdAsString = cmd; end;
  fprintf(1, 'logging to <%s>\n', logname);
  fprintf(efp, '#\n# Log file created %s\n', datestr(now));
  fprintf(efp, '# Command: %s\n', cmdAsString);
  fprintf(efp, '# Filenames: %s\n', lstname);
  fprintf(efp, '#\n\n');
end;

% looks ok, call initial hook if needed
if skip(3) ~= 0,
  % time input (optional) handled via conversion of cell arrays to arg lists
  if cat_mode,
    TimeArg = {-Inf}; % cell array of size 1
  else,
    TimeArg = cell(0); % empty cell array
  end;
  % call the command (don't record the result)
  if skip(1) ~= 0,
    res1 = feval(cmd, -Inf, TimeArg{:}, [], varargin{:});
  else
    res1 = feval(cmd, -Inf, TimeArg{:},     varargin{:});
  end;
end;  
% initialize to empty
fn = cell(length([1:skip(2):nI]),1);
fnfail = {}; % typically will stay empty
res  = []; res_accum = [];
nfile = 1; % number of files loaded so far
           % nb, nline is number of current file in total list
dots  = 0; % dots printed on this line
err = 0; % define this var
% loop over files
for nline = [1:skip(2):nI],
  ifile = ImN(:,nline); % cell of {file, file, ...}
  % display status update
  fprintf(1,'.'); dots = dots+1;
  if dots > 79; fprintf(1,'\n'); dots = 0; end;
  % note the toc within
  fprintf(efp, '# clock time = %s; dt = %.2f\n%s\n', datestr(now), toc, sprintf('%s\n', ifile{:}));
  tic
  % 
  % Get image(s) if needed
  % 
  x = cell(1,nF); % accumulated images
  for f = [1:nF], 
    if skip(1) > 0,
      [x1,err1] = loadfits (ifile{f},  skip(1));
    elseif skip(1) < 0,
      [x1,err1] = loadfitsa(ifile{f}, -skip(1));
    else,
      x1 = []; err1 = 0; % skip = 0 case
    end;
    % Check for an error
    if err1,
      % summarize
      fprintf(1,'\nFailed to load <%s>, raw file number %d\n', ...
              ifile{f}, nline);
      fprintf(efp, 'line %d failed on load of %s from %s\n', ...
              nline, ifile{f}, lsts{f});
      dots = 0;
      % record offending filenames
      if nF == 1,
        fnfail{end+1,1} = ifile{1};
      else,
        fnfail{end+1} = ifile;
      end;
      % skip processing it
      continue; % this does not concatenate a result instance
    end;
    x{f} = x1;
  end;
  % 
  % Do command 
  % 
  % filename and image inputs
  if nF == 1,
    FileArg = ifile{1}; % a string
    ImageArg = x{1}; % an array
  else,
    FileArg = ifile; % a cell-array of strings
    ImageArg = x; % a cell-array of arrays
  end;
  % time input (optional) handled via conversion of cell arrays to arg lists
  if cat_mode,
    TimeArg = {ImT(nline)}; % cell array of size 1
  else,
    TimeArg = cell(0); % empty cell array
  end;
  % call the command
  % would like a try/catch here, and to augment fnfail
  % (at least under some circumstances)
  if (skip(1) ~= 0)
    res1 = feval(cmd, FileArg, TimeArg{:}, ImageArg, varargin{:});
  else,
    res1 = feval(cmd, FileArg, TimeArg{:}, varargin{:});
  end;
  % do what we can to make the results commensurate
  if isstruct(res1),
    res1 = orderfields(res1);
  end;

  %
  % Append the result
  %
  % maintain filename list
  if nF == 1,
    fn{nfile} = ifile{1};
  else,
    fn{nfile} = ifile;
  end;
  nfile = nfile+1;
  % put individual results into res_accum, transfer over to the (larger) 
  % res every so often to avoid O(n^2) behavior.
  if isempty(res_accum),
    res_accum = res1;
  else,
    res_accum = cat(2, res_accum, res1);
  end;
  % append res_accum to res every so often
  if rem(nfile, max(10,ceil(sqrt(nI/skip(2))))) == 0,
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
  % time input (optional) handled via conversion of cell arrays to arg lists
  if cat_mode,
    TimeArg = {Inf}; % cell array of size 1
  else,
    TimeArg = cell(0); % empty cell array
  end;
  % call the command
  if skip(1) ~= 0,
    res1 = feval(cmd, Inf, TimeArg{:}, [], varargin{:});
  else
    res1 = feval(cmd, Inf, TimeArg{:},     varargin{:});
  end;
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

