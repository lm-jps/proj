function [status,fn2] = avi2mp4(fn1, opts, varargin)
%avi2mp4	convert avi to mp4
% 
% [status,fn2] = avi2mp4(fn1, opts, varargin)
% * Convert an avi movie file fn1 into an mp4 file fn2, according to
% various options, returning nonzero status on error.
% * This routine attempts *not* to throw an error itself.  So, you
% should check the status.
% * If ffmpeg is available, we use that, otherwise qt-avi2mp4, which
% is my own script taking no options.
% * Options to are passed in via a 1x1 structure opts.  The field
% names in opts must correspond to the converter argument names, and 
% the values must be the corresponding value (the value associated
% with flag options should be left empty).
% * Sometimes it is preferable to add options on the same command 
% line.  To do this, let varargin be a series of tuples containing
% the argument name, and the value. 
% * For example: 
%   >> avi2mp4('/tmp/in.avi', 'qscale', 2, 'r', 20)
% would give options -qscale 2 -r 20 to ffmpeg.
% * Some pseudo-options are supported.  They affect this routine, not
% the converter, and can be passed in thru opts or varargin.
%   - A flag 'dryrun' which does not execute the converter command, but
%   returns it thru fn2.
%   - A flag 'vv' which shows the command run, and shows the command
%   transcript (stdout and stderr).  Can be useful with the 'v' option,
%   which tells ffmpeg itself to be verbose.
%   - A string-valued option 'path' giving the path to the converter 
%   binary.
%   - A flag 'clean' which deletes the input file (only if the conversion
%   succeeds).
% 
% Inputs:
%   string fn1
%   opt struct opts = []
%   varargin
% 
% Outputs:
%    int status  -- you must request it
%    opt string fn2
%    
% See Also:

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 15 Jun 09.
% Copyright (c) 2009.  All rights reserved.

% default path to command
% (can be overridden)
cmds = {'ffmpeg', '/usr/local/bin/ffmpeg.n02', 'qt-avi2mp4'};
cmdpath_def = '';
for cmd1 = cmds,
  [s,r] = unix(sprintf('which %s', cmd1{1}));
  if s == 0,
    cmdpath_def = r(1:end-1); % strip newline
    break;
  end;
end;

% indicate error until we succeed
status = 1;
fn2 = '';

% 
% Error checking
% 
if (nargin < 1) || (nargin > 2 && rem(nargin, 2) ~= 0), 
  % 1 or 2 ok, and 4, 6, ... OK
  error('%s: Bad input arg number', mfilename); 
end;
if all(nargout ~= [1 2]), 
  error('%s: Bad output arg number', mfilename); 
end;
% variable input args
if nargin == 1 || isempty(opts), opts = struct(); end; % need a 1x1 struct
if ~isstruct(opts) || length(opts) ~= 1, 
  error('%s: opts is a structure scalar', mfilename); 
end;
% plug trailing input args into opts
Ntrail = max(0, (nargin-2)/2); % number of trailing args
for i = 1:Ntrail,
  opts.(varargin{2*i-1}) = varargin{2*i};
end;

% set up and check for command path
if isfield(opts, 'path'),
  cmdpath = opts.path;
else,
  cmdpath = cmdpath_def;
end;
if ~exist(cmdpath, 'file'),
  fprintf('%s: Error: Could not find converter executable at "%s"\n', ...
          mfilename, cmdpath);
  return;
end;

% opts, based on path
[junk,cmdname]=fileparts(cmdpath);
if strcmp(cmdname, 'ffmpeg'),
  opts.i = fn1;
  opts.qscale = 5;
  opts.y = '';  % overwrite output if present, just a flag
elseif strcmp(cmdname, 'qt-avi2mp4'),
  % nothing special
else,
  fprintf('%s: Warning: No builtin defaults for converter <%s>\n', ...
          mfilename, cmdname);
end;

% ensure input file exists
if ~exist(fn1, 'file'), 
  fprintf('%s: Error: Could not find input file in <%s>\n', mfilename, fn1);
  return;
end;

% make output file
fn2 = regexprep(fn1, '\.avi$', '.mp4');
if strcmp(fn1, fn2) == 1,
  fprintf('%s: Error: Output and input (%s) are the same (input not .avi?)\n', ...
          mfilename, fn1);
  return;
end;
% kill output file if it exists
% (we check for clean termination by looking for the file)
if exist(fn2, 'file'), 
  delete(fn2);
end;

%
% Computation
% 

% cons up argument string
argstr = '';
for tag1 = fieldnames(opts)',
  tag = tag1{1}; % extract cell contents
  % skip special tags
  if any(strcmp(tag, {'dryrun', 'vv', 'clean'})),
    continue;
  end;
  % process the (potential) argument
  tagval = opts.(tag); 
  argstr = [argstr ' ' argfmt(tag, tagval)];
end;

% did not bother to generalize this one -- for ffmpeg, outfile has to be last
if strcmp(cmdname, 'ffmpeg'),
  argstr = [argstr ' ' fn2 ];
end

% final command
cmd = [cmdpath argstr]; % argstr begins with ' '

% vv flag
if isfield(opts, 'vv'),
  fprintf('%s: Verbose: Running command:\n%s\n', mfilename, cmd);
end;
% run command
if isfield(opts, 'dryrun'),
  % dry run: return command string
  status = 0;
  fn2 = cmd;
  return;
else,
  % regular run
  [s,w] = unix(cmd);
  % if vv, echo stdout
  if isfield(opts, 'vv'),
    fprintf('%s: Verbose: Command output follows:\n', mfilename);
    disp(w);
  end;
  % exit status signaled error
  if s ~= 0,
    fprintf('%s: Error running command, returning.\n', mfilename);
    status = 1; % error
    return;
  end;
end;

% check that output file exists
if ~exist(fn2, 'file'), 
  fprintf('%s: Error: Could not find outfile in <%s>\n', mfilename, fn2);
  status = 1; % error
  return
end;

% clean flag -- only runs when the command ran OK
if isfield(opts, 'clean'),
  delete(fn1);
end;

status = 0; % ok
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Format one argument
%
function str = argfmt(tag, val)

% build arg string
if isempty(val),
  str = sprintf('-%s', tag);
elseif isnumeric(val),
  str = sprintf('-%s %g', tag, val);
else,
  str = sprintf('-%s %s', tag, val);
end;
% returns str
return;
end; % argfmt


end % hmm_fit

