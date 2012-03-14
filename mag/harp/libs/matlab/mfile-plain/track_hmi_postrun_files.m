function ctag=track_hmi_postrun_files(ddir, tag, fnlog, fnfail)
%track_hmi_postrun_files	save tracker state and log files
% 
% ctag=track_hmi_postrun_files(ddir, tag, fnlog, fnfail)
% * ddir is a destination root directory, like /tmp/Tracks/jsoc;
% tag is an identifying tag (in current usage, the T_REC slice),
% fnlog is a log filename, and fnfail is a cell array of tuples
% containing a filename and a message.
% * Copies a tracker-state file in ddir to a tag-specific location
% below ddir.
% * Moves a transient log file (typically in /tmp) to a permanent 
% location below ddir.
% * Saves a list of failure filenames to an error file below ddir.
% * Returns the tag-specific location name, which has some 
% metacharacters stripped, so other names can be harmonious.
% 
% Inputs:
%   string ddir
%   string tag
%   string fnlog  -- file name
%   cell fnfail of cell(1,2) of string
% 
% Outputs:
%   string ctag
% 
% See Also:

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 22 Feb 2011
% Copyright (c) 2011.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [4]), error ('Bad input arg number'); end;
%if all(nargout ~= [0 1 2 3 5]), error ('Bad output arg number'); end;

%
% Computation
% 
% re-encode the tag to remove any weird chars (/ is the main issue)
% RE is of form: [^xyz] where xyz are good
% note - has to be the last char because it's a [] metacharacter
% note - we use system() on this string, so it has to have NO metacharacters
not_goodchar = '[^\w@+=%.,-]'; 
ctag = regexprep(tag, not_goodchar, '_');
% kill multiple _'s
ctag = regexprep(ctag, '_+', '_');
% kill _ at end, if any
ctag = regexprep(ctag, '_$', '');

% copy tracker state in ddir to another spot
fs = sprintf('%s/track-post.mat', ddir); % source file
fn = sprintf('%s/State/track-post-%s.mat', ddir, ctag);
[junk1,junk2] = mkdir(fileparts(fn));
[st,res] = system(sprintf('cp -pf %s %s', fs, fn));
if st > 0,
  error('Could not copy tracker state file');
end;

% move log in /tmp to a permanent spot
fn = sprintf('%s/Log/%s.log', ddir, ctag);
[junk1,junk2] = mkdir(fileparts(fn));
[st,res] = system(sprintf('mv %s %s', fnlog, fn));
if st > 0,
  error('Could not copy log file');
end;

% create error file if needed
fn = sprintf('%s/track-latest.err', ddir);
if exist(fn, 'file'),
  delete(fn); % ensure it is removed if not used
end;
if length(fnfail) > 0,
  fp = fopen(fn, 'w');
  if fp < 0,
    error('Could not create skipped-file file');
  end;
  fprintf(fp, '# HMI skipped-T_REC list\n');
  fprintf(fp, '# Created while running: %s\n', tag);
  fprintf(fp, '# Generated by %s on %s from matlab\n', ...
          getenv('USER'), datestr(now));
  fprintf(fp, '#\n');
  for k = 1:length(fnfail),
    fprintf(fp, 'T_REC = %s\n', fnfail{k}{1});
    fprintf(fp, 'Reason = %s\n\n', fnfail{k}{2});
  end;
  fclose(fp);
  % copy to a standard place
  fs = sprintf('%s/Log/%s.err', ddir, ctag);
  [junk1,junk2] = mkdir(fileparts(fn));
  [st,res] = system(sprintf('cp -pf %s %s', fn, fs));
  if st > 0,
    error('Could not copy skipped-file file');
  end;
end;

return
end
