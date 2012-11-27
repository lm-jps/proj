function ctag=track_hmi_clean_runname(tag)
%track_hmi_clean_runname	sanitize a tag to use as a filename
% 
% ctag=track_hmi_clean_runname(tag)
% * tag is an identifying tag (in current usage, the T_REC slice),
% and ctag is that tag with metacharacters/punctuation found in T_REC's,
% like /[]()$, replaced with un-doubled underscores.  We want no 
% shell metacharacters in ctag.
% * It's intended that the cleaned tag be used to identify a log file or
% a movie in a standard way.
% * So for example:
%   tag = 'hmi.M_720s[2010.11.03_10:00_TAI/12h][?camera>0?]'
%   ->
%   ctag = 'hmi.M_720s_2010.11.03_10_00_TAI_12h_camera_0'
% 
% Inputs:
%   string tag
% 
% Outputs:
%   string ctag
% 
% See Also: track_hmi_postrun_files

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 22 Feb 2011
% Copyright (c) 2011.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [1]), error ('Bad input arg number'); end;
%if all(nargout ~= [0 1 2 3 5]), error ('Bad output arg number'); end;

%
% Computation
% 
% re-encode the tag to remove any weird chars (/ is the main issue)
% RE is of form: [^xyz] where xyz are good; thus, @ is allowed.
% note: "-" has to be the last char because it's a [] metacharacter
% note: we may use system() on this string, so it has to have NO metacharacters
not_goodchar = '[^\w@+=%.,-]'; 
ctag = regexprep(tag, not_goodchar, '_');
% kill multiple _'s
ctag = regexprep(ctag, '_+', '_');
% kill _ at end, if any
ctag = regexprep(ctag, '_$', '');

if isempty(ctag), ctag = 'empty_tag'; end;

return
end
