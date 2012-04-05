function [q,ok,msg]=hmiquality(trec)
%hmiquality	find quality parameter from HMI T_REC index
% 
% [q,ok,msg]=hmiquality(trec)
% * Given a T_REC index, finds the quality parameter which lives
% in the QUALITY keyword.
% * If two outputs are given, also finds ok, a boolean value which
% tells if the quality is sufficient for most operations.
% * If three outputs are given, also returns a descriptive error
% message upon error, or empty if no error.  In this case, no errors
% are thrown.
% * A typical T_REC is: 2010.07.03_13:12:00_TAI
% 
% Inputs:
%   string trec(nf);        -- a valid time index
% 
% Outputs:
%   real q;
%   bool ok;
%   string msg;
% 
% See Also:

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 17 Feb 2011
% Copyright (c) 2011.  All rights reserved.

% 
% Error checking
% 
if all(nargin  ~= [1]), error ('Bad input arg number'); end;
if all(nargout ~= [0 1 2 3]), error ('Bad output arg number'); end;

do_err_out = (nargout < 3); % otherwise, return an error message

%
% Computation
% 

% default return values in case of early exit
q = [];
ok = [];
msg = '';

% look for metadata in this data series
nrt_mode = hmi_property('get', 'nrt_mode');
if nrt_mode,
  parent = 'hmi.M_720s_nrt';
else,
  parent = 'hmi.M_720s';
end

% this mask contains the OR of:
%   QUAL_NODATA:  0x80000000 -- the MSB is set when there is no LOS observable produced
%   QUAL_ECLIPSE: 0x00000200 -- at least 1 lev1 observable had an eclipse
% See http://jsoc.stanford.edu/doc/data/hmi/Quality_Bits/QUALITY.txt -- you are looking
% for the lines for the 12-min averaged IQUV -- for lev 1.5 quality bits (not the same as lev1).
% the eclipses are seen around 2011.03.04_12:36_TAI, and possibly again in feb-mar 2012.
bad_quality = hex2dec('80000200');

% get keywords
series = sprintf('%s[%s]', parent, trec);
if do_err_out,
  % do not ask for msg argument -- will error out if QUALITY is not present
  a = rs_list(series, {'key', 'QUALITY'});
else,
  % ask for msg argument
  [a,msg] = rs_list(series, {'key', 'QUALITY'});
  if ~isempty(msg),
    % q and ok already set up
    return;
  end;
end;  

% insist on exactly one response
if a.count ~= 1,
  msg = sprintf('Response for %s for QUALITY had count = %d, needed 1', series, a.count);
  if do_err_out,
    error(msg);
  else,
    % q and ok already set up
    return;
  end;
end;

q1 = jsoc_cell2struct_keys(a.keywords);

if ~isfield(q1, 'quality'),
  msg = sprintf('Keyword response for %s for QUALITY not present', series);
  if do_err_out,
    error(msg);
  else,
    % q and ok already set up
    return;
  end;
end;

if length(q1.quality) ~= 1,
  msg = sprintf('Keyword response for %s for QUALITY wrong length', series);
  if do_err_out,
    error(msg);
  else,
    % q and ok already set up
    return;
  end;
end;

qstr = q1.quality{1}; % the string
qchr = qstr(3:end); % mask off the leading 0x
q = hex2dec(qchr); % 32-bit integer stored as a double

% insert the "ok for tracking purposes" flag
ok = bitand(q, bad_quality) == 0;

return;
end
