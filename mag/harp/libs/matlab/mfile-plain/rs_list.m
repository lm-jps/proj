function varargout = rs_list(series, varargin)
% rs_list	return result from jsoc_info op=rs_list ds=series
%
% [res,msg] = rs_list(series, params, method, retries)
% * Invoke jsoc_info, op=rs_list, using the given query.  This can 
% return lists of keywords, or segments, or links, or any combination
% thereof.
% * The params argument is an even-length cell array containing 
% pairs of parameters for rs_list, most importantly: 'key' for a keyword list,
% 'seg' for segments, and 'link' for links.
% * If method contains 'web' (the default) the CGI gateway speaking JSON 
% over http will be used.  This works remotely.  If method is 'shell', 
% the system jsoc_info command will be run and its output intercepted.
% * If the second output, msg, is not requested, errors will result if
% there is no response from JSOC.  If msg is requested, it will be 
% empty if no error, or a descriptive string if there was an error, but
% an error will not be raised for communication or parsing failures.
%
% Usage: 
%   res = rs_list('hmi.M_720s[$]', {'key', 'T_REC,QUALITY'})
%   [res,msg] = rs_list('hmi.M_720s[2010.06.01/1d]', {'key','T_REC','seg','**ALL**'})
% 
% Inputs:
%   string query
%   opt string method
% 
% Outputs:
%   struct res
%   string msg
%    
% See Also:

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 25 Oct 10.
% Copyright (c) 2010.  All rights reserved.

% 
% Error checking
% 
if nargin < 1 || nargin > 4,
  error('Need 1 to 4 input arguments');
end

% set up output args (length must pass down to jsoc_info)
varargout = cell(1, max(1,nargout));

% allow jsoc_info_request to handle argument defaults
% note: params is always {}
[varargout{:}] = jsoc_info_request('rs_list', series, varargin{:});

% always can count on the result
res = varargout{1}; 

% if status OK, validate res
must_exist = {'keywords', 'segments', 'links'};
msg = '';
EFMT = sprintf('%s: Error at %s: %%s.', mfilename, datestr(now));
if ~isempty(res) && res.status == 0,
  % OK: check fields
  if ~isfield(res, 'count'),
    msg = sprintf(EFMT, 'Missing count field in result, this should not happen.');
  elseif res.count > 0,
    % yes, it does not return other fields if count is zero
    not_exist = ~isfield(res, must_exist);
    if any(not_exist),
      missing = sprintf('%s,', must_exist{not_exist});
      msg = sprintf(EFMT, ['Missing field(s) in result: ' ...
                          missing ...
                          ' this should not happen.']);
    end;
  end;
end;

% allow error() exit if true
do_err_out = (nargout < 2);

% act on msg if it is nonempty
if ~isempty(msg),
  if do_err_out,
    error(msg);
  elseif isempty(varargout{2}),
    varargout{2} = msg;
  end;
end;

return
end
