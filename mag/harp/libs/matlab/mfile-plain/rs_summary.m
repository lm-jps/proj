function varargout = rs_summary(series, varargin)
% rs_summary	find: jsoc_info op=rs_summary ds=series
%
% [res,msg] = rs_summary(series, method, retries)
% * Invoke jsoc_info, op=rs_summary, using the given series.
% * See jsoc_info_request for more on the method and retries arguments.
% * If the second output, msg, is not requested, errors will result if
% there is no response from JSOC.  If msg is requested, it will be 
% empty if no error, or a descriptive string if there was an error, but
% an error will not be raised for communication or parsing failures.
% * So, if you request the msg output, you *must* check it to be sure it
% is empty, before using res.
%
% Usage: 
%   res = rs_summary('hmi.M_720s')
%   [res,msg] = rs_summary('hmi.M_720s','shell')
% 
% Inputs:
%   string series
%   opt string method
%   opt real retries
% 
% Outputs:
%   struct res
%   string msg
%    
% See Also: jsoc_info_request


% 
% Error checking
% 
if nargin < 1 || nargin > 3,
  error('Need 1 to 3 input arguments');
end

% set up output args (length must pass down to jsoc_info)
varargout = cell(1, max(1,nargout));

% allow jsoc_info_request to handle argument defaults
% note: params is always {}
[varargout{:}] = jsoc_info_request('rs_summary', series, {}, varargin{:});

return
end
