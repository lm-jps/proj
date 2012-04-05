function varargout = series_struct(series, varargin)
% [res,msg] = series_struct(series, method, retries)
% * Invoke jsoc_info, op=series_struct, using the given series.
% * See jsoc_info_request for more on the method and retries arguments.
% * If the second output, msg, is not requested, errors will result if
% there is no response from JSOC.  If msg is requested, it will be 
% empty if no error, or a descriptive string if there was an error, but
% an error will not be raised for communication or parsing failures.
% * So, if you request the msg output, you *must* check it to be sure it
% is empty, before using res.
%
% Usage: 
%   res = series_struct('hmi.M_720s')
%   series_struct('hmi.M_720s', 'verbose')
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

% verbosity
if nargin >= 2, 
  verbose = ~isempty(strfind(varargin{1}, 'verbose'));
else,
  verbose = false;
end;

% set up output args (length must pass down to jsoc_info)
varargout = cell(1, max(1,nargout));

% allow jsoc_info_request to handle argument defaults
% note: params is always {}
[varargout{:}] = jsoc_info_request('series_struct', series, {}, varargin{:});

if verbose && isfield(varargout{1}, 'status') && varargout{1}.status == 0,
  results = varargout{1};
  
  % Keywords {'name': , 'type': , 'units': , 'note': }
  fprintf('\nKeywords:\n');
  for k=1:length(results.keywords)    
      fprintf ('  %-30s\t%s\n',results.keywords{k}.name, results.keywords{k}.note);    
  end

  % Segments
  fprintf('\nSegments:\n');
  for k=1:length(results.segments)    
      fprintf ('  %s\n',results.segments{k}.name);
      fprintf('\n');
      disp(results.segments{k})
  end

  % Links
  fprintf('\nLinks:\n');
  for k=1:length(results.links)    
      fprintf ('  %s\n',results.links{k}.name); % turmon: link-> links
  end

  % DB Indexes
  fprintf('DB Index:\n');
  disp(results.dbindex)

  % Intervals
  fprintf('\nInterval:\n');
  disp(results.Interval)

  fprintf('\n');
end;

return
end

