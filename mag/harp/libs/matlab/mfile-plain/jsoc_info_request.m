function [res,msg] = jsoc_info_request(op, ds, params, method, retries)
% jsoc_info_request	result of: jsoc_info query op=op, ds=ds
%
% [res,msg] = jsoc_info_request(op, ds, params, method, retries)
% * Return the results of calling jsoc_info, op=op, ds=ds, using 
% the given parameters encoded as a cell array of strings: 
%    params = {key1, value1, key2, value2,...}.
% The resulting JSON response is parsed into a result res; an
% error message may be returned in msg.
% * This is intended to be an omnibus routine through which most
% queries to JSOC can be funneled.
% * By setting up the right params, this routine can return lists of 
% keywords, or segments, or links, or any combination thereof.
% * If method contains 'web' (the default) the CGI gateway at 
% jsoc2.stanford.edu speaking JSON over http will be used.  
% This works remotely.  If method contains 'shell', the system 
% jsoc_info command will be run and its output intercepted.
% Otherwise, the jsoc_host set up by hmi_property will be used.
% * If method contains 'verbose' (not the default), a printed summary
% of the request is shown.
% * There are three sorts of exit condition:
%  1: The request completed without error. A non-empty output will be
%     returned in res, and msg = '' if it was asked for.
%  2: The request completed but it found an error; e.g., a non-existent
%     data series.  In this case, res will be nonempty, and msg will 
%     contain a message.
%  3: The request did not complete, e.g. because of network failure.
%     In this case, res will be empty, and msg will contain a message.
% Note that, if the request completed, res will always be non-empty,
% because it will contain a status field.  You can check if msg is empty
% to see whether the request completed correctly.
% * If the second output, msg, is not requested, errors will be thrown
% in cases 2 and 3 above.  If msg is requested, it will be empty 
% under case 1, or a descriptive string otherwise.  If msg is requested,
% an error will not be raised for communication or parsing failures.
% * So, if you request the msg output, you *must* check it to be sure it
% is empty, before using res.
% * By default, the routine re-tries several times upon failure.  
% Specify retries for a different behavior.  The numbers given indicate
% delays-between-retries in seconds.  Negative numbers mean a message
% is printed if the corresponding retry is used. [] means no retries.
%
% Examples: 
%   jsoc_info_request('rs_summary', 'hmi.M_720s')
%   jsoc_info_request('rs_list', 'hmi.M_720s[$]', {'key', 'T_OBS'})
% 
% Inputs:
%   string op
%   string ds
%   opt cell query{1,2*NP} of string = {}
%   opt string method from hmi_property('get','jsoc_host');
%   opt retries = hmi_property('get', 'jsoc_retries');
% 
% Outputs:
%   struct res
%   string msg
%    
% See Also:
%   The JSOC wiki page for jsoc_info, at:
%   http://jsoc.stanford.edu/jsocwiki/AjaxJsocConnect

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 25 Oct 2011.
% Copyright (c) 2011.  All rights reserved.

% parameter defaults
% abs() = between-try delay in s, [] for no retries, <0 for announcement
%  (set to announce for long delays)
default_retries = hmi_property('get', 'jsoc_retries');

% 
% Error checking
% 
if nargin < 2 || nargin > 5,
  error('Need 2 ... 5 input arguments');
end
% input defaults
if nargin < 3, params = {}; end;
if nargin < 4, method = ''; end;
if nargin < 5, retries = default_retries; end;
% arg checks
if ~iscell(params) || rem(length(params), 2) ~= 0,
  error('params must be a cell array of even length');
end;
if ~all(cellfun(@ischar, params)),
  error('params must be a cell array of strings');
end;  
if ~isempty(method) && ~ischar(method),
  error('if supplied, method must be empty, or a string');
end;
if ~isnumeric(retries) || (~isvector(retries) && ~isempty(retries)),
  error('if supplied, retries must be a numeric vector, or empty');
end;

% service endpoints for shell AND web
base_service = 'jsoc_info';
shell_service = base_service;
web_service = ['cgi-bin/ajax/' base_service];

% parameter defaults
if strfind(method, 'shell'),
  host = 'shell';
elseif strfind(method, 'web'),
  % default host for generic 'web' -- use property manager's default
  % (expect: 'http://hmiteam:hmiteam@jsoc2.stanford.edu')
  host = hmi_property('get', 'jsoc_host', 'default');
else,
  % not explicitly given -- use property manager
  host = hmi_property('get', 'jsoc_host');
end;
use_shell = strcmp(host, 'shell');

% verbosity -- default to quiet
verbose = ~isempty(strfind(method, 'verbose'));

% allow error() exit if true
do_err_out = (nargout < 2);

%
% Computation
% 
all_params = {'op', op, 'ds', ds, params{:}};
% uniform format for error messages
% takes 3 string args: short error message, request, longer message (or empty)
EFMT = sprintf('%s: Error at %s: %%s.  Request:\n%%s\n%%s', mfilename, datestr(now));

% the loop exits with valid values for msg and res
max_try = length(retries)+1;
for ntry = 1:max_try,
  msg = ''; % no error message => OK
  res = []; % ensure there is no existing res
  % semantics of url/shell accessors:
  %   status = 1 implies success
  %   if success, json_content contains the json
  %   if not, json_content contains a descriptive message, if any
  %   req contains the access command used
  if ~use_shell,
    endpoint = sprintf('%s/%s', host, web_service);
    [json_content,status,req] = urlread_jsoc(endpoint, 'get', all_params);
    if status == 0, json_content = ['urlread_jsoc failed: ' json_content]; end;
  else,
    % note, to conform with semantics above, status = 1 for success
    [json_content,status,req] = shell_jsoc(shell_service, all_params);
  end;
  if verbose && ntry == 1,
    fprintf('%s: request was: %s\n', mfilename, req);
  end;
  % if request failed, set up error message; res already empty
  if status == 0, 
    msg = sprintf(EFMT, 'JSOC request failed', req, json_content);
  end;
  % if request returned OK, parse the JSON response
  if isempty(msg),
    % this sets up res
    [res,msg1] = parse_json(json_content);
    if ~isempty(msg1),
      % retry: bad JSON
      res = []; % bad JSON => ensure res is not junk
      msg = sprintf(EFMT, 'Could not parse json reply', req, msg1);
      % currently written to be fatal, subject to change
      if do_err_out, error(msg); else, return; end;
    end;
  end;
  % if OK, check the status wrapped in JSON (status type is double)
  % status decoder ring: 
  %   0 if OK, 1 if series not found, -1 if the backend process was terminated
  % hence, do not retry if status = 1
  if isempty(msg),
    if res.status > 0, 
      % make this return immediately: we have our answer
      msg = sprintf(EFMT, 'JSOC returned error within JSON status', ...
                    req, res.status);
      % preserve res in this case; it's valid
      if do_err_out, error(msg); else, return; end;
    elseif res.status < 0,
      % not fatal, allow retry
      msg = sprintf(EFMT, 'JSOC returned cancellation error (status < 0)', ...
                    req, res.status);
      % preserve res in this case; it's valid
    end;
  end
  % response found to be OK: break loop (res is OK)
  if isempty(msg), break; end;
  % retry mechanism: a bit baroque
  if ntry < max_try,
    % more tries remain
    msg1 = 'Retrying after pause.';
    delay = retries(ntry);
  else,
    % the last try failed
    msg1 = 'Giving up.';
    delay = -0.001; % tiny delay, with announcement
  end;
  if delay < 0,
    fprintf('%s: Query %d of %d failed.  %s\n', mfilename, ntry, max_try, msg1);
  end;
  pause(abs(delay)); % can have delay < 0
end;
% at loop exit, msg and res are valid already

if ~isempty(msg) && do_err_out,
  error(msg);
end;
return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% shell_jsoc: run shell command and scoop up JSON from stdout
%
function [json_content,status,req] = shell_jsoc(cmd, pars)

% cons up arg list
args = '';
for i = 1:2:length(pars),
  key = shell_encode(pars{i});
  val = shell_encode(pars{i+1});
  args = [args ' ' key '=' val];
end
% record command to be run
req = sprintf('%s %s', cmd, args);

% shell out to a system command, e.g. jsoc_info
[shell_status,stdout] = system(req);
if shell_status ~= 0,
  json_content = sprintf('system call failed: %.80s', stdout);
  status = 0; % return convention: status = 0 means failure
  return
end;
status = 1; % OK

% check for and remove header
json_header = 'Content-type: application/json';
len_header = length(json_header);
% valid output always starts with json_header.  But, if something was put on
% stderr, the two streams will be mixed.  To avoid confusion in this
% case, we do not unconditionally strip the first length(json_header) chars.
% (remove it in a way that is efficient even for long result strings)
json_content = stdout;
if (length(stdout) >= len_header) &&  ...
      strcmp(stdout(1:len_header), json_header),
  % this idiom avoids errors if json_content == json_header
  json_content(1:len_header) = [];
end;
return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% shell_encode: escape shell metacharacters, esp. []
%
function s2 = shell_encode(s1)

s2 = regexprep(s1, '([\[\]{}])','\\$1');
return
end

