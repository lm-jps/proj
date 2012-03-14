function [res,msg] = rs_list(query, method)
% rs_list	return result from jsoc_info op=rs_list
%
% [res,msg] = rs_list(query, method)
% * Invoke jsoc_info, op=rs_list, using the given query.  This can 
% return lists of keywords, or segments, or links, or any combination
% thereof.
% * If method contains 'web' (the default) the CGI gateway speaking JSON 
% over http will be used.  This works remotely.  If method is 'system', 
% the system jsoc_info command will be run and its output intercepted.
% * If method contains 'verbose' (not the default), a printed summary
% of the result will be shown.
% * This routine queries jsoc2.stanford.edu, and re-tries twice on 
% failure.
% * If the second output, msg, is not requested, errors will result if
% there is no response from JSOC.  If msg is requested, it will be 
% empty if no error, or a descriptive string if there was an error, but
% an error will not be raised for communication or parsing failures.
%
% Examples: 
%   res = rs_list('hmi.M_720s[$]&key=T_OBS')
%   res = rs_list('hmi.M_720s[$]&key=T_OBS','web');
%   [res,msg] = rs_list('hmi.M_720s[2010.06.01/1d]&key=T_OBS','web');
% 
% Inputs:
%   string query
%   opt string method = 'web'
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
if nargin < 1 || nargin > 2,
  error('Need 1 or 2 input arguments');
end
% input defaults
if nargin < 2, method = 'web'; end;
% parameter defaults
% abs() = between-try delay in s, [] for no retries, <0 for announcement
%  (set to announce for long delays)
retries = [0.1 10 -60 -60 -60]; 
host = 'jsoc2.stanford.edu';


%
% Computation
% 

% Check for web or local (default)
use_ajax = ~isempty(strfind(method, 'web'));
verbose = ~isempty(strfind(method, 'verbose'));

% allow error() exit if true
do_err_out = (nargout < 2);
res = [];

if use_ajax,
  % url escape, replacing spaces with &
  query2 = regexprep(query, ' ', '&');
  query_url = sprintf('http://%s/cgi-bin/ajax/jsoc_info?op=rs_list&ds=%s',...
                      host, query2);
  max_try = length(retries)+1;
  for ntry = 1:max_try,
    [json_content,status] = urlread_jsoc(query_url);
    if status == 1, 
      break; % success
    end;
    % retry mechanism: a bit baroque
    if ntry < max_try,
      % more tries remain
      msg = 'Retrying after pause.';
      delay = retries(ntry);
    else,
      % the last try failed
      msg = 'Giving up.';
      delay = -0.001; % tiny delay, with announcement
    end;
    if delay < 0,
      fprintf('%s: Query %d of %d failed.  %s\n', mfilename, ntry, max_try, msg);
    end;
    pause(abs(delay)); % can have delay < 0
    % (go on to retry as permitted)
  end;
  if status == 0,
    % failed
    msg = sprintf('urlread_jsoc failed (%s)', query_url);
    res = json_content;
    if do_err_out, error(msg); else, return; end;
  end;

else,
  % shell out to a system command
  line = strcat('jsoc_info op=rs_list ds=', query);
  % escape shell metacharacters, esp. []
  line_e = regexprep(line, '([\[\]{}])','\\$1');
  [status,res] = system(line_e);
  if status ~= 0,
    msg = sprintf('jsoc_info system call (%s) failed: %s', line_e, res);
    if do_err_out, error(msg); else, return; end;
  end;
  % a valid message always starts with this...
  % remove it in a way that is efficient even for long result strings.
  json_header = 'Content-type: application/json';
  if length(res) >= length(json_header) &&  ...
        strcmp(res(1:length(json_header)), json_header) == 1,
    json_content = res(length(json_header)+1:end);
  else,
    json_content = json_header;
  end;
end

% parse the JSON response
[res,msg1] = parse_json(json_content);
if ~isempty(msg1),
  fprintf('rs_list: failed parse.\n  query = %s\n  reply length = \n', ...
          query, length(json_content));
  fprintf('rs_list: reply follows:\n%s', json_content);
  % failed
  msg = sprintf('could not parse json reply (%s)', msg1);
  if do_err_out, error(msg); else, return; end;
end;
% Note: status type is double
if res.status > 0, 
  msg = sprintf('Error response from JSOC (status = %g)', res.status);
  if do_err_out, error(msg); else, return; end;
end

%
% Display the response...this is getting into the weeds
%
if verbose,
  % Keywords {'name': , 'type': , 'units': , 'note': }
  fprintf('\nKeywords:\n');
  for k = 1:length(res.keywords),
    fprintf ('  %s\n',res.keywords{k}.name);  
  end
  % Records found
  fprintf('\nTotal records found: %d\n',res.count);
  for j = 1:res.count,
    for k = 1:length(res.keywords),
      if ischar(res.keywords{k}.values{j})
        fprintf('%s\t',res.keywords{k}.values{j});
      else
        fprintf('%f\t',res.keywords{k}.values{j});            
      end
    end
    fprintf('\n');
  end
  % Segments
  fprintf('\nSegments:\n');
  for k=1:length(res.segments)    
    fprintf ('  %s\n',res.segments{k}.name);
    disp(res.segments{k})
  end
  % Links
  fprintf('\nLinks:\n');
  for k=1:length(res.links)    
    fprintf ('  %s\n',res.link{k}.name);
  end
end; % if verbose

msg = ''; % signals OK

return

