% rs_summary.m
%
% Usage  : rs_summary series_name
% Example: rs_summary ('hmi.lev0e')
%          results = rs_summary('hmi.lev0e', 'web_access')
%
% Note: if (2nd parameter == 'web_access'), go through web server
%       default call jsoc_info directly

function results = rs_summary(series_name, web_access)

% Check series_name
if (nargin <1)
    fprintf ('Series name not specified.\n\n');
    return;
end

% Check web_access or local (defaul)
if ((nargin > 1) & (web_access == 'web_access'))
    web_access = true;
else
    web_access = false;   
end


if (web_access == true)
    
    try 
        line = strcat('http://jsoc2.stanford.edu/cgi-bin/ajax/jsoc_info?op=rs_summary&ds=',series_name);
        json_content = urlread_jsoc(line);  
    catch
        disp(lasterror);
        fprintf('Fail to get a response from JSOC\n');
        return;      
    end
    
else
    
    line = strcat('jsoc_info op=rs_summary ds=',series_name);
    % escape shell metacharacters, esp. []
    line_e = regexprep(line, '([\[\]{}])','\\$1');
    [status results] = system(line_e);
    if status > 0,
        error('system call <%s> failed', line_e);
    end;
    % a valid message always starts with this...
    % remove it in a way that is efficient even for long result strings.
    json_header = 'Content-type: application/json';
    if length(results) >= length(json_header) &&  ...
        strcmp(results(1:length(json_header)), json_header) == 1,
        json_content = results(length(json_header)+1:end);
    else,
        json_content = json_header;
    end;
end

% turmon: change API to handle old matlabs
[results,err] = parse_json(json_content);
if ~isempty(err),
   error('could not parse json reply: %s', err);
end;

if nargout == 0, 
  disp(results);
end;

return
