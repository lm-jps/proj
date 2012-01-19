% series_struct.m
%
% Usage  : series_struct series_name
% Example: series_struct hmi.lev0e  
%          a = series_struct ('hmi.lev0e', 'web_access')
%

function results = series_struct(series_name, web_access)

% Check series_name
if (nargin <1)
    fprintf ('Series name not specified.\n\n');
    return;
end

% Check web_access or local (default)
if ((nargin > 1) & (web_access == 'web_access'))
    web_access = true;
else
    web_access = false;   
end

% default to *not* printing the response
verbose = false;

if (web_access == true)

    try
        url_string = strcat('http://jsoc2.stanford.edu/cgi-bin/ajax/jsoc_info?op=series_struct&ds=',series_name);
        json_content = urlread_jsoc(url_string);
  
        % turmon: change API to handle old matlabs
        [results,err] = parse_json(json_content);
        if ~isempty(err),
           error('could not parse json reply: %s', err);
        end;
    catch
        disp(lasterror);
    return;
end

if (results.status > 0) % Note: status type is double
    fprintf ('Fail to get a response from JSOC\n');
    return;
end

else

    line = strcat('jsoc_info op=series_struct ds=',series_name);
    [status results] = system(line);
    if status > 0,
      error('system call to <%s> failed', line);
    end;
    json_content = regexprep(results,'Content-type: application/json',''); % remove this line

    % turmon: change API to handle old matlabs
    [results,err] = parse_json(json_content);
    if ~isempty(err),
       error('could not parse json reply: %s', err);
    end;

end



if verbose,
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
