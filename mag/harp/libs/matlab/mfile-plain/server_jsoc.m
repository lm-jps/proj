function [json_msg,status,req] = server_jsoc(params)
%server_jsoc	query drms server for drms metadata
%
% [json_msg,status,req] = server_jsoc(params)
% * Return the results of calling the query_engine server with
% the given parameters encoded as a cell array of strings: 
%    params = {name1, value1, name2, value2,...}.
% Example names are 'ds' (data series) and 'op' (operation).
% * The JSON response string is json_msg, which must be decoded
% by the caller.  The request sent to the server is req.
% * Status string is empty for success, contains a message if error.
% * A server is started if needed, with its outputs to a log
% file supplied by this routine.  Be sure to stop the server
% before exit with either:
%   >> server_jsoc({'op', 'exit'});
% or, from the shell:
%   % kill <server_pid>
% * Queries can be made locally against a remote server.  Start
% a query_engine (manually) on n02.stanford.edu; say it ends up on
% port 48442.  Then, from the local system, create a tunnel:
%  % ssh -qTnN -L 41000:n02.stanford.edu:48442 solarport.stanford.edu
% You can point queries at port 41000 on the local machine, and they
% will go through the tunnel to run on n02.  
% * Another idiom combines the tunnel and the server, with:
%  % ssh -qTn -L 42000:n02.stanford.edu:54000 n02.stanford.edu "query_engine port=54000" &
% creating the tunnel and the server at once.  This relies on the
% proxy through solarport being set up in .ssh/config.
%
% Examples: 
%   server_jsoc({'op', 'rs_summary', 'ds', 'hmi.M_720s'})
%   server_jsoc({'ds', 'hmi.M_720s[$]', 'key', 'T_REC,T_OBS'})
% The second example omits the op because rs_list is the default.
% 
% Inputs:
%   cell params{1,2*NP} of string
% 
% Outputs:
%   string json_msg
%   string status -- empty if OK
%   opt string req
%    
% See Also:
%   The JSOC wiki page for jsoc_info, at:
%     http://jsoc.stanford.edu/jsocwiki/AjaxJsocConnect
%   The query_engine DRMS module documentation

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 25 Sep 2013.
% Copyright (c) 2013.  All rights reserved.

persistent port err_fn;
if isempty(port) port = 0; end;
% (leaving err_fn empty is OK)

% imports for socket services
import com.mathworks.mlwidgets.io.InterruptibleStreamCopier;
import java.net.Socket;
import java.io.*;

% 
% Error checking
% 
% (error out if called incorrectly)
if ~iscell(params) || rem(length(params), 2) ~= 0,
  error('params must be a cell array of even length');
end;
if ~all(cellfun(@ischar, params)),
  error('params must be a cell array of strings');
end;  
% could set this as an input parameter
host = '';
if isempty(host),
  host = 'localhost';
end;
% if host is local, we can start one ourselves
host_is_local = strcmp(host, 'localhost');
  

% default values for outputs
status = 'Unknown error';
json_msg = '';

% cons argument list, which is also an output
req = '';
sep = '';
for i = 1:2:length(params),
  req = [req sep params{i} '=' params{i+1}];
  sep = ' ';
end
% must be newline-terminated
%  (otherwise, write on socket seems to not complete)
req(end+1) = sprintf('\n'); 
clear sep;

% pseudo-param "port"
if length(params) > 0 && strcmp(params{1}, 'port'),
  if ~isempty(params{2}),
    port = str2double(params{2});
  end;
  fprintf('%s: Using port = %d\n', mfilename, port);
  return;
end;

% start local server if possible and needed
if host_is_local && port == 0,
  key = sprintf('%s-%06d-%03d', ...
                getenv('LOGNAME'), floor(1e6*rem(now,1)), fix(1000*rand));
  % FIXME: logfile, stderr, stdout
  port_fn = sprintf('/tmp/query-%s.port', key);
  log_fn  = sprintf('/tmp/query-%s.log',  key);
  err_fn  = sprintf('/tmp/query-%s.err',  key); % persistent!
  % start a server -- let user know.
  fprintf('%s: starting query server\n', mfilename);
  system(sprintf('rm -f %s %s', log_fn, err_fn));
  % could reword to use sh explicitly for this (via sh -c)
  % otherwise, can get another SHELL from env; MATLAB_SHELL and SHELL
  % in environment seem to sometimes not work (?)
  ccmd = 'query_engine';
  cmd = sprintf('%s portfile=%s logfile=%s >& %s &', ...
                ccmd, port_fn, log_fn, err_fn);
  % fprintf('Running: %s\n', cmd);
  % [s,r]=system('echo server_jsoc shell is $SHELL')
  [s,r] = system(cmd);
  if s ~= 0,
    % FIXME: better error handling
    fprintf('query engine error (%s): check %s\n', r, err_fn);
  else,
    % engine is in background, so wait until it creates the portfile
    %  (10 tries x 0.1 base ~= 5 sec total possible wait time)
    for attempt = 1:10,
      pause(0.1*attempt);
      if exist(port_fn, 'file'),
        port = load(port_fn);
        delete(port_fn);
        break;
      end;
    end;
    % (will have port = 0 if the file never appeared)
  end;
  if port ~= 0,
    fprintf('%s: query server ok, log in %s\n', mfilename, log_fn);
  end;
end
if port == 0,
  fprintf('%s: Could not start query engine\n', mfilename);
  status = 'Could not start query engine';
  return;
end;

% fprintf('%s: req = <%s>\n', mfilename, req);

% indicate failure
reply = '';
input_socket = [];
try,
  % connect to host:port -- throws if unable to connect
  input_socket = Socket(host, port);
  % fprintf('Connected to server\n');

  % stream for outgoing request
  printStream = java.io.PrintStream(input_socket.getOutputStream);
  % stream for incoming reply
  inputStream   = input_socket.getInputStream;
  % put the request out
  printStream.print(req);

  % read data from the socket
  %  bytes_available = inputStream.available;
  % (this is taken from urlread_jsoc)
  % (this StreamCopier is unsupported and may change at any time)
  byteArrayOutputStream = java.io.ByteArrayOutputStream;
  isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
  isc.copyStream(inputStream, byteArrayOutputStream);
  printStream.close;
  inputStream.close;
  byteArrayOutputStream.close;
  reply = native2unicode(typecast(byteArrayOutputStream.toByteArray','uint8'),'UTF-8');
    
  input_socket.close;
  status = ''; % OK
catch,
  % perhaps make this catch more restrictive?
  if ~isempty(input_socket)
    input_socket.close;
  end
  s = lasterror;
  fprintf('%s: socket exception: %s\n', mfilename, s.message);
  status = sprintf('Socket error: %.20s', s.message); % not OK
end;

% if we ran the exit operator, or experienced an error 
% (such as a socket error), record the port as not connected
if ~isempty(status) || length(strfind(req, 'op=exit')) > 0,
  port = 0;
  if ~isempty(err_fn),
    d = dir(err_fn);
    if length(d) == 1 && d.bytes == 0,
      delete(err_fn); % if empty, delete
    end;
  end;
end;
    
% check for and remove header
json_header = 'Content-type: application/json';
len_header = length(json_header);
% valid output always starts with json_header.  But, if something was put on
% stderr, the two streams will be mixed.  To avoid confusion in this
% case, we do not unconditionally strip the first length(json_header) chars.
% (remove it in a way that is efficient even for long result strings)
json_msg = reply;
if (length(reply) >= len_header) &&  ...
      strcmp(reply(1:len_header), json_header),
  % this idiom avoids errors if json_content == json_header
  json_msg(1:len_header) = [];
end;
return
end


