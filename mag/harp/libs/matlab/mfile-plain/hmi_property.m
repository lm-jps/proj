function res=hmi_property(mode, key, val)
%hmi_property	set and get HMI properties
% 
% res=hmi_property(mode, key, val)
% * Operates as a property setter, or a property getter.
% Supports three types of call.  The ordinary setter and getters:
%    hmi_property('set', propertyname, val)
%    res = hmi_property('get', propertyname)
% plus the overall reset-to-defaults of all properties:
%    hmi_property('reset')
% * For the getter/setter, you can use val = 'default' to get or 
% set the property to the original default value.
% * The function behaves as if 'reset' was called once at the beginning
% of a session.
% * Currently supports:
%     key = 'nrt_mode'     --  val = 0 or 1
%                              default = 0
%     key = 'jsoc_method'  --  val = web, shell, or server
%                              default = web
%     key = 'jsoc_address' --  val = '(default)'  (the literal string)
%                              This value allows jsoc_info_request to choose
%                              the correct service.  For fine-grained tweaking,
%                              set up your own value.
%     key = 'at_stanford'  --  val = 0 or 1
%                              default determined by inspecting hostname
% * Note, keys have to be usable as Matlab `struct' field names.
%    
% Inputs:
%   string mode = 'reset', 'get', or 'set'
%   string key
%   (anything) val
% 
% Outputs:
%   (anything) res
% 
% See Also:

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 17 Feb 2011
% Copyright (c) 2011.  All rights reserved.

% preserve properties across invocations -- initially equals []
persistent hmi_props

if isempty(getenv('HOST')),
  % put this back in the environment so we don't re-execute it
  [st,hn] = system('hostname');
  if st == 0,
    hn(end) = []; % strip newline
  else,
    hn = 'localhost';
  end;
  setenv('HOST', hn);
end;
  
at_stanford = ~isempty(regexpi(getenv('HOST'), 'stanford'));

%% Set up default property values
%
hmi_default = struct();
% can toggle nrt_mode externally to use NRT metadata
hmi_default.nrt_mode = 0;
% strategy with jsoc access is to put an easy top-level switch 
% here (values: shell/web/remote), and to put the specific 
% access points into jsoc_info_request.  Using "default" as the
% host value selects that method.  However, you can still use a
% custom value for the host through this routine.
hmi_default.at_stanford = at_stanford;
if at_stanford,
  % (turns out shell is slower, even at Stanford, so use web)
  % hmi_default.jsoc_method = 'shell';
  hmi_default.jsoc_method = 'web';
  hmi_default.jsoc_address = '(default)';
else
  hmi_default.jsoc_method = 'web';
  hmi_default.jsoc_address = '(default)';
end;
% retry list: 
%    one retry is made for each entry,
%    [] for no retries, 
%    abs() = between-try delay in s, 
%    < 0 for announcement (now, is set to announce for long delays)
hmi_default.jsoc_retries = [0.1 10 -60 -60 -60]; 

% if not ever initialized, take the default
if isempty(hmi_props),
  hmi_props = hmi_default;
end;

% 
% Error checking
% 
if all(nargin  ~= [1 2 3]), error ('Bad input arg number'); end;
if all(nargout ~= [0 1]), error ('Bad output arg number'); end;

%
% Computation
% 
res = []; % ensure something is set
%% Allow to get or set specific properties
%
if strcmp(mode, 'reset'),
  %% reset-all
  hmi_props = hmi_default;
elseif strcmp(mode, 'get'),
  %% getter
  if ~isfield(hmi_props, key),
    error('Property %s has not been set', key);
  elseif nargin > 2 && strcmp(val, 'default'),
    % return the default setting (not the current setting)
    res = hmi_default.(key);
  else,
    % return the current setting
    res = hmi_props.(key);
  end;
elseif strcmp(mode, 'set'),
  %% setter
  if ~isfield(hmi_props, key),
    error('Property %s was not initialized', key);
  elseif strcmp(val, 'default'),
    % use the default
    hmi_props.(key) = hmi_default.(key);
  else,
    hmi_props.(key) = val;
  end;
else,
  error('Mode must be get or set');
end;
return;
end

