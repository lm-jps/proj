function [im,err] = jsoc_imload(series, trec, seg)
%jsoc_imload	load a jsoc image segment at a given trec from series
% 
% [im,err] = jsoc_imload(series, trec, seg)
% * Load a jsoc image from segment `seg' within the given series, at 
% the given trec.  For instance:
%    [im, err] = jsoc_imload('hmi.M_720s', '2011.01.01_12:36_TAI', 'magnetogram')
% * As a shortcut, you can combine either the trec or the seg, or both, into
% the series, like this:
%    [im, err] = jsoc_imload('hmi.M_720s[2011.01.01_12:36_TAI]{magnetogram}')
%    [im, err] = jsoc_imload('hmi.M_720s[2011.01.01_12:36_TAI]', [], 'magnetogram')
%    [im, err] = jsoc_imload('hmi.M_720s{magnetogram}', '2011.01.01_12:36_TAI')
% * Also, if seg is not given, it is taken to be '***ALL***', so the above can
% also be given as:
%    [im, err] = jsoc_imload('hmi.M_720s[2011.01.01_12:36_TAI]')
% * If there was no error, err will be empty.  
% * If there was an error, im will be empty and err will be an explanatory string.  
% Note that no matlab error will be raised, so that the caller can handle the 
% error as desired without a try/catch block.
% * It should be that if the err output is not asked for, a matlab error will
% be raised anyway.
% * We should be able to return an image stack if a multi-image request is given,
% e.g. for trec's matching multiple images, but this is not implemented.  
% It would be easy to put this in, perhaps not in full generality.
% 
% Inputs:
%   string series   -- data series
%   opt string trec
%   opt string seg
% 
% Outputs:
%   real im(...)
%   string err
% 
% See Also: 

% 
% Error checking
% 
if nargin > 3, error('Bad input arg number'); end
% if all(nargout ~= [2 3 4]), error ('Bad output arg number'); end  
if nargin < 3 || isempty(seg),
  % this will be over-ridden by a segment embedded in the series
  seg = '***ALL***';
end;
if nargin < 2,
  % in this case it must be given in the series
  trec = '';
end;

%
% Computation
% 
im = [];

% FIXME: juggle multi-args correctly
% FIXME: allow multi-segment load


% Parse the data series string., 
% We allow it to encode T_REC and the segment name 
% if it is given as: series[TREC]{segname}
% Either the [] part or the {} part are optional, and can be 
% broken out as separate arguments
% The re is: ^ + (series) + ([time] | nil) + ({segment} | nil) $
% where each of the three subgroups is captured as a token.
% If there was no match, there was a problem with the series name
% as it was given.  If time or segment came out as empty, they
% must be non-empty in the argument list.
% Besides alphanumerics, in series we allow the literal .
% Besides digits, in time we allow ., :, and _TAI.
%   We also allow @, / and dhms for time intervals (/10d@1h, etc.)
% In segment, we allow alphanumerics plus * as in **ALL**.
% We may not be allowing for all cases here.
re = '^([\w.]+)(\[[\d:._TAI/@dhms]+\]|)(\{[\w*]+\}|)$';
m = regexp(series, re, 'tokens');
if isempty(m),
  err = 'Badly formed data series input';
  return;
end;
% if the "trec" part of the match is a digit and not a time,
% dump it back into the series (indexing for HARPs)
if ~isempty(m{1}{2}),
  % if it is just [ + digits + ]
  if regexp(m{1}{2}, '^\[\d+\]$'),
    m{1}{1} = [m{1}{1} m{1}{2}];
    m{1}{2} = '';
  end;
end;
% because there was a match, the series must be non-empty
series1 = m{1}{1};
% get the trec
if ~isempty(m{1}{2}),
  trec1 = m{1}{2}(2:end-1); % strip off []
elseif ~isempty(trec),
  trec1 = trec;
else,
  err = 'Could not get T_REC'
  return;
end;
% get the segment
if ~isempty(m{1}{3}),
  seg1 = m{1}{3}(2:end-1); % strip {}
elseif ~isempty(seg),
  seg1 = seg;
else,
  err = 'Could not get seg';
  return;
end;

obj1 = sprintf('%s[%s]', series1, trec1);
% query JSOC to locate the segment
a = rs_list(obj1, {'key', 'T_REC', 'seg', seg1});

%% the T_REC indexes as a cell array of strings:
% val_trec_s = a.keywords{1}.values;
%% times as doubles a la matlab
% temp = jsoc_cell2struct_keys(a.keywords);
% val_trec = temp.t_rec'; 

if a.status > 0,
  err = sprintf('Failed to load %s, not found by rs_list\n', obj1);
  return;
end;
% exhaustively test the response
if ~isfield(a, 'segments') || ~iscell(a.segments) || length(a.segments) ~= 1,
  err = sprintf('Failed to load %s, no segments\n', obj1);
  return;
end;
s = a.segments{1};
if ~isfield(s, 'values') || ~iscell(s.values) || length(s.values) ~= 1,
  err = sprintf('Failed to load %s, no segment values\n', obj1);
  return;
end;
v = s.values{1};
if ~ischar(v) || isempty(v),
  err = sprintf('Failed to load %s, segment values bad\n', obj1);
  return;
end;
if v(1) ~= '/',
  % 2/2011, offline gives v = 'NoDataDirectory'
  err = sprintf('Failed to load %s, segment seems offline\n', obj1);
  return;
end;
% this is the SUMS filename
file1 = v;
at_stanford = ~isempty(regexpi(getenv('HOST'),'stanford'));
if at_stanford,
  % load as a plain file
  [im,err1] = loadfitsa(file1);
else,
  % load thru http, with optional cache (loadfitsa understands URLs)
  file1 = sprintf('http://jsoc.stanford.edu/%s', file1);
  [im,err1] = loadfitsa_jsoc(file1);
end;
% Check for an error
if err1,
  % summarize
  err = sprintf('Failed to load %s from %s\n', obj1, file1);
  return;
end;

err = '';
return;

end
 
