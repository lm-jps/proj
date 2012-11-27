function s = jsoc_cell2struct_keys(r)
%jsoc_cell2struct_keys	convert cell array returned by jsoc to struct array
% 
% s = jsoc_cell2struct_keys(r)
% * Given a return value r from a JSOC query, which is packaged as a
% cell array, convert it to a struct array indexed by keyword.
% This is a convenience function.
% 
% Inputs:
%   cell r{nk} of (name,values) pairs
% 
% Outputs:
%   struct s mapping names to values
%    
% See Also:

% Written by Michael Turmon (turmon@jpl.nasa.gov) on 25 Jun 09.
% Copyright (c) 2009.  All rights reserved.

% (Pass on error checking)

%
% Computation
% 
nk = length(r);
s = struct([]);
for k = 1:nk,
  name1 = r{k}.name;
  key1 = name2key(name1);
  % TODO: make this more robust/predictable by using more than first?
  type1 = deduce_type(r{k}.values{1});
  switch type1,
   case 'double',
    s(1).(key1) = str2double(r{k}.values);
   case 'time',
    block = char(r{k}.values); % convert to matrix to enable indexing
    s(1).(key1) = datenum(block(:,1:end-4),'yyyy.mm.dd_HH:MM:SS');
   case 'char',
    s(1).(key1) = r{k}.values;
   otherwise,
    warning('unknown type %s', type1);
  end;
end
return
end

% map JSOC metadata name (like FITS keyword name) into matlab struct key
function key = name2key(name)

key = lower(strrep(name, '-', '__'));
return;
end

% determine type from string
function typ = deduce_type(s)

if length(s) > 4 && strcmp(s(end-3:end), '_TAI'),
  typ = 'time';
elseif strcmp(s, 'MISSING') || ~isnan(str2double(s))
  typ = 'double';
else
  typ = 'char';
end
return;
end
  
