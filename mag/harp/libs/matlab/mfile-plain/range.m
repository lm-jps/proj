function m=range(A,B,C)
%range	min and max of a vector or matrix
% 
% m=range(A,B,C)
% * If B not supplied, finds the range of A as follows.  If A a vector, 
% returns a column vector of the min (in m(1,:)) and max (in m(2,:)) of
% the vector.  If A a matrix, returns the per-column mins and per-column
% maxes in a 2-row matrix.
% * If B supplied, restricts A to the range [B,C].  Elements below
% or above the interval are pulled to its edge.  Here B,C can be 
% scalars or sized just as A is.
% * The action of this function is intended to mirror max(), min().
% 
% Inputs:
%   real A[];
%   opt real B = []; 
%   opt real C = [];  -- B <= C
% 
% Outputs:
%   real m[];
% 
% See Also: range2, max, min

% 
% Error checking
% 
if all(nargin  ~= [1 3]), error ('Bad input arg number'); end
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end  

%
% Computation
% 
if nargin < 2,
  m = [min(A) ; max(A) ];
else
  m = max(min(A,C),B);
end;
return;
