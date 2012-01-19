function cm=prism2(n,bw)
%prism2	prism-like colormap for nominal data
% 
% cm=prism2(n)
% * This returns an n-length colormap of distinct colors, intended
% to be used to make plots of categorical or nominal data.
% * Black, gray, and white aren't among them by default; they can be
% put in by giving the second argument, which gives the grayscale 
% values of the first few colors.  For example, if bw = [1 0 0.5],
% the first three colors will be white, black, and 50% gray.
% 
% Inputs:
%   int n;
%   opt real bw(np) = [];
% 
% Outputs:
%   real cm(3,n)
% 
% See Also:  prism, colorcube

if nargin < 2, bw=[]; end;

cm0 = colorcube(27);
cm1 = cm0([1,14,3,7,4,13,6,11,17,20],:);
% add in grays
cm1 = [repmat(bw(:), [1 3]); cm1];
m = size(cm1, 1); % number of colors we have
cm = cm1(1 + mod(0:n-1,m),:);
return;

