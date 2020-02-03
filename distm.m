%DISTM Compute square Euclidean distance matrix
% 
%   D = DISTM(A,B)
% 
% INPUT
%   A,B   matrices; B is optional, default B = A 
%
% OUTPUT
%   D     Square Euclidean distance matrix
%
% DESCRIPTION  
% Computation of the square Euclidean distance matrix D between two
% sets A and B. If A has M objects and B has N objects, then D is 
% [M x N].

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% Simplified by Cris Luengo to not depend on other PRTOOLS functions.
% See http://www.prtools.org/ for the full code.

function D = distm(A,B)

if nargin < 2
	B = A;
end
[ma,ka] = size(A); 
[mb,kb] = size(B);
if (ka ~= kb)
	error('Feature sizes should be equal')
end

% The order of operations below is good for the accuracy.
D = ones(ma,1)*sum(B'.*B',1);
D = D + sum(A.*A,2)*ones(1,mb);
D = D - 2 .* (+A)*(+B)';

J = find(D<0);                  % Check for a numerical inaccuracy. 
D(J) = zeros(size(J));          % D should be nonnegative.

if ((nargin < 2) && (ma == mb)) % take care of symmetric distance matrix
	D = (D + D')/2;              
	D([1:ma+1:ma*ma]) = zeros(1,ma);
end
