function dists=euclideanDist(X,C)
% Vectorized euclidean dist between X and C, where X and C are coordinates
% of points, each point per row.

%Debug
%{
X=[0 0 0; 0 0 1; 1 2 3];
C=X;
%}
% Check inputs
if size(X,2)<3
    X=[X zeros(size(X,1),1)];
end
% Self distance
if nargin<2
    C=X;
end
% Calculate the sum of squares for all input vectors and
% for all cluster centers / models.
%
% Matrix dimensions:
%   X  [m  x  n]
%  XX  [m  x  1]
%   C  [k  x  n]
%  CC  [1  x  k]
XX = sum(X.^2, 2);
CC = sum(C.^2, 2)';

% Calculate the dot product between each input vector and
% each cluster center / model.
%
% Matrix dimensions:
%   X  [m  x  n]
%   C  [k  x  n]
%   C' [n  x  k]
%  XC  [m  x  k]
XC = X * C';

% Calculate the Euclidean distance between all input vectors in X
% and all clusters / models in C using the following equation:
%
%   z = sqrt(||x||^2 - 2xc' + ||c||^2)
%
%  Step 1: Subtract the column vector XX from every column of XC.
%  Step 2: Add the row vector CC to every row of XC.
%
% Matrix dimensions:
%     XX  [m  x  1]
%     XC  [m  x  k]
%     CC  [1  x  k]
%  dists  [m  x  k]
%
dists = real(sqrt(bsxfun(@plus, CC, bsxfun(@minus, XX, 2*XC))));

end