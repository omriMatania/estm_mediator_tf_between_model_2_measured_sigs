function [x_list] = convert_mtrx2list(x,dim)
% convert_mtrx2list convert matrix to list.
%   X_LIST = CONVERT_MTRX2LIST(X) converts the matrix X to the list
%   X_LIST.
%
%   CONVERT_MTRX2LIST(X,DIM) convert the matrix along the dimension DIM of 
%   X.

if nargin < 2
    dim = 1;
end % of if

if dim == 2
    x = x.' ;
end

size_x = size(x) ;
M = size_x(2) ; % number of examples

x_list = cell(1, M) ; % pre-allocation

for m = 1 : 1 : M
    x_list{m} = x(:, m) ;
end % of for

end % of conver_mtrx2list

