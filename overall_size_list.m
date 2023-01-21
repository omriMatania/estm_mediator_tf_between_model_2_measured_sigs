function [M] = overall_size_list(x_list)
% overall_size_list calculates the sum of the vectors of x_list.

M = 0 ;
for ii = 1 : 1 : length(x_list)
    x = x_list{ii} ;
    size_x = size(x) ;
    M_ii = size_x(2) ;
    M = M + M_ii ;
end % of for

end % overall_size_list

