function [K_list, L_list] = set_K_and_L_lists(N)
% set_K_and_L_lists generates two lists of number of poles and zeros.

ii = 1 ; K_list = [0] ;
while ii < N
    K_list = [K_list, ii] ; ii = 2 * ii ;
end % of while
if K_list(end) < N - 1
    K_list = [K_list, N - 1] ;
end % of if

L_list = K_list ;

end % set_K_and_L_lists