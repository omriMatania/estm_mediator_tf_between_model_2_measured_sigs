function [tf_estm, delays_estm, x_new] = estm_tf_delays_recursivly(...
    x_list, ys_f, ups_initial_guess)
% estm_tf_delays_recursivly finds initial guess for the delays recursivly.
%   [TF_ESTM, DELAYS_ESTM, X_NEW] = estm_tf_delays_recursivly(X_LIST, YS_F,
%   PRINT_FLAG) are the estimated transfer function and the delyas between
%   X_list and YS_F.
%
%   PRINT_FLAG indicates if to printing the things along the function.
%
%   TF_ESTM is always in the frequency domain.

if length(x_list) == 1
    
    tf_estm = [] ;
    delays_estm = zeros(overall_size_list(x_list), 1) ;
    x_new = x_list{1} ;
    
    return
    
end % of if

% divide X and Y to two groups
len_x_list = length(x_list) ;
x_list_left = x_list(1 : ceil(len_x_list / 2)) ;
x_list_right = x_list(ceil(len_x_list / 2) + 1 : end) ;
y_left = ys_f(:, 1 : overall_size_list(x_list_left)) ;
y_right = ys_f(:, overall_size_list(x_list_left) + 1 : end) ;

% calculate recursivly the delays between X and Y
[~, delays_left, x_new_left] = estm_tf_delays_recursivly(x_list_left, y_left, ...
    ups_initial_guess) ;

[~, delays_right, x_new_right] = estm_tf_delays_recursivly(x_list_right, y_right, ...
    ups_initial_guess) ;

% combine the two groups
[tf_estm, groups_delay, x_new] = combine_groups(x_new_left, x_new_right, ...
    ys_f, ups_initial_guess) ;

delays_estm = [delays_left ; delays_right] + groups_delay ;

end % of estm_tf_delays_recursivly

%   -------------------------------------------------------------------   %

function [best_tf_estm, best_groups_delay_estm, x_new] = combine_groups(x_left, ...
    x_right, y, ups_initial_guess)
% combine_groups combines x_left and x_right

size_x1 = size(x_left) ;
N = size_x1(1) ;
M1 = size_x1(2) ;
size_x2 = size(x_right) ;
M2 = size_x2(2) ;

min_mse = inf ;

delta_delay = 1 / ups_initial_guess ;

for delay = 0 : delta_delay : N - (1 / delta_delay)
    
    x_rigth_delay = calc_x_delays(x_right, delay * ones(M2, 1)) ;
    x = [x_left, x_rigth_delay] ;

    [x_ATF, y, tf_estm] = rectiy_tf_h(x, y) ;
    
    mse = calc_xy_mse(x_ATF, y, ones(N, 1), 'f') ;

    if min_mse > mse 
        
        min_mse = mse ;
        best_tf_estm = tf_estm ;
        best_delay = delay ;
        
    elseif isnan(mse)
    
    end % of if

end % of for

best_groups_delay_estm = [zeros(M1, 1) ; best_delay * ones(M2, 1)] ;

x_rigth_delay = calc_x_delays(x_right, best_delay * ones(M2, 1)) ;
x_new = [x_left, x_rigth_delay] ;

end % of combine_groups