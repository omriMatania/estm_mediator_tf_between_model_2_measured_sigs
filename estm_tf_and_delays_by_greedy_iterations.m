function [tf_estm, delays_estm, errs_mse] = ...
    estm_tf_and_delays_by_greedy_iterations(xs_f, ys_f, K_model, ...
    L_model, max_itr, min_decline_ratio, delays_num_pnts, ys_train_interp)
% estm_tf_and_delays_by_greedy_iterations estimates the transfer function and the
% delays between xs and ys using greedy iterations. 
%   [TF_ESTM, DELAYS_ESTM, ERRS_MSE] = estm_tf_and_delays_by_greedy_iterations(XS,
%   YS,MAX_ITR) estimates the transfer function and the delayes between XS
%   and YS using greedy estimation. In each delay step, the function finds 
%   the optimal delays (minimum of the MSE) between XS and YS. In each 
%   transfer function step, the function finds the optimal transfer 
%   function (minimum MSE) between XS and YS.
%   XS and YS should be given in the frequency domain. ERRS_MSE are the
%   errors along the iterations.
%   MAX_ITR bounds the maximal number of iterations by MAX_ITR. Each 
%   iteration comtains a delays step and a transfer function step.
%
%   ESTM_TF_AND_DELAYS_BY_ITERATIONS(... , K, L) set the order of the poles
%   coefficients to be K and the order of the zeros coefficients to be L.
%
%   ESTM_TF_AND_DELAYS_BY_ITERATIONS(... , MIN_DECLINE_RATIO) bounds the
%   the minimal improvment between following iterations. If the MSE of the
%   r step divided by the MSE of the r-1 step is less than 
%   MIN_DECLINE_RATIO then the function stops to iterate.
%
%   ESTM_TF_AND_DELAYS_BY_ITERATIONS(... , INIT_DELAYS) beggins the 
%   iterations with INIT_DELAYS as the delays.
%
%   ESTM_TF_AND_DELAYS_BY_ITERATIONS(... , INIT_TF) beggins the iterations
%   with INIT_TF as the trasfer function. If INIT_DELAYS and INIT_TF are
%   both given, then the first step to follow is delays step.  

if nargin < 9
    ys_train_interp = ys_f ; 
end % of if

size_x = size(xs_f) ;
N = size_x(1) ; % number of samples per vector.
M = size_x(2) ; % number of examples.

errs_mse = [calc_xy_mse(xs_f, ys_f)] ;

delays_estm = zeros(M, 1);
decline_ratio = 0 ;
num_itr = 0 ;
xs_AD = xs_f ;

while max_itr > num_itr && min_decline_ratio > decline_ratio  
    
    % ----------------------------------------------------------------------- %
    % Step 4 of the algorithm: Transfer function step (see Fig. 2).
    [xs_ATF, new_tf_estm] = tf_step(xs_AD, ys_f, K_model, L_model) ;
    % ----------------------------------------------------------------------- %

    % ----------------------------------------------------------------------- %
    % Step 5 of the algorithm: Delays step (see Fig. 2).
    [xs_ATF_AD, new_delays_estm] = delays_step(xs_ATF, ys_train_interp, delays_num_pnts) ;
    % ----------------------------------------------------------------------- %
    
    delays_estm = mod(delays_estm + new_delays_estm, N * delays_num_pnts) ;
    tf_estm = new_tf_estm ;
    
    xs_AD = calc_x_delays(xs_f, delays_estm, delays_num_pnts) ;
    
    errs_mse = [errs_mse; calc_xy_mse(xs_ATF_AD, ys_f)] ;
    decline_ratio = errs_mse(end) / errs_mse(end-1) ;
    num_itr = num_itr + 1 ;
    
end % of while

end % of estm_tf_and_delays_by_iterations