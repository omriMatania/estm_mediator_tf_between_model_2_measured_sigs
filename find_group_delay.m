function [best_tf_estm, best_delayes_estm, xs_AD] = find_group_delay( ...
    xs_f, ys_f, init_delays, K_model, L_model, A_12_long, v_11_long)
% find_group_delay finds delay of the group (XS_F)
%   [BEST_TF_ESTM, BEST_DELAY_ESTM, X_AD_ATF] = FIND_GROUP_DELAY(XS_F, 
%   YS_F, INIT_DELAYS, K_MODEL, L_MODEL) are the best transfer function and
%   delay between XS_F and YS_F. X_AD_ATF is XS_F after the estimated delay 
%   and after the transfer funciton.
%
%   INIT_DELAYS is the initial delyas between XS_F and YS_F.
%
%   K_MODEl and L_MODEl are the number of coefficients in the denominator 
%   and numerator of the transfer function after Z transform, respectivaly.
%   (tf = b(s) / a(s)).

size_x = size(xs_f) ;
N = size_x(1) ; % number of samples in x.
M = size_x(2) ; % number of exsamples in x.

min_mse = inf ;
mse_vctr = zeros(N, 1) ;

xs_AD = calc_x_delays(xs_f, init_delays) ;

delays_mtrx = repmat(ones(M, 1).', N, 1) ;
num_mtrx = repmat([0 : 1 : N - 1].' / N, 1, M) ;
exp_mtrx = exp(2 * pi * 1i * delays_mtrx .* num_mtrx) ;

[A_delay_0, v_delay_0] = calc_A_v(K_model, L_model, xs_AD, ys_f) ;
A_11 = A_delay_0(1 : K_model + 1, 1 : K_model + 1) ;
A_22 = A_delay_0(K_model + 2 : end, K_model + 2 : end) ;
v_21 = v_delay_0(K_model + 2 : end) ;

A_delay = zeros((K_model + 1) + L_model, (K_model + 1) + L_model) ;
A_delay(1 : K_model + 1, 1 : K_model + 1) = A_11 ;
A_delay(K_model + 2 : end, K_model + 2 : end) = A_22 ;

for delay = 0 : 1 : N - 1
    
    A_delay_12 = A_12_long(1 : K_model + 1, 1 + delay : L_model + delay) ;
    
    A_delay(1 : K_model + 1, K_model + 2 : end) = A_delay_12 ;
    A_delay(K_model + 2 : end, 1 : K_model + 1) = A_delay_12.' ;
    
    v_delay_11 = v_11_long(N + 1 - delay : N + K_model + 1 - delay) ;
    v_delay = [v_delay_11; v_21] ;

    tf_estm = estm_tf_min_mse_known_A_v(N, K_model, L_model, ...
        A_delay, v_delay) ;
    
    tf_mtrx = repmat(tf_estm, 1, M) ;
    xs_ATF = tf_mtrx .* xs_AD ;
    
    diff = xs_ATF - ys_f ;
    diff_abs_squre = (diff) .* conj(diff) ;
    sum_1 = sum(diff_abs_squre, 1) ;
    mse = sum(sum_1, 2) ; 

    mse_vctr(delay + 1) = mse ;
    
    if min_mse > mse 
        
        min_mse = mse ;
        best_tf_estm = tf_estm ;
        best_delayes_estm = delay * ones(M, 1) + init_delays ;
    
    end % of if
    
    xs_AD = xs_AD .* exp_mtrx ;

end % of for

xs_AD = calc_x_delays(xs_f, best_delayes_estm) ;

end % of find_group_delay