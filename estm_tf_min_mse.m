function tf_estm = estm_tf_min_mse(xs_f, ys_f, K_model, L_model)
% estm_tf_min_mse estimtes the tf with the lowest MSE. 
%   S = estm_tf_min_mse(XS,YS,K_MODEL,L_MODEL) finds the transfer function with the
%   lowest MSE between XS_F and YS_F under the constrains of K_MODEL + 1 poles and
%   L_MODEL + 1 zeros.

size_x = size(xs_f) ;
N = size_x(1) ; % number of samples

% solve linear equations
[A, v] = calc_A_v(K_model, L_model, xs_f, ys_f) ;
w = (A \ v) ;

% extracted a and b coefficients
b_estm = w(1 : K_model + 1) ;
a_estm = [1; w(K_model + 2 : (K_model + 1) + L_model)] ;

% conver a and b coefficients to tf in the frequency domain
tf_estm = convert_a_b_coeff_2_tf_f_fast(N, b_estm, a_estm) ;

end % of estm_tf_min_mse