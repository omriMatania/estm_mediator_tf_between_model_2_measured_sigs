function tf_estm = estm_tf_min_mse_known_A_v(N, K_model, L_model, A, v)
% estm_tf_min_mse_known_A_v estimtes the transfer function with the lowest MSE. 

% solve linear equations
w = (A \ v) ;

% extracted a and b coefficients
b_estm = w(1 : K_model + 1) ;
a_estm = [1; w(K_model + 2 : (K_model + 1) + L_model)] ;

% conver a and b coefficients to tf in the frequency domain
tf_estm = convert_a_b_coeff_2_tf_f_fast(N, b_estm, a_estm) ;

end % of estm_tf_min_mse_known_A_v

