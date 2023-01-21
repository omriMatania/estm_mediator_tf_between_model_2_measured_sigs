function [xs_ATF, tf_estm] = tf_step(xs_f, ys_f, K_model, L_model)
% tf_step   estimtes the tf with the lowest MSE. 
%   S = tf_step(XS,YS,K_MODEL,L_MODEL) finds the transfer function with the
%   lowest MSE between XS_F and YS_F under the contrains of at most 
%   K_MODEL + 1 poles and L_MODEL + 1 zeros.

tf_estm = estm_tf_min_mse(xs_f, ys_f, K_model, L_model);

xs_ATF = calc_x_ATF(xs_f, tf_estm);

end % of tf_step

