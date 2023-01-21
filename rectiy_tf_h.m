function [x_ATF, ys_f, tf_estm] = rectiy_tf_h(xs_f, ys_f)
% rectiy_tf_h between X and Y.
%   [X_ATF, YS_F, TF_ESTM] = rectiy_tf_h(XS_F, YS_F) are X and Y after the
%   estimation of the transfer function TF_ESTM that minimzes the MSE
%   between the X and Y. The hypothesis class is all the transfer
%   functions with N coordinates, where N is the number of smples in each
%   vector.
%
%   TF_ESTM is always in the frequency domain.

size_x = size(xs_f) ;
M = size_x(2) ;

xx = sum(xs_f .* conj(xs_f), 2) ;
xy = sum(xs_f .* conj(ys_f), 2) ;

% estimate the transfer function
h_real = real(xy) ./ xx ;
h_img = imag(conj(xy)) ./ xx ;
tf_estm = h_real + 1i * h_img ;

% calc x after the transfer function
tf_mtrx = repmat(tf_estm, 1, M) ;
x_ATF = tf_mtrx .* xs_f ;

end % of rectiy_tf_h