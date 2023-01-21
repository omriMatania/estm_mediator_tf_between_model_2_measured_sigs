function [A_12_long, v_11_long] = generate_long_mtrx_vctr(xs_f, ys_f, ...
    init_delays)
% generate_long_mtrx_vctr generates long matrix A_12_long and vector v_11_long. 

size_x = size(xs_f) ;
N = size_x(1) ; % number of samples in x.
M = size_x(2) ; % number of exsamples in x.

xs_AD = calc_x_delays(xs_f, init_delays) ;

delays_mtrx = repmat(ones(M, 1).', N, 1) ;
num_mtrx = repmat([0 : 1 : N - 1].' / N, 1, M) ;
exp_mtrx_m1 = exp(-2 * pi * 1i * delays_mtrx .* num_mtrx) ;

[A_N_m1, v_N_m1] = calc_A_v(N-1, N-1, xs_AD, ys_f) ;
[A_N_m1_d_m1, ~] = calc_A_v(N-1, N-1, xs_AD .* exp_mtrx_m1, ys_f) ;

A_12_long = [A_N_m1(1 : N, N + 1 : end), A_N_m1_d_m1(1 : N, N + 1 : end)] ;
v_11_long =  [v_N_m1(1 : N); v_N_m1(1 : N)] ;

end % of generate_long_mtrx_vctr