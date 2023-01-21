function [A, v] = calc_A_v(K, L, xs, ys)
% calc_A_v calculates A and v for the linear equitions of Aw = v.
%   [A, v] = CALC_A_v(K, L, XS, YS) is the metrix A and the vector v of 
%   the linear system Aw = v.

size_x = size(xs) ;
N = size_x(1) ; % number of samples in each vector.

A = zeros((K + 1) + L, (K + 1) + L) ; % pre-allocation.

xx = sum(xs .* conj(xs), 2) ;
xy = sum(xs .* conj(ys), 2) ;
yx = sum(ys .* conj(xs), 2) ;
yy = sum(ys .* conj(ys), 2) ;

xx_f = fft(xx) ; % DFT of xx
xy_f = fft(xy) ;
yx_f = fft(yx) ;
yy_f = fft(yy) ;

for alpha = 0 : 1 : K
    
    a_clmn = [xx_f(end - alpha + 1 : min(N, N - alpha + 1 + K)) ; ...
        xx_f(1 : K + 1 - alpha)] ;
    A(1 : K + 1, alpha + 1) = a_clmn ;
    a_clmn = - 1 * [yx_f(end - alpha + 2 : min(N, N - alpha + 2 + (L - 1))) ; ...
        yx_f(max(2 - alpha, 1) : L + 1 - alpha)] ;
    A(K + 2 : end, alpha + 1) = a_clmn ;
    
end % of for

for alpha = 1 : 1 : L
    
    a_clmn = -1 * [xy_f(end - alpha + 1 : min(N, N - alpha + 1 + K)) ; ...
        xy_f(1 : K + 1 - alpha)] ;
    A(1 : K + 1, (K + 1) + alpha) = a_clmn ;
    a_clmn = [yy_f(end - alpha + 2 : min(N, N - alpha + 2 + (L - 1))) ; ...
        yy_f(max(1 - alpha, 1) : L + 1 - alpha)] ;
    A(K + 2 : end, (K + 1) + alpha) = a_clmn ;
    
end % of for

A = real(A) ;
v = real([xy_f(1 : K + 1); -yy_f(2 : L + 1)]) ;

end % of calc_A_v