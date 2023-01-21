function [xs_AD] = calc_x_delays(xs_f, delays, num_pnts)
% calc_x_delays XS_F after DELAYS 
%   XS_AD = CALC_X_DELAYS(XS_F,DELAYS) is XS_F after the DELAYS.

size_x = size(xs_f) ;
N = size_x(1) ; % number of samples.
M = size_x(2) ; % number of exsamples.

if nargin < 3
    num_pnts = 1 ; 
end % of if

if num_pnts > 1
    
    delays_mtrx = repmat(delays.', N * num_pnts, 1) ;
    num_mtrx = repmat([0 : 1 : (N * num_pnts) - 1].' / (N * num_pnts), 1, M) ;

    xs_f = Nyquist_interp(xs_f, N * num_pnts, 'frequency') ;
    
    xs_AD = xs_f .* exp(2 * pi * 1i * delays_mtrx .* num_mtrx) ;

    xs_AD_old = xs_AD ;
    xs_AD = zeros(N, M) ;
    
    xs_t = ifft(xs_AD_old) ;

    for m = 1 : 1 : M
        x_t = xs_t(:, m) ;
        x_t = x_t([1 : num_pnts : N * num_pnts]) ;
        xs_AD(:, m) = fft(x_t) ;
    end % of for
    
else
    delays_mtrx = repmat(delays.', N, 1) ;
    num_mtrx = repmat([0 : 1 : N - 1].' / N, 1, M) ;
    xs_AD = xs_f .* exp(2 * pi * 1i * delays_mtrx .* num_mtrx) ;
end % of if

end % of calc_x_delays