function [xs_delays, estm_delays] = delays_step(xs_f, ys_f, delays_ups_num_pnts)
% delays_step ESTM_DELAYS corresponding to minimum MSE
%   [XS_DELAYS, ESTM_DELAYS] = delays_step(XS_F, YS_F) is a greedy delays 
%   step claculting the ESTM_DELAYS corresponding to minimum MSE between 
%   XS_F and YS_F.
%
%   DELAYS_UPS_NUM_PNTS specifies the up sampling number op points in the
%   delay step process.

size_x = size(xs_f) ;
N = size_x(1) ;

if nargin < 3
    delays_ups_num_pnts = 1 ;
end % of if

xs_f = Nyquist_interp(xs_f, N * delays_ups_num_pnts, 'frequency') ;

if size(ys_f, 1) < size(xs_f, 1)
    ys_f = Nyquist_interp(ys_f, N * delays_ups_num_pnts, 'frequency') ;
end % of if

mses_as_function_of_d = -real(ifft(xs_f .* conj(ys_f))) ;
[~, ds_estm] = min(mses_as_function_of_d) ;
estm_delays = (ds_estm - 1).';

xs_delays = calc_x_delays(xs_f, estm_delays) ;
 
xs_t = ifft(xs_delays) ;
xs_t = xs_t([1 : delays_ups_num_pnts : N * delays_ups_num_pnts], :) ; 

xs_delays = fft(xs_t) ;

end % of delays_step
