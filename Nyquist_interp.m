function [xs_interp] = Nyquist_interp(xs, new_num_pnts, domain)
% Nyquist_interp implements Nyquist interpolation of XS. 
%   S = Nyquist_interp(XS, NEW_NUM_PNTS) interpolates XS to contain 
%   NEW_NUM_PNTS by Nyquist theorem.
%   
%   Nyquist_interp(XS, NEW_NUM_PNTS, DOMAIN) gets as input and outs as 
%   output XS and XS_INTERP in the domain DOMAIN.
%
%   Nyquist_interp(XS, NEW_NUM_PNTS, DOMAIN) specifies the domain of XS 
%   and XS_INTERP. Available options are:
%
%   'frequency' - XS and XS_INTERP are presented in the frequency domain. 
%   'time'      - XS and XS_INTERP are presented in the time domain.

if nargin < 3
    domain = 'time' ;
end % of if

size_x = size(xs) ;
N = size_x(1) ;
M = size_x(2) ;

if strcmp(domain, 'time')
    xs_f = fft(xs) ;
elseif strcmp(domain, 'frequency')
    xs_f = xs ;
end % of if

if new_num_pnts == N 
    xs_interp_f = xs_f ;
elseif mod(N, 2) == 0
    xs_interp_f = (new_num_pnts / N) * [xs_f(1 : N / 2, :) ; 0.5 * xs_f(N / 2 + 1, :) ; ...
        zeros(new_num_pnts - N - 1, M) ; 0.5 * conj(xs_f(N / 2 + 1, :)) ; ...
        xs_f(N / 2 + 2 : end, :)] ;
elseif mod(N, 2) == 1
    xs_interp_f = (new_num_pnts / N) * [xs_f(1 : floor(N / 2) + 1, :) ; ...
        zeros(new_num_pnts - N, M) ; xs_f(floor(N / 2) + 2 : end, :)] ;
end % of if

if strcmp(domain, 'time')
    xs_interp = real(ifft(xs_interp_f)) ;
elseif strcmp(domain, 'frequency')
    xs_interp = xs_interp_f ;
end % of if

end % of Nyquist_interp

