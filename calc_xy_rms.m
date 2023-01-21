function xy_rms = calc_xy_rms(xs,ys,tf,domain)
% calc_xy_rms the rms between xs and ys. 
%   S = calc_xy_rms(XS,YS) is the mean rms between the normalized vectors
%   (x / rms(x)) and (y / rms(y)) in the metrices XS and YS, 
%   correspondingly.
%   
%   CALC_XY_RMS(XS,YS,TF) calculates the rms when the vectors in XS are
%   propagated through the transfer function TF.
%
%   S = CALC_XY_RMS(XS, YS, DOMAIN) specifies the domain of XS and YS. 
%   Available options are:
%
%   'frequency' - XS and YS are presented in the frequency domain. 
%   'time'      - XS and YS are presented in the time domain.
%
%   TF is always in the frequency domain.

if nargin == 2 || isempty(tf)
    
    size_xs = size(xs);
    N = size_xs(1);
    tf = ones(N, 1);

end % of if

if nargin <= 3
    domain = 'frequency';
end % of if

size_xs = size(xs);
N = size_xs(1); % number of samples per example
M = size_xs(2); % number of examples

if strcmp(domain, 'frequency')
    
    xs_atf = calc_x_ATF(xs, tf); % x_atf - x after transfer function. 
    xs_atf_rms = rms(xs_atf, 1); % rms of x_atf.
    xs_atf_rms = repmat(xs_atf_rms , N, 1);
    xs_atf = xs_atf ./ xs_atf_rms;
    ys_rms = repmat(rms(ys, 1), N, 1);
    ys = ys ./ ys_rms;
    xy_rms = mean(rms(xs_atf - ys, 1), 2);

elseif strcmp(domain, 'time')
    
    xy_rms = 0 ;
    
    for m = 1 : 1 : M
    
        xs_atf = real(ifft(tf .* fft(xs(:, m)))) ;
        xs_atf = xs_atf / rms(xs_atf, 1) ;
        ys = ys(:, m) ;
        ys = ys / rms(ys, 1) ;
        xy_rms = xy_rms + (1/M) * rms(ys - xs_atf, 1) ; 
    
    end % of for
    
end % of if

end % of calc_xy_rms