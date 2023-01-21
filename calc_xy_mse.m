function xy_mse = calc_xy_mse(xs_f, ys_f, tf, domain)
% calc_xy_mse the MSE between xs and ys. 
%   S = calc_xy_mse(XS,YS) is the MSE between the metrices xs and ys.
%   
%   CALC_XY_RMS(X,Y, TF) calculates the MSE when the vectors in X are
%   propagated through the transfer function TF.
%
%   S = CALC_XY_RMS(X, DOMAIN) specifies the domain of XS and YS. Available
%   options are:
%
%   'frequency' - XS and YS are presented in the frequency domain. 
%   'time'      - XS and YS are presented in the time domain.
%
%   TF is always in the frequency domain.

if nargin == 2 || isempty(tf)
    
    size_xs = size(xs_f);
    N = size_xs(1);
    tf = ones(N, 1);

end % of if

if nargin <= 3
    domain = 'frequency';
end % of if

if strcmp(domain, 'time')
    
    xs_f = fft(xs_f, [], 1); % DFT on the columns of xs.
    ys_f = fft(ys_f, [], 1);
    
end  % of if

size_xs = size(xs_f);
N = size_xs(1); % number of samples per example.
M = size_xs(2); % number of examples.
    
tf_mtrx = repmat(tf, 1, M) ;
xs_ATF = tf_mtrx .* xs_f ;
xy_mse = (1/(M*N)) * sum(sum(abs(xs_ATF - ys_f) .^ 2, 1), 2) ;

end % of calc_xy_mse