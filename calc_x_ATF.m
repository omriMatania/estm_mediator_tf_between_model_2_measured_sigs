function [xs_atf] = calc_x_ATF(xs, tf, domain)
% calc_x_ATF X after the transfer function.
%   S = calc_x_ATF(XS, TF) is the vectors of XS after the transfer
%   function TF.
%
%   calc_x_ATF(XS, TF, DOMAIN) calculates the vectors of XS after the 
%   transfer function in the domain DOMAIN. Available options are:
%
%   'frequency' - XS and XS_ATF are in the frequency domain. 
%   'time'      - XS and XS_ATF are in the time domain.
%   
%   TF is always in the frequency domain.

if nargin < 3
    domain = 'frequency';
end % of if

if strcmp(domain, 'time')
    xs = fft(xs,[],1); % DFT on the columns of xs
end % of if

size_x = size(xs);
M = size_x(2); % number of examples

tf_mtrx = repmat(tf, 1, M);
xs_atf = tf_mtrx .* xs;

if strcmp(domain, 'time')
    xs_atf = real(ifft(xs_atf,[],1)); % IDFT on the columns of xs_atf
end % of if

end % of calc_x_ATF

