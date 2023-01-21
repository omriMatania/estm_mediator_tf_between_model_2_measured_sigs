function [y] = convert_a_b_coeff_2_tf_f_fast(num_pnts, b, a)
% convert_a_b_coeff_2_tf_f_fast converts the coefficients a and b to a transfer funciton. 

x = zeros(num_pnts, 1) ;
x(1) = 1 ;
y_t = filter(b, a, x) ;
y = fft(y_t) ;

end % of convert_a_b_coeff_2_tf_f_fast

