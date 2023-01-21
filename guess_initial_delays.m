function [init_delays] = guess_initial_delays(xs_f, ys_f, ups_initial_guess)
% guess_initial_delays guesses initial delays between xs_f and ys_f using a 
% recursive procedure implemented in estm_tf_delays_recursivly.  
N = size(xs_f, 1) ;

x_list = convert_mtrx2list(xs_f) ;

[~, init_delays] = estm_tf_delays_recursivly(x_list, ys_f, ups_initial_guess) ;

init_delays = mod(init_delays, N) ;

end % of guess_init_delays

