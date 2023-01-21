function [tf_estm] = approximate_the_mediator_tf(xs_f, ys_f, K_model, L_model, max_itr, ...
    min_decline_ratio, delays_num_pnts, K_L_init_delays)
% approximate_the_mediator_tf estimates the transfer function that minimizes
% the MSE between xs_f and ys_f under the constrains of K_model and L_model.
% There is not guarantee that the iterative greedy procedure will convergence 
% to a global minimum and hence the function just approximates the mediator
% transfer funciton.

xs_train = xs_f ;
ys_train = ys_f ;

% for each K and L the function does gready tf and dleayes steps iterations 

xs_train_AD = calc_x_delays(xs_train, K_L_init_delays) ;

% greedy tf and delays steps iterations.
tf_estm = estm_tf_and_delays_by_greedy_iterations(xs_train_AD, ...
    ys_train, K_model, L_model, max_itr, min_decline_ratio, ...
    delays_num_pnts) ;

end % of approximate_the_mediator_tf