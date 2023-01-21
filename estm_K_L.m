function [K_min, L_min, K_min_L_min_init_delays] = estm_K_L(...
    num_of_types_of_sigs, xs_f, ys_f, ups_initial_guess, K_vctr, L_vctr, max_itr, ...
    min_decline_ratio, delays_num_pnts)
% estm_K_L estimates K and L (corresponding to the number of zeros and poles, respectively)
% corresponding to the minimal MSE on the validation set after applying EMTF on the training set.

% ----------------------------------------------------------------------- %
% Step 1 of the algorithm: Calculate initial delays guess (see Fig. 2).
init_delays = guess_initial_delays(xs_f, ys_f, ups_initial_guess) ;
% ----------------------------------------------------------------------- %

[A_12_long, v_11_long] = generate_long_mtrx_vctr(xs_f, ys_f, init_delays) ;
min_mse_val = inf ;
ys_f_interp = Nyquist_interp(ys_f, size(ys_f, 1) * delays_num_pnts, 'frequency') ;
for L_ind = 1 : 1 : length(L_vctr)
    L_model = L_vctr(L_ind) ; 

    % for each K and L the function do gready tf and dleayes steps iterations
    for K_ind = 1 : 1 : length(K_vctr)
        
        K_model = K_vctr(K_ind) ;
        mse_vals = zeros(num_of_types_of_sigs, 1) ;
        
        % ----------------------------------------------------------------------- %
        % Step 3 of the algorithm: Calculate modified group delay (see Fig. 2).
        % In this version of the code Step 3 is before Step 2 for faster
        % runnint time. There is a little effect on the results.
        [~, K_L_init_delays, xs_f_AD] = find_group_delay(xs_f, ys_f, ...
            init_delays, K_model, L_model, A_12_long, v_11_long) ;
        % ----------------------------------------------------------------------- %
        
        for out_ind = 1 : 1 : num_of_types_of_sigs
            
            % ----------------------------------------------------------------------- %
            % Step 2 of the algorithm: Dividing into training and validation sets (see Fig. 2).
            [xs_train, xs_val, ys_train, ys_val, train_inds, ~] = ...
                train_test_split_differnet_sigs(xs_f_AD, ys_f, ...
                num_of_types_of_sigs, out_ind) ;
            % ----------------------------------------------------------------------- %
            
            ys_train_interp = ys_f_interp(:, train_inds) ;

            % greedy tf and delays steps iterations. (Steps 4 + 5)
            tf_estm = estm_tf_and_delays_by_greedy_iterations(xs_train, ...
                ys_train, K_model, L_model, max_itr, min_decline_ratio, ...
                delays_num_pnts, ys_train_interp) ;

            % calculated the MSE on the validation set
            xs_val_ATF = calc_x_ATF(xs_val, tf_estm) ;
            [xs_val_ATF_AD, ~] = delays_step(xs_val_ATF, ys_val, delays_num_pnts) ;
            
            % ----------------------------------------------------------------------- %
            % Step 6 of the algorithm: Calculate MSE in validation set (see Fig. 2).
            mse_val = calc_xy_mse(xs_val_ATF_AD, ys_val) ;
            mse_vals(out_ind) = mse_val ;
            % ----------------------------------------------------------------------- %

        end % of for
        
        % ----------------------------------------------------------------------- %
        % Step 7 of the algorithm: Choose the K and L degree with the lowest MSE (see Fig. 2).
        if min_mse_val > mean(mse_vals)
            min_mse_val = mean(mse_vals) ;
            K_min = K_model ; 
            L_min = L_model ;
            K_min_L_min_init_delays = K_L_init_delays ;
        end % of if
        % ----------------------------------------------------------------------- %

    end % of for

end % of for

end % of estm_K_L