%{
These codes implement the new suggested algorithm in the following study on
the demonstration dataset:

Matania O., Zamir O., Bortman J. "A new tool for model examination: Estimation
of the mediator transfer function between the model and measured signals", 
Journal of Sound and Vibration, (2023). https://doi.org/10.1016/j.jsv.2023.117560

You should set the transfer function number (between 1-24) and the noise
amplitude and also the data path and processed data path. More information
is described in the ReadMe file at:
https://github.com/omriMatania/estm_mediator_tf_between_model_2_measured_sigs

For any question do not hesitate to contact Omri Matania via the email
omrimatania@gmail.com.

You can use the codes and the data for any academic purpose, however you
are requested to cite the study:

Matania O., Zamir O., Bortman J. "A new tool for model examination: Estimation
of the mediator transfer function between the model and measured signals", 
Journal of Sound and Vibration, (2023). https://doi.org/10.1016/j.jsv.2023.117560
%}

clear all ; close all ;
rng('default')

% Set transfer function number (1-24) and noise amplitude

tf_num = 1 ;
noise_amp = 0.25 ;

% Set the directory of the data and the processed data
% You can download them from:
% https://drive.google.com/drive/folders/17syVgYJh8yAM0Q_-ALHQIqCoAvY6kbvq?usp=sharing 
data_path = 'D:\data\estm_mediator_tf_between_model_2_measured_sigs\data\' ;
processed_data_path = 'D:\data\estm_mediator_tf_between_model_2_measured_sigs\processed_data\' ;

% ----------------------------------------------------------------------- %

% general parameters:
N = 1024 ;
num_winds = 100 ;
delays_ups_num_pnts = 10 ;
ups_initial_guess = 1 ;
max_itr = 20 ;
epsilon = 0 ;
min_decline_ratio = 0.99999 ;
p_test = 0.3 ;
p_val = 0.3 ;
num_of_types_of_sigs = 5 ;

[K_list, L_list] = set_K_and_L_lists(N) ; % number of poles and zeros

% Pre-allocation
tic()
selcted_K_mtrx = zeros(num_of_types_of_sigs, 1) ;
selcted_L_mtrx = zeros(num_of_types_of_sigs, 1) ;
NMSE_after_EMTF_list = zeros(num_of_types_of_sigs, 1) ;

% load the data
[xs_t, xs_f, ys_t, ys_f, tf_real, NMSE_internat_deviation] = ...
    load_dataset(data_path, processed_data_path, N, tf_num, noise_amp, num_winds) ;

% Decrease the number of signal by half for faster running time
[xs_t] = decrease_num_sigs(xs_t) ; [xs_f] = decrease_num_sigs(xs_f) ;
[ys_t] = decrease_num_sigs(ys_t) ; [ys_f] = decrease_num_sigs(ys_f) ;

% Calculate NMSE before EMTF
xs_test_delays = delays_step(xs_f, ys_f, delays_ups_num_pnts) ;
NMSE_before_EMTF = calc_xy_rms(xs_test_delays, ys_f) ;

% Calculate NMSE after EMTF
for out_ind = 1 : 1 : num_of_types_of_sigs
    
    [xs_train, xs_test, ys_train, ys_test] = ...
        train_test_split_differnet_sigs(xs_f, ys_f, num_of_types_of_sigs, out_ind) ;
    
    % Steps 1 - 7 of the algorithm (see Fig. 2).
    [K_min, L_min, K_min_L_min_init_delays] = estm_K_L(...
        num_of_types_of_sigs - 1, xs_train, ys_train, ups_initial_guess, ...
        K_list, L_list, max_itr, min_decline_ratio, delays_ups_num_pnts) ;

    % ----------------------------------------------------------------------- %
    % Step 7 of the algorithm: Repeat steps 3-6 on the training set with the chosen K and L (see Fig. 2).
    [best_tf_estm] = approximate_the_mediator_tf(xs_train, ys_train, K_min, L_min, max_itr, ...
        min_decline_ratio, delays_ups_num_pnts, K_min_L_min_init_delays) ;
    % ----------------------------------------------------------------------- %

    xs_test_ATF = calc_x_ATF(xs_test, best_tf_estm) ;
    xs_test_ATF_AD = delays_step(xs_test_ATF, ys_test, delays_ups_num_pnts) ;
    NMSE_after_EMTF = calc_xy_rms(xs_test_ATF_AD, ys_test) ;

    selcted_K_mtrx(out_ind) = K_min ;
    selcted_L_mtrx(out_ind) = L_min ;
    NMSE_after_EMTF_list(out_ind) = NMSE_after_EMTF ;

end % of for
NMSE_after_EMTF = mean(NMSE_after_EMTF_list(:, tf_num)) ;

disp(['NMSE internal deviation = ', num2str(NMSE_internat_deviation)])
disp(['NMSE before applying EMTF = ', num2str(NMSE_before_EMTF)])
disp(['NMSE after applying EMTF = ', num2str(NMSE_after_EMTF)])