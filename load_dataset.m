function [xs_t, xs_f, ys_t, ys_f, tf, NMSE_internat_deviation] = ...
    load_dataset(data_path, processed_data_path, N, tf_num, noise_amp, num_winds)
% load_dataset loads the demonstration dataset.

[xs_t, xs_f] = load_sigs_from_experiment(data_path, processed_data_path, N, num_winds) ;

% generate the mesaured signals after the transfer function.
[ys_t, ys_f, ys_f_group_1, ys_f_group_2, tf] = genrete_Y(data_path, ...
    tf_num, xs_t, noise_amp) ;

NMSE_internat_deviation = calc_xy_rms(ys_f_group_1, ys_f_group_2) ;

end % of load_dataset

%   -------------------------------------------------------------------   %
function [xs_t, xs_f] = load_sigs_from_experiment(data_path, processed_data_path, N, num_winds)
% load_sigs_from_experiment loads the signals from the experiments

sig_load = 10 ;
sig_rps = 40 ;
T = (38 / 17) * 1 / sig_rps ;
sig_num = 50*floor(60/T) ;

path = [data_path,'experimental signals_load_',num2str(sig_load), ...
    '_speed_',num2str(sig_rps),'\'] ;

xs_t = zeros(N, sig_num) ; % pre-allocation
xs_f = zeros(N, sig_num) ;

% the settings of the data
data_settings = [num2str(N),'_',num2str(num_winds),'_',num2str(sig_rps), ...
    '_', num2str(sig_load)] ;

try % if the datda with the specific settings was generated befire
    load([processed_data_path,'xs_t_',data_settings,'.mat'])
    load([processed_data_path,'xs_f_',data_settings,'.mat'])
catch
    counter = 1 ;
    for ii = 1 : 1 : 30
        load([path,'signal number ',num2str(ii),'.mat'])
        N_pnts = N * num_winds;
        if N ~= sig_properties.RFs
            acc_sig_cyc = resample(acc_sig_cyc, N, sig_properties.RFs, 10);
        end % of if
        N_segments = floor(length(acc_sig_cyc) / N_pnts) ;
        for mm = 1 : 1 : N_segments
            acc_sig = acc_sig_cyc(((mm-1)*N_pnts+1):(mm*N_pnts)) ;
            [SA_sig, ~] = calc_SA(acc_sig, N, 1, num_winds) ;
            xs_t(:, counter) = SA_sig ;
            xs_f(:, counter) = fft(SA_sig) ;
            counter = counter + 1 ;        
        end % of for
    end % of for
    
    xs_t = xs_t(:, 1 : (counter-1)) ;
    xs_f = xs_f(:, 1 : (counter-1)) ;
    
    save([processed_data_path,'xs_t_',data_settings,'.mat'], 'xs_t')
    save([processed_data_path,'xs_f_',data_settings,'.mat'], 'xs_f')

end % of try

end % of load_sigs_from_experiment

%   -------------------------------------------------------------------   %

function [ys_t, ys_f, ys_f_group_1, ys_f_group_2, tf] = genrete_Y(...
    data_path, tf_num, xs_t, noise_amp)
% genrete_Y generates Y based on xs_t and the transfer function.

size_x = size(xs_t) ;
N = size_x(1) ;
M = size_x(2) ;

% pre-allocation
ys_t = zeros(N, M) ;
ys_f = zeros(N, M) ;
ys_f_group_1 = zeros(N, M) ;
ys_f_group_2 = zeros(N, M) ;

% load the transfer function
tfs_path = [data_path, 'measured transfer functions\H', num2str(tf_num)] ;
load(tfs_path) ;
tf = interp_tf(H.', N) ;
tf = tf / rms(tf, 1) ;

% generate the delasy between X and Y
delays = rand(M, 1) ;

rms_x_1 = rms(xs_t(:, 1), 1) ;

for m = 1 : 1 : M
    xm_t = xs_t(:, m) ;
    num_pnts = 100 ;
    xm_t_long = Nyquist_interp(xm_t, N * num_pnts) ;
    v_t_long = circshift(xm_t_long, floor(delays(m) * N * num_pnts)) ;
    v_t = v_t_long([1 : num_pnts : N * num_pnts]) ;

    v_t = real(ifft(tf .* fft(v_t))) ;
    ys_t(:, m) = v_t + noise_amp * rms_x_1 * randn(size(xm_t)) ;
    ys_f(:, m) = fft(ys_t(:, m)) ;
    
    ys_f_group_1(:, m) = fft(v_t + noise_amp * rms_x_1 * randn(size(xm_t))) ;
    ys_f_group_2(:, m) = fft(v_t + noise_amp * rms_x_1 * randn(size(xm_t))) ;

end % of for

end % of genrete_Y