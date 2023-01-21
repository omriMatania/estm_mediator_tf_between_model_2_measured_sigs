function [xs_training, xs_test, ys_training, ys_test, training_inds, test_inds] = ...
    train_test_split_differnet_sigs(xs_f, ys_f, num_of_types_of_sigs, out_ind)
% train_test_split_differnet_sigs splits xs_f and ys_f into training and
% test sets.

num_sigs_per_type = size(xs_f, 2) / num_of_types_of_sigs ;
            
test_inds = [num_sigs_per_type * (out_ind - 1) + 1 : 1 : num_sigs_per_type * out_ind] ;

training_inds = setdiff([1 : 1 : size(xs_f, 2)], test_inds) ;

xs_training = xs_f(:, training_inds) ;
xs_test = xs_f(:, test_inds) ;
ys_training = ys_f(:, training_inds) ;
ys_test = ys_f(:, test_inds) ;

end % of train_test_split_differnet_sigs

