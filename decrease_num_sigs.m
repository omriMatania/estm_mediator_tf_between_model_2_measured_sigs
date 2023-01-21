function [new_xs] = decrease_num_sigs(xs)
% decrease_num_sigs decreases the number of signals in xs for fuster
% running time.

num_records_per_type = 6 ;
num_sigs_types = 5 ;
decreasing_ratio = 1 / 2 ;

num_sigs_per_type = size(xs, 2) / num_sigs_types ;
new_xs = zeros(size(xs, 1), decreasing_ratio * size(xs, 2)) ;

ind = 1 ;
for ii = 1 : 1 : num_sigs_types
    
    xs_type = xs(:, (ii-1) * num_sigs_per_type + 1 : ii * num_sigs_per_type) ;
    num_sigs_per_record = size(xs_type, 2) / num_records_per_type ;
    
    for jj = 1 : 1 : num_records_per_type
        
        xs_record = xs_type(:, (jj-1) * num_sigs_per_record + 1 : ...
            jj * num_sigs_per_record) ;
        xs_record = xs_record(:, 1 : decreasing_ratio * size(xs_record, 2)) ;
        new_xs(:, ind : ind + size(xs_record, 2) - 1) = xs_record ;
        ind = ind + size(xs_record, 2) ;
        
    end % of for
    
end % of for

end % of decrease_num_sigs