%%source: https://uk.mathworks.com/matlabcentral/answers/1988-function-ecdf-break-down-for-large-datasets
function [v_f,v_x, ia,ic] = a_ecdf(v_data)
    nb_data = numel(v_data);
    [v_sorted_data, v_i] = sort(v_data);
    [v_unique_data, ia,ic] = unique(v_data);
    nb_unique_data = numel(v_unique_data);
    v_data_ecdf = zeros(1,nb_unique_data);
    for index = 1:nb_unique_data
        current_data = v_unique_data(index);
        v_data_ecdf(index) = sum(v_sorted_data <= current_data)/nb_data;
    end
    v_x = [v_unique_data(1) v_unique_data];
    v_f = [0 v_data_ecdf];
end