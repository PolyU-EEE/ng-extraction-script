clc;close all;clear
%% Input
file_path = 'Fig.4.xlsx';
rawdata1 = readtable(file_path, 'Sheet', '4(a)', 'VariableNamingRule', 'preserve');
lambda = table2array(rawdata1(:,1));        % Wavelength [nm]
Transmission = table2array(rawdata1(:,2));  % Transmission [dB]
%% Find peak parameters 
end_num = [2447 5928  6490  8440  11816 12895 15360  17574 18858 20440 21305  22045 22632 23143 23461 23485     23535 length(Transmission)];  
windows = [221   149   85   51    41    47    17     15    15     13    9     7     7     5     2      3        2      2];    
Min_dB  = [-7.8 -7.8  -7.5  -7    -6    -7    -7     -5.1  -4.76  -7.67 -7.67 -7.67 -6.79 -7    -7     -12.25   -11.9  -11.27]; 
Transmission_1 = zeros(size(Transmission));
peak_indices_all = [];
peak_heights_all = [];
start_index = 1;
%% Find peaks
for i = 1:length(end_num)
    if i == length(end_num)
        end_index = length(Transmission); 
    else
        end_index = end_num(i);
    end
    current_window = windows(i);
    half_window = floor(current_window / 2);
    expand_start = max(1, start_index - half_window);
    expand_end = min(length(Transmission), end_index + half_window);
    expanded_data = Transmission(expand_start:expand_end);
    sm_expanded = movmean(expanded_data, current_window, 'Endpoints', 'shrink');
    actual_start_in_expanded = start_index - expand_start + 1;
    actual_end_in_expanded = actual_start_in_expanded + (end_index - start_index);
    current_sm_data  = sm_expanded(actual_start_in_expanded:actual_end_in_expanded);
    current_index = start_index:end_index;
    Transmission_1(current_index) = current_sm_data ;
    if i ~=15
        [peak_heights, peak_locs] = findpeaks(-current_sm_data,"MinPeakDistance",2*windows(i),"MinPeakHeight",-Min_dB(i) );
    else
        [peak_heights, peak_locs] = findpeaks(-current_sm_data,"MinPeakDistance",3*windows(i),"MinPeakHeight",-Min_dB(i) );
    end 
    peak_heights = -peak_heights;
    global_peak_indices = peak_locs + start_index - 1;
    peak_indices_all = [peak_indices_all;global_peak_indices];
    peak_heights_all = [peak_heights_all;peak_heights ];
    start_index = end_index + 1;
end
%% ng value calculation
L = 1e6;        % L is 1 millimeter [nm]
ngref = 4.4;
lambda_peak = lambda(peak_indices_all);
lambda_1 = lambda_peak(1:end-1);
lambda_2 = lambda_peak(2:end);
Delta_lambda = diff(lambda_peak);
ng = lambda_1.*lambda_2./(Delta_lambda.*L) + ngref;  
end_num2 = [1539.38  1540.58 1556.76 1563.54 1567.56 1571.56 1572.956 1573   1573.12  1573.5]; 
windows2 = [   1      4      2        3      3       3      3         1      4        3];   
segment_indices = ones(1, length(end_num2) + 1);
segment_indices(1) = 1;
for i = 1:length(end_num2)
    [~, idx] = min(abs(lambda_1 - end_num2(i)));
    segment_indices(i + 1) = idx;
end
lambda_ng = cell(length(end_num2),2);  
for seg = 1:(length(segment_indices) - 1)
    start_idx = segment_indices(seg);
    end_idx = segment_indices(seg + 1);
    window_size = windows2(seg);
    seg_lambda = lambda_1(start_idx:end_idx);
    seg_ng = ng(start_idx:end_idx);
    if window_size == 1
        lambda_ng{seg,1} = seg_lambda;
        lambda_ng{seg,2} = seg_ng;
    else
        num_points = length(seg_lambda);
        num_full_windows = floor(num_points / window_size);
        Temp1 = seg_lambda(1:(num_full_windows * window_size)); 
        Temp_lam = mean(reshape(Temp1, window_size, num_full_windows),1);                     
        Temp1 = seg_ng(1:(num_full_windows * window_size));  
        Temp_ng = mean(reshape(Temp1, window_size, num_full_windows),1);             
        if num_full_windows * window_size < num_points
            remaining_start = num_full_windows * window_size +1;
            remaining_end = num_points;
            Temp_lam(end+1) = mean(seg_lambda(remaining_start:remaining_end));
            Temp_ng(end+1) = mean(seg_ng(remaining_start:remaining_end));
        end
        lambda_ng{seg,1} = Temp_lam(:);
        lambda_ng{seg,2} = Temp_ng(:);
    end
    if seg == (length(segment_indices) - 1)
        lambda_m = vertcat(lambda_ng{:,1});
        ng_m = vertcat(lambda_ng{:,2});
        lambda_m(end+1) = mean(lambda_1(377)); 
        ng_m(end+1) = mean(ng(377)); 
    end
end
%% Output
figure('Name',"Fig.4")
plot(lambda_m, ng_m,'LineWidth', 3,'Color','#4497C4')
hold on ;
scatter(lambda_m, ng_m, 80, ...
        'LineWidth', 3, ...
        'MarkerEdgeColor', '#4497C4', ...   
        'MarkerFaceColor', 'w');        
xlim([1475 1580]); 
xlabel('Wavelength (nm)'); 
ylabel('Measurement group index')
box on
