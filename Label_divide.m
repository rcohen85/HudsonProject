% load txt file
data = readtable('E:\S1099Hokkaido01_S01_RH403_WAV\SPICE_detector\clusters_below80K_3\cc10\S1099_below80K_Odontoceti_LOSpecific_labels.txt');

% save file
 saveDir = 'E:\S1099Hokkaido01_S01_RH403_WAV\SPICE_detector\clusters_below80K_3\cc10\';

% Specify the number of the column in which the information to be used for division is written
target_column = 3;

% Get a list of values in a column
unique_values = unique(data{:, target_column});

% Create and save a txt file for each value
for i = 1:height(unique_values)
    value = unique_values(i);
    
    % set filtering condition
    % filter = (data{:, target_column} == value); % if the value is double
    filter = strcmp(data{:, target_column}, value);
    
    % Extract rows that match the filtering condition
    filtered_data = data(filter, :);
       
    % Create split file name
    % file_name = sprintf('cclabel_cluster%s_new.txt', num2str(value)); % if the value is double
    file_name = ['cc2label_cluster', char(value), '_new.txt'];
    %file delimiter(区切り文字) is 'commma' or ','
    writetable(filtered_data, fullfile(saveDir,file_name), 'delimiter', 'comma');    
end
