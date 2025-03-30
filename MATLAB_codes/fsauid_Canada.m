% Clear the workspace and close all figures
clear;
close all;

% List of provinces and their CSV file names
provinces = {'AB', 'BC', 'MB', 'NB', 'NL', 'NS', 'NT', 'NU', 'ON', 'PE', 'QC', 'SK', 'YT'};
num_provinces = length(provinces);

% Initialize a variable to store all exposure data
all_exposure_data = [];

% Loop through each province and concatenate the exposure data
for i = 1:num_provinces
    % Get the province code
    province_code = provinces{i};
    
    % Construct the file name for the current province
    filename = [province_code, '.csv'];
    
    % Read the data from the CSV file into a table
    try
        exposure_data = readtable(filename);
    catch
        warning(['Could not read file: ', filename]);
        continue;
    end
    
    % Concatenate the data from each province
    all_exposure_data = [all_exposure_data; exposure_data];
end

% Ensure we have data
if isempty(all_exposure_data)
    error('No data loaded from CSV files.');
end

% Find the unique fsauid values
unique_fsauid = unique(all_exposure_data.fsauid);

% Initialize variables to store the results
num_unique_fsauid = length(unique_fsauid);
average_lon = zeros(num_unique_fsauid, 1);
average_lat = zeros(num_unique_fsauid, 1);
sum_number = zeros(num_unique_fsauid, 1);

% Loop through each unique fsauid
for i = 1:num_unique_fsauid
    % Get the current fsauid value as a string
    current_fsauid = unique_fsauid{i};
    
    % Find the rows corresponding to the current fsauid using strcmp
    indices = strcmp(all_exposure_data.fsauid, current_fsauid);
    
    % Extract the relevant rows
    fsauid_data = all_exposure_data(indices, :);
    
    % Calculate the average longitude and latitude for the current fsauid
    average_lon(i) = mean(fsauid_data.lon);
    average_lat(i) = mean(fsauid_data.lat);
    
    % Sum the 'number' column for the current fsauid
    sum_number(i) = sum(fsauid_data.number);
end

% Create a new table to store the results
fsauid_summary = table(unique_fsauid, average_lon, average_lat, sum_number, ...
    'VariableNames', {'fsauid', 'AverageLon', 'AverageLat', 'SumNumber'});

% Save the results to a CSV file
writetable(fsauid_summary, 'fsauid_summary.csv');

%% Find closest fsauid to certain location
% Earth's radius in kilometers
R = 6371;

% Input latitude and longitude (replace these with your values)
input_lat = 49.2827;  % Example: Latitude of Vancouver
input_lon = -123.1207; % Example: Longitude of Vancouver

% Initialize variables to store the minimum distance and corresponding fsauid
min_distance = inf;
closest_fsauid = '';
closest_rows = [];  % To store the rows corresponding to the closest fsauid

% Loop through each unique fsauid and calculate the distance
for i = 1:num_unique_fsauid
    % Get the latitude and longitude of the current fsauid
    lat = average_lat(i);
    lon = average_lon(i);
    
    % Convert latitude and longitude from degrees to radians
    lat1 = deg2rad(input_lat);
    lon1 = deg2rad(input_lon);
    lat2 = deg2rad(lat);
    lon2 = deg2rad(lon);
    
    % Calculate the difference in latitudes and longitudes
    dlat = lat2 - lat1;
    dlon = lon2 - lon1;
    
    % Haversine formula to calculate the distance
    a = sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    distance = R * c; % Distance in kilometers
    
    % Check if this is the minimum distance
    if distance < min_distance
        min_distance = distance;
        closest_fsauid = unique_fsauid{i};
        
        % Store the rows corresponding to this closest fsauid
        closest_rows = find(strcmp(fsauid_summary.fsauid, closest_fsauid));
    end
end

% Display the result
disp(['The closest fsauid to the given location is: ', closest_fsauid]);
disp(['The distance to the closest fsauid is: ', num2str(min_distance), ' km']);
disp('The rows corresponding to this fsauid are: ');
disp(fsauid_summary(closest_rows, :)); % Display the rows corresponding to the closest fsauid

disp('Analysis complete. The results have been saved to fsauid_summary.csv.');
