clc; clear all

% Load data
filename = 'NewYork_Shif_TypSum_MPC_2004_TypOcc_TypBehav_NoTES_03032024_163856';
load([filename '.mat']);

% Hardware data 
hdata_name = fieldnames(HardwareData);
hdata = eval(['HardwareData.' hdata_name{1}]);

% Save the hdata table to a CSV file
csvFileName = [hdata_name{1} '.csv']; % Construct CSV file name
writetable(hdata, csvFileName); % Save table to CSV