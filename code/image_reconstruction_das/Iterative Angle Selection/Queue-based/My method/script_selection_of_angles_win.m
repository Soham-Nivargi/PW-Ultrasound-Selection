%script to select angles iteratively
clear all;
close all;
clc;

addpath(genpath('../../../src'));


%-- Parameters
acquisition_type = 1;       %-- 1 = simulation || 2 = experiments
phantom_type = 2;           %-- 1 = resolution & distorsion || 2 = contrast & speckle quality
data_type = 1;              %-- 1 = IQ || 2 = RF
flag_display = 0;

%-- Parsing parameter choices
switch acquisition_type    
    case 1
        acquisition = 'simulation';
        acqui = 'simu';
        flag_simu = 1;
    case 2
        acquisition = 'experiments';
        acqui = 'expe';
        flag_simu = 0;
    otherwise       %-- Do deal with bad values
        acquisition = 'simulation';
        acqui = 'simu';
        flag_simu = 1;
end
switch phantom_type    
    case 1
        phantom = 'resolution_distorsion';
    case 2
        phantom = 'contrast_speckle';
    otherwise       %-- Do deal with bad values
        phantom = 'resolution';
end
switch data_type    
    case 1
        data = 'iq';
    case 2
        data = 'rf';
    otherwise       %-- Do deal with bad values
        data = 'iq';        
end

%-- Create path to load corresponding files
path_dataset = ['../../../../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_dataset_',data,'.hdf5'];
path_scan = ['../../../../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_scan.hdf5'];
path_pht = ['../../../../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_phantom.hdf5'];

%-- Read the corresponding dataset and the region where to reconstruct the image
dataset = us_dataset();
dataset.read_file(path_dataset);
scan = linear_scan();
scan.read_file(path_scan);
pht = us_phantom();
pht.read_file(path_pht);

pw_indices = (1+dataset.firings)/2;%-- dataset.firings corresponding to the total number of emitted steered plane waves
flag_simu = num2str(flag_simu);
flag_display = num2str(flag_display);
beamform_curr = das_iq_one_image(scan, dataset, pw_indices);

[~, beamform_2] = das_iq_for_selection_window_og(scan,pht,dataset,{[1,2], [74,75]});

beamform_1l = das_iq_one_image(scan,dataset,1);
beamform_1r = das_iq_one_image(scan,dataset,75);


thresh = max(ssim(abs(beamform_1l), beamform_2{1}) , ...
    ssim(abs(beamform_1r), beamform_2{2}));

queue = {};
queue{end+1} = [1,(1+dataset.firings)/2,dataset.firings];


disp(['Starting image reconstruction from ',acquisition,' for ',phantom,' using ',data,' dataset'])
while ~isempty(queue)
    [beamform_curr, pw_indices, queue] = selection_angles_seq_win(scan,dataset,pht,pw_indices,thresh,...
        queue, beamform_curr);
    queue(1) = [];  % Remove first element
end
    
pw_to_be_put = {pw_indices};
reg_us = das_iq_windowed(scan,dataset,pw_to_be_put,path_scan, path_pht, flag_simu, flag_display);

dynamic_range = 60;
reg_us.show(dynamic_range);
saveas(gcf,'Results/Simulation/Threshold limit/final_result_windowed.jpg');
reg_us.write_file('Results/Simulation/Threshold limit/final_result_windowed.hdf5');
tools.exec_evaluation_contrast_speckle(path_scan,path_pht,'Results/Simulation/Threshold limit/final_result_windowed.hdf5', ...
    flag_simu,flag_display,'Results/Simulation/Threshold limit/final_result_windowed.txt');
save('Results/Simulation/Threshold limit/variables_windowed.mat', 'pw_indices');
disp('Reconstruction Done')
 