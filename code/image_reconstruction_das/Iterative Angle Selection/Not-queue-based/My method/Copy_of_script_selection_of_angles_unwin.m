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

% beamform_l = (abs(das_iq_one_image(scan,dataset,1))+abs(das_iq_one_image(scan,dataset,2)))/2;
% beamform_r = (abs(das_iq_one_image(scan,dataset,74))+abs(das_iq_one_image(scan,dataset,75)))/2;
% 
% beamform_le = abs(das_iq_one_image(scan,dataset,1));
% beamform_re = abs(das_iq_one_image(scan,dataset,75));
% 
% thresh = max(ssim(beamform_l, beamform_le) , ...
%     ssim(beamform_r, beamform_re));

% queue = {};
% queue{end+1} = [1,(1+dataset.firings)/2,dataset.firings];


disp(['Starting image reconstruction from ',acquisition,' for ',phantom,' using ',data,' dataset'])
% while length(pw_indices{1})<15
[beamform_curr, pw_indices, thresh] = Copy_of_selection_angles_seq_unwin(scan,dataset,pw_indices,14,...
        1, (dataset.firings+1)/2, dataset.firings, beamform_curr);
% end
    
dynamic_range = 60;
reg_us = us_image('DAS-IQ beamforming');
reg_us.author = 'Soham Nivargi';
reg_us.affiliation = 'IIT Bombay';
reg_us.algorithm = 'Delay-and-Sum (IQ version)';
reg_us.scan = scan;
reg_us.number_plane_waves = length(pw_indices);
reg_us.data = abs(beamform_curr);
reg_us.transmit_f_number = 0;
reg_us.receive_f_number = 1.75;
reg_us.transmit_apodization_window = 'none';
reg_us.receive_apodization_window = 'Tukey 50%';
reg_us.show(dynamic_range);
saveas(gcf,'Results/Simulation/Angle Limit/final_result_unwindowed.jpg');
reg_us.write_file('Results/Simulation/Angle Limit/final_result_unwindowed.hdf5');
tools.exec_evaluation_contrast_speckle(path_scan,path_pht,'Results/Simulation/Angle Limit/final_result_unwindowed.hdf5', ...
    flag_simu,flag_display,'Results/Simulation/Angle Limit/final_result_unwindowed.txt');
save('Results/Simulation/Angle Limit/variables_unwindowed.mat', 'pw_indices');
disp('Reconstruction Done')
 