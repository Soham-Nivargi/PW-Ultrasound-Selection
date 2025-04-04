%-- Script to be used to greedily select a subset from the
%-- provided dataset and reconstruct

%-- After choosing the specific configuration through acquisition_type, 
%-- phantom_type and data_type parameters, this script allows reconstructing
%-- images for evaluation for different choices of steered plane waves involved
%-- in the compounding scheme (specified by the pw_indices parameter)

%-- The implemented method (to be used as example) corresponds to the standard Delay
%-- And Sum (DAS) technique with apodization in reception

%-- Authors: Kaushani Majumder (kaushanim@iitb.ac.in)

%-- $Date: 2023/04/05 $  


clear all;
close all;
clc;
addpath(genpath('../src'));
addpath(genpath('tools'));

%-- Input Parameters
K = 15;                      %-- Number of images to be selected for reconstruction
init_set = 38;
                             %-- Initial set of images to start with

%-- Parameters
acquisition_type = 1;       %-- 1 = simulation || 2 = experiments
phantom_type = 2;           %-- 1 = resolution & distorsion || 2 = contrast & speckle quality
data_type = 1;              %-- 1 = IQ || 2 = RF
flag_display = 0;           %-- 0 = do not display || 1 = display intermediate results


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
path_phantom = ['../../../../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_phantom.hdf5'];
% path_reconstruted_img = ['oct30/',acquisition,'/greedy','/greedy_selection',phantom,'_',acqui,'_img_from_',data,'_K_',num2str(K),'.hdf5'];


%-- Read the corresponding dataset and the region where to reconstruct the image
dataset = us_dataset();
dataset.read_file(path_dataset);
scan = linear_scan();
scan.read_file(path_scan);
pht = us_phantom();
pht.read_file(path_phantom);

%-- Total number of plane-wave images in the dataset
N = round(dataset.firings);                %-- Total number of images available


%-- Greedy selection of images
disp('In greedy selection algorithm');
[beamform_curr, greedy_sel_iter] = greedy_selection_method(scan,pht,dataset, K, init_set);
                  
disp(['Selected set of images (sorted by iteration): ', num2str(pw)]);

dynamic_range = 60;
reg_us = us_image('DAS-IQ beamforming');
reg_us.author = 'Soham Nivargi';
reg_us.affiliation = 'IIT Bombay';
reg_us.algorithm = 'Delay-and-Sum (IQ version)';
reg_us.scan = scan;
reg_us.number_plane_waves = K;
reg_us.data = abs(beamform_curr);
reg_us.transmit_f_number = 0;
reg_us.receive_f_number = 1.75;
reg_us.transmit_apodization_window = 'none';
reg_us.receive_apodization_window = 'Tukey 50%';
reg_us.show();
saveas(gcf,'Results/simulation/final_result_unwindowed.jpg');
reg_us.write_file('Results/simulation/final_result_unwindowed.hdf5');
tools.exec_evaluation_contrast_speckle(path_scan,path_phantom,'Results/simulation/final_result_unwindowed.hdf5', ...
    flag_simu,flag_display,'Results/simulation/final_result_unwindowed.txt');
% save('Results/simulation/variables_unwindowed.mat', 'pw_indices');
disp('Reconstruction Done')

%-- Save results
% image.write_file(path_reconstruted_img);