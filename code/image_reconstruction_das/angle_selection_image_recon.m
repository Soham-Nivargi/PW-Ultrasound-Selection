%-- Script to be used as an example to manipulate the provided dataset

%-- After choosing the specific configuration through acquisition_type, 
%-- phantom_type and data_type parameters, this script allows reconstructing
%-- images for evaluation for different choices of steered plane waves involved
%-- in the compounding scheme (specified by the pw_indices parameter)

%-- The implemented method (to be used as example) corresponds to the standard Delay
%-- And Sum (DAS) technique with apodization in reception

%-- Authors: Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)
%--          Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)

%-- $Date: 2016/03/01 $  


clear all;
close all;
clc;
addpath(genpath('../src'));


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

K = 5;
%-- Create path to load corresponding files
path_dataset = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_dataset_',data,'.hdf5'];
path_scan = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_scan.hdf5'];
path_reconstruted_img = ['../../reconstructed_image/',acquisition,'/',phantom,'/trying_',phantom,'_',acqui,'_img_from_',data,'_K_',num2str(K),'.hdf5'];
path_pht = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_phantom.hdf5'];


%-- Read the corresponding dataset and the region where to reconstruct the image
dataset = us_dataset();
dataset.read_file(path_dataset);
scan = linear_scan();
scan.read_file(path_scan);

chosen_angles = [0];
contrast_list = [];

image = das_iq(scan,dataset,pw_indices, path_scan, path_pht, flag_simu, flag_display, phantom_type);

contrast_list.append()


















%-- Indices of plane waves to be used for each reconstruction
pw_indices{1} = [29 47 15 60];
% pw_indices{2} = round(linspace(1,dataset.firings,3));
% pw_indices{3} = round(linspace(1,dataset.firings,11));
% pw_indices{4} = round(1:dataset.firings);               %-- dataset.firings corresponding to the total number of emitted steered plane waves
flag_simu = num2str(flag_simu);
flag_display = num2str(flag_display);

%-- Reconstruct Bmode images for each pw_indices
disp(['Starting image reconstruction from ',acquisition,' for ',phantom,' using ',data,' dataset'])
switch data_type    
    case 1
        image = das_iq(scan,dataset,pw_indices, path_scan, path_pht, flag_simu, flag_display, phantom_type);
    case 2
        image = das_rf(scan,dataset,pw_indices);
    otherwise       %-- Do deal with bad values
        image = das_iq(scan,dataset,pw_indices);       
end
disp('Reconstruction Done')
disp(['Result saved in "',path_reconstruted_img,'"'])


%-- Show the corresponding beamformed images
dynamic_range = 60;
image.show(dynamic_range);


%-- Save results
image.write_file(path_reconstruted_img);

path_output_log = ['../../evaluation/',phantom,'/','windowed_nahi_',phantom,'_',acqui,'_evaluation_from_',data,'.txt'];
path_phantom = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_phantom.hdf5'];

flag_simu = num2str(flag_simu);
flag_display = num2str(flag_display);
switch phantom_type    
    case 1 	%-- evaluating resolution and distorsion
        tools.exec_evaluation_resolution_distorsion(path_scan,path_phantom,path_reconstruted_img,flag_simu,flag_display,path_output_log);
    case 2 	%-- evaluating contrast and speckle quality
        tools.exec_evaluation_contrast_speckle(path_scan,path_phantom,path_reconstruted_img,flag_simu,flag_display,path_output_log);
    otherwise       %-- Do deal with bad values
        tools.exec_evaluation_resolution(path_scan,path_phantom,path_reconstruted_img,flag_simu,flag_display,path_output_log);
end
disp('Evaluation Done')
disp(['Result saved in "',path_output_log,'"'])
