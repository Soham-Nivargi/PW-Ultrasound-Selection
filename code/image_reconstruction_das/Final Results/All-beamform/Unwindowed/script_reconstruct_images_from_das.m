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
    otherwise       
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

%-- Indices of plane waves to be used for each reconstruction
% pw_indices{1} = [38 34 42 22 54 30 46 26 50 18 58 14 62 10 66]; %uniform
flag_simu = num2str(flag_simu);
flag_display = num2str(flag_display);
%-- Reconstruct Bmode images for each pw_indices
disp(['Starting image reconstruction from ',acquisition,' for ',phantom,' using ',data,' dataset'])
switch data_type    
    case 1
        % disp(arr{i});
        % image1 = das_iq_windowed(scan,dataset,pw_indices, path_scan, path_pht, flag_simu, flag_display);
        image2 = das_iq_original(scan,dataset);
        % image3 = das_iq_trial_windowing(scan, pht, dataset, pw_indices);
    case 2
        image = das_rf(scan,dataset,pw_indices);
    otherwise       %-- Do deal with bad values
        image = das_iq(scan,dataset,pw_indices);       
end


image2.show();
saveas(gcf,'Results/simulation/final_result_unwindowed.jpg');
image2.write_file('Results/simulation/final_result_unwindowed.hdf5');
tools.exec_evaluation_contrast_speckle(path_scan,path_pht,'Results/simulation/final_result_unwindowed.hdf5', ...
    flag_simu,flag_display,'Results/simulation/final_result_unwindowed.txt');
% save('Results/simulation/variables_unwindowed.mat', 'pw_indices');
disp('Reconstruction Done')
 