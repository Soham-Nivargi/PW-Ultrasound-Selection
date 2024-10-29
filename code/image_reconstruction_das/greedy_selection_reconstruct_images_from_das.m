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
K = 9;                      %-- Number of images to be selected for reconstruction
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
path_dataset = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_dataset_',data,'.hdf5'];
path_scan = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_scan.hdf5'];
path_phantom = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_phantom.hdf5'];
path_reconstruted_img = ['../../reconstructed_image/',acquisition,'/',phantom,'/greedy_selection_sorted_left_top_',phantom,'_',acqui,'_img_from_',data,'_K_',num2str(K),'.hdf5'];


%-- Read the corresponding dataset and the region where to reconstruct the image
dataset = us_dataset();
dataset.read_file(path_dataset);
scan = linear_scan();
scan.read_file(path_scan);
pht = us_phantom();
pht.read_file(path_phantom);

%-- Total number of plane-wave images in the dataset
N = round(dataset.firings);                %-- Total number of images available


%-- Indices of plane waves to be used for each reconstruction
% pw_indices{1} = randperm(N,K);
% pw_indices{2} = round(linspace(1,dataset.firings,3));
% pw_indices{3} = round(linspace(1,dataset.firings,11));
% pw_indices{4} = round(1:dataset.firings);               %-- dataset.firings corresponding to the total number of emitted steered plane waves
% pw_indices = cell(1,K);

%-- Greedy selection of images
disp('In greedy selection algorithm');
pw = greedy_selection_method(phantom_type, data_type,...
                      flag_display, dataset, scan, pht, K, init_set);

% pw = [35,43,28,51,38,21,29,18,46,54,34,15,24,42,57,31,47,39,60,50,61,25,12,41,32];
                  
disp(['Selected set of images (sorted by iteration): ', num2str(pw)]);

% for i = 1:K
%     pw_indices{i} = pw(1:i);
% end

pw_indices{1} = pw;


%-- Reconstruct Bmode images for each pw_indices
disp(['Starting image reconstruction from ',acquisition,' for ',phantom,' using ',data,' dataset'])
switch data_type    
    case 1
        image = das_iq_original(scan,dataset,pw_indices);
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