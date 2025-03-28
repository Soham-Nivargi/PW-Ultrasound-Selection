%-- Script to be used as an example to evaluate the quality of the reconstructed images

%-- After choosing the specific configuration through acquisition_type, 
%-- phantom_type, data_type and flag_display parameters, this script allows computing the 
%-- chosen evaluation metrics.

%-- By default, the computed metrics are saved in a text file in the folder named "evaluation"

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
flag_display = 0;           %-- 0 = do not display || 1 = display intermediate results

K = 5;

%-- Parse parameter choices
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
    case 1	%-- evaluating resolution and distorsion
        phantom = 'resolution_distorsion';
    case 2	%-- evaluating contrast and speckle quality
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
path_scan = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_scan.hdf5'];
path_phantom = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_phantom.hdf5'];
path_reconstruted_img = ['../../reconstructed_image/',acquisition,'/',phantom,'/trying_',phantom,'_',acqui,'_img_from_',data,'_K_',num2str(K),'.hdf5'];
path_output_log = ['../../evaluation/',phantom,'/','trying_',phantom,'_',acqui,'_evaluation_from_',data,'1.txt'];


%-- Perform evaluation for resolution
disp(['Starting evaluation from ',acquisition,' for ',phantom,' using ',data,' dataset'])

%-- Pass online string instances
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


