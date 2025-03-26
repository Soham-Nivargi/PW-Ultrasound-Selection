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
% path_reconstruted_img = ['../../reconstructed_image/',acquisition,'/',phantom,'/trying_',phantom,'_',acqui,'_img_from_',data,'_K_',num2str(K),'.hdf5'];
path_phantom = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_phantom.hdf5'];


%-- Read the corresponding dataset and the region where to reconstruct the image
dataset = us_dataset();
dataset.read_file(path_dataset);
scan = linear_scan();
scan.read_file(path_scan);


left_top = zeros(1, 75);
left_middle = zeros(1,75);
left_bottom = zeros(1,75);

middle_top = zeros(1, 75);
middle_middle = zeros(1,75);
middle_bottom = zeros(1,75);

right_top = zeros(1, 75);
right_middle = zeros(1,75);
right_bottom = zeros(1,75);

overall = zeros(1,75);

for i = 1:75

    path_reconstruted_img = ['../../reconstructed_image/',acquisition,'/',phantom,'/single_image','/single_image_',phantom,'_',acqui,'_img_from_',data,'_',num2str(i),'.hdf5'];

    %-- Perform evaluation for resolution
    % disp(['Starting evaluation from ',acquisition,' for ',phantom,' using ',data,' dataset'])
    
    %-- Pass online string instances
    flag_simu = num2str(flag_simu);
    flag_display = num2str(flag_display);
    switch phantom_type    
        case 1 	%-- evaluating resolution and distorsion
            tools.exec_evaluation_resolution_distorsion(path_scan,path_phantom,windowed_path,flag_simu,flag_display,path_window_output_log);
        case 2 	%-- evaluating contrast and speckle quality
            score_contrast = tools.value_contrast_ret(path_scan,path_phantom,path_reconstruted_img,flag_simu,flag_display, 1);
        otherwise       %-- Do deal with bad values
            tools.exec_evaluation_resolution(path_scan,path_phantom,path_reconstruted_img,flag_simu,flag_display,path_output_log);
    end

    middle_top(i) = score_contrast(1);
    middle_middle(i) = score_contrast(2);
    middle_bottom(i) = score_contrast(3);

    left_top(i) = score_contrast(4);
    left_middle(i) = score_contrast(5);
    left_bottom(i) = score_contrast(6);

    right_top(i) = score_contrast(7);
    right_middle(i) = score_contrast(8);
    right_bottom(i) = score_contrast(9);

    overall(i) = mean(score_contrast);

    % disp('Evaluation Done')
    % disp(['Result saved in "',path_output_log,'"'])
end
% chosen_angles = [0];
% contrast_list = [];
% 
% image = das_iq(scan,dataset,pw_indices, path_scan, path_pht, flag_simu, flag_display, phantom_type);
% 
% contrast_list.append()
% Initialize variables
angles = [38]; % Start with the initial angle 0°
CNR_values = []; % Store calculated CNR for each angle
max_iterations = 10; % Number of angles to select
explore_range = 10; % Range of degrees to consider around each chosen angle

% Initial CNR calculation for the 0° angle
CNR_values(end+1) = overall(angles(end));

for iter = 2:max_iterations
    % Fit a linear regression model based on the angles and their CNR values
    model = fitlm(angles, CNR_values); % Simple linear regression model
    
    % Predict CNR for angles within a small range in both directions
    pos_candidate = angles(end) + explore_range; % Next positive angle
    neg_candidate = angles(end) - explore_range; % Next negative angle
    pos_mirror = -pos_candidate; % Mirror of the positive candidate
    neg_mirror = -neg_candidate; % Mirror of the negative candidate
    
    % Predict the CNR values for all candidates
    CNR_pred_pos = predict(model, pos_candidate);
    CNR_pred_neg = predict(model, neg_candidate);
    CNR_pred_pos_mirror = predict(model, pos_mirror);
    CNR_pred_neg_mirror = predict(model, neg_mirror);
    
    % Introduce a small randomization factor to avoid getting stuck
    random_offset = randi([-5, 5]);
    pos_candidate = pos_candidate + random_offset;
    neg_candidate = neg_candidate - random_offset;
    pos_mirror = pos_mirror + random_offset;
    neg_mirror = neg_mirror - random_offset;
    
    % Store all candidates and their predicted CNRs in a table
    candidates = [pos_candidate, neg_candidate, pos_mirror, neg_mirror];
    CNR_preds = [CNR_pred_pos, CNR_pred_neg, CNR_pred_pos_mirror, CNR_pred_neg_mirror];
    
    % Select the candidate with the highest predicted CNR
    [~, max_idx] = max(CNR_preds);
    next_angle = candidates(max_idx);
    
    % Calculate the CNR for the chosen angle and add it to the lists
    angles(end+1) = next_angle;
    CNR_values(end+1) = overall(next_angle);
    
    % Display the selected angle and corresponding CNR value
    fprintf('Iteration %d: Selected angle = %.2f°, CNR = %.4f\n', iter, next_angle, CNR_values(end));
end

% Display final chosen angles and corresponding CNR values
disp('Final chosen angles and corresponding CNR values:');
disp(table(angles', CNR_values', 'VariableNames', {'Angle', 'CNR'}));


















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
