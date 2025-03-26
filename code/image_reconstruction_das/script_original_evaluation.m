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

K = 9;

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

K=9;
%-- Create path to load corresponding files
path_scan = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_scan.hdf5'];
path_phantom = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_phantom.hdf5'];
% path_reconstruted_img = ['../../reconstructed_image/',acquisition,'/',phantom,'/single_image','/single_image_',phantom,'_',acqui,'_img_from_',data,'_75','.hdf5'];
path_reconstruted_img1 = ['Results/uniform_simulation/new_window/unwindowed_15', '.hdf5'];
path_output_log1 = ['Results/uniform_simulation/new_window/unwindowed_15', '.txt'];

% path_reconstruted_img2 = ['Results/uniform_experiment/yes_window/windowed', '.hdf5'];
% path_output_log2 = ['Results/uniform_experiment/yes_window/windowed', '.txt'];

% path_reconstruted_img1 = ['Results/uniform_experiment/no_window/uniform','_K_',num2str(9),'.hdf5'];
% path_output_log1 = ['Results/uniform_experiment/no_window/uniform','_K_',num2str(9),'.txt'];

% path_window_output_log = ['tukey_regionwise_2/trying_window_wala_image', '.txt'];
% path_unwindow_output_log = ['tukey_regionwise_2/trying_non_window_wala_image', '.txt'];
% windowed_path = ['tukey_regionwise_2/trying_half_window_wala_image', '.hdf5'];
% unwindowed_path = ['tukey_regionwise_2/trying_non_window_wala_image', '.hdf5'];
% left_top = zeros(1, 75);
% left_middle = zeros(1,75);
% left_bottom = zeros(1,75);
% 
% middle_top = zeros(1, 75);
% middle_middle = zeros(1,75);
% middle_bottom = zeros(1,75);
% 
% right_top = zeros(1, 75);
% right_middle = zeros(1,75);
% right_bottom = zeros(1,75);
% 
% overall = zeros(1,75);
scan = linear_scan();
scan.read_file(path_scan);
pht = us_phantom();
pht.read_file(path_phantom);
image = us_image();
image.read_file(path_reconstruted_img1);
disp(image.number_plane_waves);
% for i = 1:75
%     if i==38
%         path_reconstruted_img = ['../../reconstructed_image/',acquisition,'/',phantom,'/single_image','/single_image_',phantom,'_',acqui,'_img_from_',data,'_',num2str(i),'.hdf5'];
%         cat = 1;
%     else
%         path_reconstruted_img = ['../../reconstructed_image/',acquisition,'/',phantom,'/double_image','/double_image_',phantom,'_',acqui,'_img_from_',data,'_',num2str(i),'.hdf5'];
%         cat = 2;
%     end
    %-- Perform evaluation for resolution
    disp(['Starting evaluation from ',acquisition,' for ',phantom,' using ',data,' dataset'])
    % score = tools.contrast_score(scan, pht, image, flag_simu, flag_display, 5);
    % disp(score);
    %-- Pass online string instances
    flag_simu = num2str(flag_simu);
    flag_display = num2str(flag_display);
    switch phantom_type    
        case 1 	%-- evaluating resolution and distorsion
            tools.exec_evaluation_resolution_distorsion(path_scan,path_phantom,windowed_path,flag_simu,flag_display,path_window_output_log);
        case 2 	%-- evaluating contrast and speckle quality
            % score_contrast = tools.value_contrast_ret(path_scan,path_phantom,path_reconstruted_img,flag_simu,flag_display, cat);
            tools.exec_evaluation_contrast_speckle(path_scan,path_phantom,path_reconstruted_img1,flag_simu,flag_display,path_output_log1);
            % tools.exec_evaluation_contrast_speckle(path_scan,path_phantom,path_reconstruted_img2,flag_simu,flag_display,path_output_log2);
            % disp(score_contrast);
        otherwise       %-- Do deal with bad values
            tools.exec_evaluation_resolution(path_scan,path_phantom,path_reconstruted_img,flag_simu,flag_display,path_output_log);
    end
    
    % middle_top(i) = score_contrast(1);
    % middle_middle(i) = score_contrast(2);
    % middle_bottom(i) = score_contrast(3);
    % 
    % left_top(i) = score_contrast(4);
    % left_middle(i) = score_contrast(5);
    % left_bottom(i) = score_contrast(6);
    % 
    % right_top(i) = score_contrast(7);
    % right_middle(i) = score_contrast(8);
    % right_bottom(i) = score_contrast(9);
    % 
    % overall(i) = mean(score_contrast);

    disp('Evaluation Done');
    % disp(['Result saved in "',path_output_log1,'"'])
% end
% % left_top(38) = sum(left_top)/74;
% % left_middle(38) = sum(left_middle)/74;
% % left_bottom(38) = sum(left_bottom)/74;
% % middle_top(38) = sum(middle_top)/74;
% % middle_middle(38) = sum(middle_middle)/74;
% % middle_bottom(38) = sum(middle_bottom)/74;
% % right_top(38) = sum(right_top)/74;
% % right_middle(38) = sum(right_middle)/74;
% % right_bottom(38) = sum(right_bottom)/74;
% [~, idx_lt] = max(left_top);
% [~, idx_lm] = max(left_middle);
% [~, idx_lb] = max(left_bottom);
% 
% [~, idx_mt] = max(middle_top);
% [~, idx_mm] = max(middle_middle);
% [~, idx_mb] = max(middle_bottom);
% 
% [~, idx_rt] = max(right_top);
% [~, idx_rm] = max(right_middle);
% [~, idx_rb] = max(right_bottom);
% 
% [~, idx_ov] = max(overall);
% 
% 
% 
% 
% 
% figure();
% plot(left_top-left_top(38));
% xlabel('PW Index');
% ylabel('CNR');
% title('left top');
% saveas(gcf, 'double/left_top.jpg');
% 
% figure();
% plot(left_middle-left_middle(38));
% xlabel('PW Index');
% ylabel('CNR');
% title('left middle');
% saveas(gcf, 'double/left_middle.jpg');
% 
% figure();
% plot(left_bottom-left_bottom(38));
% xlabel('PW Index');
% ylabel('CNR');
% title('left bottom');
% saveas(gcf, 'double/left_bottom.jpg');
% 
% figure();
% plot(middle_top-middle_top(38));
% xlabel('PW Index');
% ylabel('CNR');
% title('middle top');
% saveas(gcf, 'double/middle_top.jpg');
% 
% figure();
% plot(middle_middle-middle_middle(38));
% xlabel('PW Index');
% ylabel('CNR');
% title('middle middle');
% saveas(gcf, 'double/middle_middle.jpg');
% 
% figure();
% plot(middle_bottom-middle_bottom(38));
% xlabel('PW Index');
% ylabel('CNR');
% title('middle bottom');
% saveas(gcf, 'double/middle_bottom.jpg');
% 
% figure();
% plot(right_top-right_top(38));
% xlabel('PW Index');
% ylabel('CNR');
% title('right top');
% saveas(gcf, 'double/right_top.jpg');
% 
% figure();
% plot(right_middle-right_middle(38));
% xlabel('PW Index');
% ylabel('CNR');
% title('right middle');
% saveas(gcf, 'double/right_middle.jpg');
% 
% figure();
% plot(right_bottom-right_bottom(38));
% xlabel('PW Index');
% ylabel('CNR');
% title('right bottom');
% saveas(gcf, 'double/right_bottom.jpg');
% 
% figure();
% plot(overall-overall(38));
% xlabel('PW Index');
% ylabel('CNR');
% title('Overall');
% saveas(gcf, 'double/overall.jpg');