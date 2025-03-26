%script to select angles iteratively
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

%-- Create path to load corresponding files
path_dataset = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_dataset_',data,'.hdf5'];
path_scan = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_scan.hdf5'];
path_pht = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_phantom.hdf5'];


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
if max(abs(beamform_curr(:)))>=1
    beamform_curr1 = beamform_curr/max(abs(beamform_curr(:)));
else
    beamform_curr1 = beamform_curr;
end

% beamform_l = das_iq_one_image(scan,dataset,1)+das_iq_one_image(scan,dataset,2);
% beamform_r = das_iq_one_image(scan,dataset,74)+das_iq_one_image(scan,dataset,75);
% 
% if max(abs(beamform_l(:)))>=1
%     beamform_l1 = beamform_l/max(abs(beamform_l(:)));
% else
%     beamform_l1 = beamform_l;
% end
% 
% if max(abs(beamform_r(:)))>=1
%     beamform_r1 = beamform_r/max(abs(beamform_r(:)));
% else
%     beamform_r1 = beamform_r;
% end
% 
% beamform_le = das_iq_one_image(scan,dataset,1);
% beamform_re = das_iq_one_image(scan,dataset,75);
% 
% if max(abs(beamform_le(:)))>=1
%     beamform_le1 = beamform_le/max(abs(beamform_le(:)));
% else
%     beamform_le1 = beamform_le;
% end
% 
% if max(abs(beamform_re(:)))>=1
%     beamform_re1 = beamform_re/max(abs(beamform_re(:)));
% else
%     beamform_re1 = beamform_re;
% end
% 
% thresh = max(ssim(abs(beamform_l1), abs(beamform_le1)) , ...
%     ssim(abs(beamform_r1), abs(beamform_re1)));

queue = {};
queue{end+1} = [1,(1+dataset.firings)/2,dataset.firings];


disp(['Starting image reconstruction from ',acquisition,' for ',phantom,' using ',data,' dataset'])
while ~isempty(queue)
    [beamform_curr, pw_indices, queue] = selection_angles_seq_win(scan,dataset,pht,pw_indices,0.99,...
        queue, beamform_curr);
    queue(1) = [];  % Remove first element
end
    
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
saveas(gcf,'Results/Sampling Strategies/Paper_modifications/windowed/simulation/queue_final_result.jpg');
reg_us.write_file('Results/Sampling Strategies/Paper_modifications/windowed/simulation/queue_final_result.hdf5');
tools.exec_evaluation_contrast_speckle(path_scan,path_pht,'Results/Sampling Strategies/Paper_modifications/windowed/simulation/queue_final_result.hdf5', ...
    flag_simu,flag_display,'Results/Sampling Strategies/Paper_modifications/windowed/simulation/queue_final_result.txt');

disp('Reconstruction Done')
 