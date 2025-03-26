%-- Function to perform greedy selection on the given dataset to select 'K'
%-- images for reconstruction. Each greedy step includes image
%-- reconstruction for the candidate subset and evaluation of the performance criteria
%-- on the reconstructed image.

%-- Input parameter: K = number of images to select

%-- Input parameters include the specific configuration through acquisition_type, 
%-- phantom_type and data_type parameters

%-- The implemented method (to be used as example) corresponds to the standard Delay
%-- And Sum (DAS) technique with apodization in reception

%-- Authors: Kaushani Majumder (kaushanim@iitb.ac.in)

%-- $Date: 2023/04/05 $  

function reg_us = greedy_selection_window(phantom_type, acquisition_type,...
                      flag_simu, flag_display, dataset, scan, pht, K, init_set)
                  
    dbstop if error;
    
    N = numel(init_set);
    greedy_sel = init_set;
    greedy_sel_iter = zeros(1,K);
    greedy_sel_iter(1:N) = init_set;
    remaining_images = round(1:dataset.firings);
    remaining_images(init_set) = [];
    score_graph = [];

    pw_indices{1} = greedy_sel;

    iter_beamform = reshape(zeros([numel(scan.z_axis) numel(scan.x_axis)  1]), ...
        [numel(scan.z_axis) numel(scan.x_axis) length(pw_indices)]);
    contrast_reg = zeros([1 9]);

    
    switch acquisition_type
        case 1
            [iter_beamform, max_ind, contrast_reg] = das_iq_for_greedy(scan,pht,dataset,pw_indices, ...
                flag_simu,flag_display,iter_beamform, contrast_reg, 1);
        case 2
            [iter_beamform, max_ind, contrast_reg] = das_iq_for_greedy_exp(scan,pht,dataset,pw_indices, ...
                flag_simu,flag_display,iter_beamform, contrast_reg, 1);
        otherwise       %-- Do deal with bad values
            image = das_iq(scan,dataset,pw_indices);       
    end

    log_file = fopen('greedy_selection_log.txt', 'w');

    %-- Greedy selection iterations
    for iter = (N+1):K
        % fprintf(log_file,['Selecting image ',num2str(iter), '\n']);
        ind = 1;
        score_max = -Inf;
        pw_indices = {};
        for i = 1:numel(remaining_images)
            pw_indices{i} = remaining_images(i);
        end
        % disp(pw_indices{38});
        % pw_indices = {};
        % fprintf(log_file,['Candidate set: ',num2str(pw_indices{1}), '\n']);
        clear image score;
        switch acquisition_type    
            case 1
                [iter_beamform, max_ind, contrast_reg] = das_iq_for_greedy(scan,pht,dataset,pw_indices, ...
                    flag_simu,flag_display,iter_beamform, contrast_reg, 0);
            case 2
                [iter_beamform, max_ind, contrast_reg] = das_iq_for_greedy_exp(scan,pht,dataset,pw_indices, ...
                    flag_simu,flag_display,iter_beamform, contrast_reg, 0);
            otherwise       %-- Do deal with bad values
                image = das_iq(scan,dataset,pw_indices);       
        end

        score_max = mean(contrast_reg, 2);
        % ind = 
        % switch phantom_type    
        %     case 1 	%-- evaluating resolution and distorsion
        %         scores = recon_evaluation_resolution_distorsion(scan,pht,image,flag_display);
        %         score = scores.resolution;
        %     case 2 	%-- evaluating contrast and speckle quality
        %         scores = tools.recon_evaluation_contrast_speckle(scan,pht,image,flag_display);
        %         score = scores.contrast;
        %     otherwise       %-- Do deal with bad values
        %         scores = recon_evaluation_resolution_distortion(scan,pht,image,flag_display);
        %         score = scores.resolution;
        % end
        % if (score > score_max)
        %     score_max = score;
        %     ind = i;
        % end
        % fprintf(log_file,['Current Score: ', num2str(score), '\n']);
        % % end
        % fprintf(log_file,['Current image selected: ', num2str(remaining_images(ind)), '\n']);
        score_graph = union(score_graph, score_max);
        greedy_sel = union(greedy_sel, remaining_images(max_ind));
        greedy_sel_iter(iter) = remaining_images(max_ind);
        remaining_images(max_ind) = [];
        
        fprintf(log_file,['Current selection: ',num2str(greedy_sel_iter(iter)), '\n']);
    end
    dynamic_range = 60;
    reg_us = us_image('DAS-IQ beamforming');
    reg_us.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
    reg_us.affiliation = 'Norwegian University of Science and Technology (NTNU)';
    reg_us.algorithm = 'Delay-and-Sum (IQ version)';
    reg_us.scan = scan;
    reg_us.number_plane_waves = K;
    reg_us.data = abs(iter_beamform);
    reg_us.transmit_f_number = 0;
    reg_us.receive_f_number = 1.75;
    reg_us.transmit_apodization_window = 'none';
    reg_us.receive_apodization_window = 'Tukey 50%';
    reg_us.show(dynamic_range);
    % % saveas(gcf, filename1);
    plot(score_graph);
    fprintf(log_file,'Greedy selection completed');
    fclose(log_file);
                  
end