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

function [iter_beamform, greedy_sel_iter] = greedy_selection_method(scan,pht,dataset, K, init_set)
                  
    dbstop if error;
    
    N = numel(init_set);
    greedy_sel = init_set;
    greedy_sel_iter = zeros(1,K);
    greedy_sel_iter(1:N) = init_set;
    remaining_images = round(1:dataset.firings);
    remaining_images(init_set) = [];

    iter_beamform = das_iq_paper_unwindowed(scan, dataset, greedy_sel);
    
    log_file = fopen('greedy_selection_log.txt', 'w');

    %-- Greedy selection iterations
    for iter = (N+1):K
        % fprintf(log_file,['Selecting image ',num2str(iter), '\n']);
        pw_indices = {};
        for i = 1:numel(remaining_images)
            pw_indices{i} = remaining_images(i);
        end
        [iter_beamform,max_ind] = das_iq_for_greedy_unwin(scan,pht,dataset,pw_indices,iter_beamform);
        % fprintf(log_file,['Current Score: ', num2str(score), '\n']);     
        % score_graph = union(score_graph, score_max);
        greedy_sel = union(greedy_sel, remaining_images(max_ind));
        greedy_sel_iter(iter) = remaining_images(max_ind);
        remaining_images(max_ind) = [];
        fprintf(log_file,['Current selection: ', num2str(greedy_sel_iter(iter)), '\n']);
        % fprintf(log_file,['Current set: ',num2str(greedy_sel), '\n']);
    end
    % plot(score_graph);
    fprintf(log_file,'Greedy selection completed');
    fclose(log_file);
                  
end