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

function greedy_sel_iter = greedy_selection_method(phantom_type, data_type,...
                      flag_display, dataset, scan, pht, K, init_set)
                  
    dbstop if error;
    
    N = numel(init_set);
    greedy_sel = init_set;
    greedy_sel_iter = zeros(1,K);
    greedy_sel_iter(1:N) = init_set;
    remaining_images = round(1:dataset.firings);
    remaining_images(init_set) = [];
    score_graph = [];
    
    log_file = fopen('greedy_selection_log.txt', 'w');

    %-- Greedy selection iterations
    for iter = (N+1):K
        % fprintf(log_file,['Selecting image ',num2str(iter), '\n']);
        ind = 1;
        score_max = -Inf;
        for i = 1:numel(remaining_images)
            pw_indices{1} = round(union(greedy_sel, remaining_images(i)));
            % fprintf(log_file,['Candidate set: ',num2str(pw_indices{1}), '\n']);
            clear image score;
            switch data_type    
                case 1
                    image = das_iq_original(scan,dataset,pw_indices);
                case 2
                    image = das_rf(scan,dataset,pw_indices);
                otherwise       %-- Do deal with bad values
                    image = das_iq(scan,dataset,pw_indices);       
            end
            switch phantom_type    
                case 1 	%-- evaluating resolution and distorsion
                    scores = recon_evaluation_resolution_distorsion(scan,pht,image,flag_display);
                    score = scores.resolution;
                case 2 	%-- evaluating contrast and speckle quality
                    scores = tools.recon_evaluation_contrast_speckle(scan,pht,image,flag_display);
                    score = scores.contrast;
                otherwise       %-- Do deal with bad values
                    scores = recon_evaluation_resolution(scan,pht,image,flag_display);
                    score = scores.resolution;
            end
            if (score > score_max)
                score_max = score;
                ind = i;
            end
            % fprintf(log_file,['Current Score: ', num2str(score), '\n']);
        end
        
        score_graph = union(score_graph, score_max);
        greedy_sel = union(greedy_sel, remaining_images(ind));
        greedy_sel_iter(iter) = remaining_images(ind);
        remaining_images(ind) = [];
        fprintf(log_file,['Current selection: ', num2str(greedy_sel_iter(iter)), '\n']);
        % fprintf(log_file,['Current set: ',num2str(greedy_sel), '\n']);
    end
    plot(score_graph);
    fprintf(log_file,'Greedy selection completed');
    fclose(log_file);
                  
end