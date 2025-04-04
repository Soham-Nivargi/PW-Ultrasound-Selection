function [beamform_curr, pw_indice, queue] = selection_angles_seq_win(scan,pht,dataset,pw_indice,thresh,queue,beamform_curr)
    
    mid = floor(mean(queue{1}));

    if ismember(mid, pw_indice)
        return
    end
    
    pw_indices{2} = [queue{1}(1),mid];
    pw_indices{3} = [queue{1}(2), mid];
    pw_indices{1} = [queue{1},mid];
    %-- receive apodization: 
    %-- dynamically expanding receive aperture with Tukey 25% apodization

    [beamform_res, beamform_comp] = das_iq_for_selection_window_og(scan,pht,dataset,pw_indices);
    
    temp_us = us_image('DAS-IQ beamforming');
    temp_us.author = 'Soham Nivargi';
    temp_us.affiliation = 'IIT Bombay';
    temp_us.algorithm = 'Delay-and-Sum (IQ version)';
    temp_us.scan = scan;
    temp_us.transmit_f_number = 0;
    temp_us.receive_f_number = 1.75;
    temp_us.transmit_apodization_window = 'none';
    temp_us.receive_apodization_window = 'Tukey 25%';

    beamform_curr = beamform_res{1};


    left_sim = ssim(beamform_comp{1}, beamform_comp{2});
    right_sim = ssim(beamform_comp{1}, beamform_comp{3});
    pw_indice = [pw_indice, mid];

    if left_sim < thresh
        if(right_sim<thresh)
            if left_sim>right_sim
                queue{end+1} = pw_indices{2};
                queue{end+1} = pw_indices{3};
            else
                queue{end+1} = pw_indices{3};
                queue{end+1} = pw_indices{2};
            end
        else
            queue{end+1} = pw_indices{2};
        end
    else
        if(right_sim<thresh)
            queue{end+1} = pw_indices{3};
        end
    end
end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     pw_indices{2} = [queue{1}(1),mid];
%     pw_indices{3} = [queue{1}(2),mid];
%     pw_indices{4} = [queue{1}(1),queue{1}(2),mid];
%     pw_indices{1} = [pw_indices{1}, mid];
% 
%     beamform_left = das_iq_one_image(scan,dataset,queue{1}(1));
%     beamform_right = das_iq_one_image(scan,dataset,queue{1}(2));
%     beamform_mid = das_iq_one_image(scan,dataset,mid);
% 
%     beamform_next{1} = (abs(beamform_mid)+abs(beamform_left))/2;
% 
%     beamform_next{2} = (abs(beamform_mid)+abs(beamform_right))/2;
% 
%     beamform_comp = (abs(beamform_left)+abs(beamform_mid)+abs(beamform_right))/3;
% 
%     left_sim = ssim(beamform_comp, beamform_next{1});
%     right_sim = ssim(beamform_comp, beamform_next{2});
% 
%     beamform_curr = beamform_curr + beamform_mid;
% 
%     if left_sim < thresh
%         if(right_sim<thresh)
%             if left_sim>right_sim
%                 queue{end+1} = [queue{1}(1),mid];
%                 queue{end+1} = [mid,queue{1}(2)];
%             else
%                 queue{end+1} = [mid,queue{1}(2)];
%                 queue{end+1} = [queue{1}(1),mid];
%             end
%         else
%             queue{end+1} = [queue{1}(1),mid];
%         end
%     else
%         if(right_sim<thresh)
%             queue{end+1} = [mid,queue{1}(2)];
%         end
%     end
% end