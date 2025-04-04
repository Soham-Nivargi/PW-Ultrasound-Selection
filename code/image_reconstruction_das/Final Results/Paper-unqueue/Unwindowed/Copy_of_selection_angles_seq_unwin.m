function [beamform_curr, pw_indices] = Copy_of_selection_angles_seq_unwin(scan,dataset,pw_indices,left, right,beamform_curr, max_allowed_pw_indices)

    mid = floor((left+right)/2);
    
    if ismember(mid, pw_indices) || length(pw_indices) == max_allowed_pw_indices
        return
    end

    pw_indice{2} = [left,mid];
    pw_indice{3} = [right,mid];
    pw_indice{4} = [left,right,mid];
    pw_indices = [pw_indices, mid];
    
    beamform_left = das_iq_one_image(scan,dataset,left);
    beamform_right = das_iq_one_image(scan,dataset,right);
    beamform_mid = das_iq_one_image(scan,dataset,mid);

    beamform_next{1} = (abs(beamform_mid)+abs(beamform_left))/2;

    beamform_next{2} = (abs(beamform_mid)+abs(beamform_right))/2;
    
    beamform_comp = (abs(beamform_left)+abs(beamform_mid)+abs(beamform_right))/3;

    left_sim = ssim(beamform_comp, beamform_next{1});
    right_sim = ssim(beamform_comp, beamform_next{2});
    
    beamform_curr = beamform_curr + beamform_mid;

            if left_sim>right_sim
                % queue{end+1} = [queue{1}(1),left_mid,queue{1}(2)];
                % queue{end+1} = [queue{1}(2),right_mid,queue{1}(3)];
                [beamform_curr, pw_indices] = Copy_of_selection_angles_seq_unwin(scan,dataset,pw_indices,...
                    left, mid, beamform_curr, max_allowed_pw_indices);
                [beamform_curr, pw_indices] = Copy_of_selection_angles_seq_unwin(scan,dataset,pw_indices,...
                    mid, right, beamform_curr, max_allowed_pw_indices);
            else
                % queue{end+1} = [queue{1}(2),right_mid,queue{1}(3)];
                % queue{end+1} = [queue{1}(1),left_mid,queue{1}(2)];
                [beamform_curr, pw_indices] = Copy_of_selection_angles_seq_unwin(scan,dataset,pw_indices,...
                    mid, right, beamform_curr, max_allowed_pw_indices);
                [beamform_curr, pw_indices] = Copy_of_selection_angles_seq_unwin(scan,dataset,pw_indices,...
                    left, mid, beamform_curr, max_allowed_pw_indices);
            end
        
end