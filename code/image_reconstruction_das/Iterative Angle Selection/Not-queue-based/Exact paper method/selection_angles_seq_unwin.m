function [beamform_curr, pw_indices] = selection_angles_seq_unwin(scan,dataset,pw_indices,thresh,left,right,beamform_curr)
    % image = 0;
    mid = floor((left+right)/2);
    
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

    if left_sim < thresh
        if(right_sim<thresh)
            if left_sim>right_sim
                % queue{end+1} = [queue{1}(1),left_mid,queue{1}(2)];
                % queue{end+1} = [queue{1}(2),right_mid,queue{1}(3)];
                [beamform_curr, pw_indices] = selection_angles_seq_unwin(scan,dataset,pw_indices,thresh,...
                    left, mid, beamform_curr);
                [beamform_curr, pw_indices] = selection_angles_seq_unwin(scan,dataset,pw_indices,thresh,...
                    mid, right, beamform_curr);
            else
                % queue{end+1} = [queue{1}(2),right_mid,queue{1}(3)];
                % queue{end+1} = [queue{1}(1),left_mid,queue{1}(2)];
                [beamform_curr, pw_indices] = selection_angles_seq_unwin(scan,dataset,pw_indices,thresh,...
                    mid, right, beamform_curr);
                [beamform_curr, pw_indices] = selection_angles_seq_unwin(scan,dataset,pw_indices,thresh,...
                    left, mid, beamform_curr);
            end
        else
            [beamform_curr, pw_indices] = selection_angles_seq_unwin(scan,dataset,pw_indices,thresh,...
                    left, mid, beamform_curr);
        end
    else
        if(right_sim<thresh)
            [beamform_curr, pw_indices] = selection_angles_seq_unwin(scan,dataset,pw_indices,thresh,...
                    mid, right, beamform_curr);
        end
    end
end