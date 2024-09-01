function score = recon_evaluation_contrast_speckle(scan,pht,image,flag_display)

    %-- Function to evaluate the CNR score and the speckle noise

    %-- Convert input argument received as string
%     flag_display = uint8(str2num(flag_display));   	%-- convert string back to int
    
    %-- Perform testing for contrast
    testing_contrast = us_contrast();
    testing_contrast.pht = pht;
    testing_contrast.scan = scan;
    testing_contrast.image = image;
    testing_contrast.flagDisplay = flag_display;
    testing_contrast.evaluate();
    
    %-- Perform testing for speckle quality
    testing_speckle = us_speckle_quality();
    testing_speckle.pht = pht;
    testing_speckle.scan = scan;
    testing_speckle.image = image;
    testing_speckle.flagDisplay = flag_display;
    testing_speckle.evaluate();
    
    %-- Final output scores
    score.contrast = mean(testing_contrast.score,2);
    
    if ~isempty(find(testing_speckle.score==0,1))
        score.speckle = 0;
    else
        score.speckle = 1;
    end

end