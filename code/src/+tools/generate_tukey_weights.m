function weight_matrices = generate_tukey_weights(rows, cols, centers, D_rows, D_cols, alpha, distance_type)
    num_centers = length(centers);
    weight_matrices = zeros(rows, cols, num_centers);
    
    % Create coordinate grids
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Compute weight matrices for each center
    for i = 1:num_centers
        cx = centers{i}(1);
        cy = centers{i}(2);
        
        
        % Tukey window function
        weights = zeros(size(X)); % Initialize with ones
        temp = tukeywin(floor(D_rows)-20, alpha) * tukeywin(floor(D_cols), alpha)';
        mask1 = (abs(X - cx) < ((1-alpha) * D_cols / 2)) & (abs(Y - cy) < ((1-alpha) * D_rows / 2));
        mask3 = (1-((abs(X - cx) < ((1-alpha) * D_cols / 2)) & (abs(Y - cy) < ((1-alpha) * D_rows / 2)))) & ...
        ((abs(X - cx) < D_cols / 2) & (abs(Y - cy) < D_rows / 2));
        mask = mask1 | mask3;

        % Compute tapered cosine for start and end regions
        weights(mask1 | mask3) = temp;
        
        % Store in weight matrix
        weight_matrices(:,:,i) = weights;
    end
    
    % Normalize so weights sum to 1 at each pixel
    sum_weights = sum(weight_matrices, 3);
    sum_weights(sum_weights == 0) = 1; % Avoid division by zero
    weight_matrices = weight_matrices ./ sum_weights;
end
