function weight_matrices = generate_weight_matrix(rows, cols, centers, sigma, distance_type)
    % Get the number of centers
    num_centers = length(centers);
    
    % Initialize the weight matrix (3D array)
    weight_matrices = zeros(rows, cols, num_centers);
    
    % Create coordinate grids
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Compute weight matrices for each center
    for i = 1:num_centers
        cx = centers{i}(1);
        cy = centers{i}(2);
        
        % Compute distance from the center
        if strcmp(distance_type, 'manhattan')
            D = abs(X - cx) + abs(Y - cy);  % Manhattan Distance
        else
            D = sqrt((X - cx).^2 + (Y - cy).^2);  % Euclidean Distance
        end
        
        % Apply Gaussian or exponential kernel
        weight_matrices(:,:,i) = exp(-D.^2 / (2 * sigma^2));  % Gaussian kernel
        % weight_matrices(:,:,i) = exp(-D / sigma);  % Exponential alternative
    end
    
    % Normalize across the third dimension so that weights sum to 1 at each pixel
    sum_weights = sum(weight_matrices, 3);
    % sum_weights(sum_weights == 0) = 1; % Prevent division by zero
    weight_matrices = weight_matrices ./ sum_weights;
end
