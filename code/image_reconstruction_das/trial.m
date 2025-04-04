% Define image size
rows = 1000;
cols = 1000;
transition_size = 50; % Number of rows for smooth transition

% Create base intensities
image = ones(rows, cols) * 0.5; % Set entire image to lower intensity
image(1:500, :) = 0.3; % Upper half at intensity 1

% Generate smooth transition mask using cosine or sigmoid
transition_region = (500-transition_size/2):(500+transition_size/2);
blend = 0.5 * (1 + cos(linspace(-pi, pi, length(transition_region))))'; % Cosine blending

% Apply blending to transition rows
for col = 1:cols
    image(transition_region, col) = blend * 0.3 + (1 - blend) * 0.5;
end

% Display the image
imshow(image, []);
