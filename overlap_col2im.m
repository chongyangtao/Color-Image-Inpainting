function I = overlap_col2im(X, Mask, bb, overlap, im_size)
% Recover overlapping (bb x bb) patches into an image reconstruction
%
% X: (bb^2 x N) matrix containing the vectorised signals
% Mask: (bb^2 x N) matrix containing the vectorised mask. Each column is the mask of the corresponding column of X
% bb: patch size
% overlap: number of pixels between consecutive patches
% im_size: tuple of original image dimensions
%


% Get image dimensions
N = im_size(1);
M = im_size(2);

% Initialise the image
I = zeros(N, M);

% Initialise the matrix that stores the cumulative signal used to calculate
% each pixel
S = zeros(N, M);

% Calculate how many patches fit horisontally and vertically
row = (N - bb)/overlap;
col = (M - bb)/overlap;

% Iterate through patches and keep a counter to know what column of X to
% access next
counter = 0;
for i = 0:row
    for j = 0:col
        counter = counter + 1;
        % Get x and y range of current patch
        x_patch = [1:bb]+i*overlap;
        y_patch = [1:bb]+j*overlap;
        % Get existing patch
        existing_patch = I(x_patch, y_patch);
        % Get current patch and calculate its signal strength
        current_patch = reshape(X(:,counter), bb, bb);
        signal = sum(sum(Mask(:,counter)~=0));
        % Update cumulative signal strength
        S(x_patch, y_patch) = S(x_patch, y_patch) + signal;
        % Update patch
        I(x_patch, y_patch) = existing_patch + current_patch*signal;
    end
end
% If a pixel was not covered by a single patch with signal, set its signal
% value to 1, in order to prevent division by zero
S(S == 0) = 1;

% Divide each pixel by its corresponding cumulative signal coverage
I = I./S;
