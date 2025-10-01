function struct = FISH_ComputeMaskGradient(struct)
% FISH_ComputeMaskGradient - Computes the gradient of the mask and fills in NaN values.
%
% Syntax: sst = FISH_ComputeMaskGradient(sst)
%
% Input:
%   sst - A structure containing the following fields:
%         mask - A 2D array representing the mask where NaN values are to be processed.
%
% Output:
%   sst - The same input structure with an additional field:
%         mask_grad - A 2D array representing the gradient of the mask with NaN values filled.
%
% Example:
%   sst.mask = [1 NaN 3; 4 5 NaN; NaN 7 8];
%   sst = FISH_ComputeMaskGradient(sst);
%   disp(sst.mask_grad);

% Step 1: Create a temporary copy of the mask and replace NaNs with 0s
mask_tmp = struct.mask;
mask_tmp(isnan(mask_tmp)) = 0;

% Step 2: Compute the gradient of the mask
[~, ~, mask_grad, ~] = computeGradient(mask_tmp);

% Step 3: Replace 0s in the gradient with NaNs
mask_grad(mask_grad == 0) = nan;

% Step 4: Fill in NaNs in the gradient using nearest neighbor interpolation
mask_grad = inpaint_nans(mask_grad, 4);

% % Step 5: Create a temporary mask to handle original NaNs
% mask_tmp(isnan(sst.mask)) = 1;
% mask_tmp(~isnan(sst.mask)) = nan;
% 
% % Step 6: Multiply the filled gradient with the temporary mask to retain NaN positions
% mask_grad = mask_grad .* mask_tmp;

% Step 7: Add the computed gradient to the output structure
struct.mask_grad = mask_grad;

end
