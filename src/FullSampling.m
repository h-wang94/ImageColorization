function [idxs, ab] = FullSampling(image_lab)

%Fully sampled
sz = size(image_lab);
[r, c] = ind2sub(sz, 1:sz(1)*sz(2));
idxs = [r; c];

ab = reshape(image_lab(:,:,2:3), size(image_lab,1)*size(image_lab,2), 2);
ab = ab';

end

