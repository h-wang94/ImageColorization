function [lab_out] = CopyClosestFeatureColor(colors, target, bsft)
%Color transfer technique that assigns for each target pixel the color the
%source pixel with smallest distance in feature space.

lab_out = zeros(size(target.image, 1), size(target.image, 2), 3);
lab_out(:,:,1) = target.luminance*100;

im_sz = size(target.luminance);
for i = 1:(size(target.image,1)*size(target.image,2))
    [r, c] = ind2sub(im_sz,i);
    lab_out(r, c, 2:3) = colors(:, bsft(i));
end

end