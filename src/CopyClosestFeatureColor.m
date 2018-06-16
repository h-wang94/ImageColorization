function [lab_out, best_dists] = CopyClosestFeatureColor(samples, target)
%Color transfer technique that assigns for each target pixel the color the
%source pixel with smallest distance in feature space.
colors = samples.ab;

lab_out = zeros(size(target.image, 1), size(target.image, 2), 3);
lab_out(:,:,1) = target.luminance*100;

[r, c] = size(target.luminance);

%best samples for each target
[bsft, best_dists] = knnsearch(samples.fv', target.fv');

for i = 1:(size(target.image,1)*size(target.image,2))
    %Index conversion
    [r, c] = ind2sub(size(target.luminance),i);
    lab_out(r, c, 2:3) = colors(:, bsft(i));
end
   

end