function [spLabels, spLabelsLin, spCentroids, nSP] = SuperpixelExtraction(image, n_superpixels, type)
%Segment input image into superpixels.

[phi, ~, disp_img, spLabels] = turbosuperpixels(image, n_superpixels);

if (strcmp(type, 'turbo'))
  [~, ~, ~, spLabels] = turbosuperpixels(image, n_superpixels);
  nSP = max(max(spLabels));
else
  [spLabels, nSP] = superpixels(image, n_superpixels);
end

spLabelsLin = reshape(spLabels, size(image, 1)*size(image, 2), 1);

%Centroids computation
spCentroids = zeros(2, nSP);
for i = 1:nSP
    [rs, cs] = find(spLabels == i);
    spCentroids(:, i) = [mean(rs) mean(cs)]';
end

end

