function [sp_labels, sp_labels_lin, sp_centroids, n_sp] = SuperpixelExtraction(image, n_superpixels)
%Segment input image into superpixels.
%TODO: Change algorithm for Turbopixels
[sp_labels, n_sp] = superpixels(image, n_superpixels);

sp_labels_lin = reshape(sp_labels, size(image, 1)*size(image, 2), 1);

%centroids computation
sp_centroids = zeros(2, n_sp);
for i = 1:n_sp
    [rs, cs] = find(sp_labels == i);
    sp_centroids(:, i) = [mean(rs) mean(cs)]';
end

end

