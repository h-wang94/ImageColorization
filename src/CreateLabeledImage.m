function [im_labeled] = CreateLabeledImage(sp_labels, image_sps, image_size)
%
im_labeled = zeros(image_size);
for i = 1:length(sp_labels)
  mask = (image_sps == i);
  im_labeled = im_labeled + sp_labels(i)*mask;
end

im_labeled(1,1) = -1; %For comparison with classification

end

