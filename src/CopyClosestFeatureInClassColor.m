function [lab_out] = CopyClosestFeatureInClassColor(colors, target, ...
  neighbor_idxs, neighbor_classes, labels)
%Performs a classification on each pixel and assigns either the color of 
%the the class centroid or the color of the closest pixel in feature space
%from the class it was assigned during classification.

lab_out = zeros(size(target.image));
lab_out(:,:,1) = target.luminance*100;

%% Transfer

for i = 1:(size(target.image,1)*size(target.image,2))
  % Instances from chosen class
  [~, majority_instances] = find(neighbor_classes(i,:) == labels(i));

  % Color assignment
  [r, c] = ind2sub(size(target.luminance),i);
  lab_out(r, c, 2:3) = mean(colors(:, neighbor_idxs(i,majority_instances)), 2);
end

end

