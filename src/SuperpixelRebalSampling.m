function [ valid_superpixels, new_spClusters ] = SuperpixelRebalSampling(spClusters, ...
  nSuperpixels, nClusters)
%Sample superpixels in order to balance the classes.

new_spClusters = spClusters;

figure(130);
hg = histogram(spClusters, nClusters + 2);
hg = hg.Values(3:end); %Removes -1 and 0.
low_class = min(hg);

valid_superpixels = 1:nSuperpixels;
for i = 1:nClusters
  class_members_idxs = find(new_spClusters == i);
  class_samples_idxs = randsample(class_members_idxs, low_class, false);

  class_removals = setdiff(class_members_idxs, class_samples_idxs);

  %Indexing scheme resolves the sorting issue (?).
  new_spClusters = new_spClusters(setdiff(1:length(new_spClusters), class_removals));
  valid_superpixels = valid_superpixels(setdiff(1:length(valid_superpixels), class_removals));
end

end

