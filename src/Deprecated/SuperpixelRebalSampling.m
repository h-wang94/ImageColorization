function [ valid_superpixels, new_spClusters ] = SuperpixelRebalSampling(spClusters)
%Sample superpixels in order to balance the classes.

new_spClusters = spClusters;
rangeClusters = max(spClusters) - min(spClusters) + 1;
nClusters = max(spClusters);

figure(130);
hg = histogram(spClusters, rangeClusters);
hg = hg.Values((1+end-nClusters):end); %Removes -1 and 0.
low_class = min(hg);

valid_superpixels = 1:length(spClusters);
for i = 1:nClusters
  class_members_idxs = find(new_spClusters == i);
  class_samples_idxs = randsample(class_members_idxs, low_class, false);

  class_removals = setdiff(class_members_idxs, class_samples_idxs);

  %Indexing scheme resolves the sorting issue (?).
  new_spClusters = new_spClusters(setdiff(1:length(new_spClusters), class_removals));
  valid_superpixels = valid_superpixels(setdiff(1:length(valid_superpixels), class_removals));
end

end

