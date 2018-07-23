function ShowDistancesFSvsCS(samples, source, kSP, nClusters)
%
if (kSP)
%Superpixel
  %Compute average color of superpixels
  sp_ab = zeros(2, length(source.validSuperpixels));
  for spi = 1:length(sp_ab)
    vsi = source.validSuperpixels(spi);
    vsi_idxs = find(source.lin_sp == vsi);

    sp_ab(:,spi) = mean(samples.ab(:,vsi_idxs),2);
  end
  
  %Pairwise distance in Chrominance and Feature spaces.
  cdist = pdist(sp_ab');
  fdist = pdist(source.fv_sp');
  CD = squareform(cdist);
  FD = squareform(fdist);

  %For each spixel, show distances to all other spixels separated by class
  figure; 
  for i = 1:size(CD,1)
    for ci = 1:nClusters
      spci_idxs = source.sp_clusters == ci;
      scatter(FD(spci_idxs,i), CD(spci_idxs,i), '.');
      hold on;
    end
    hold off;
    xlabel('Feature distance');
    ylabel('Color distance');
    drawnow
%     pause;
  end

else
%Pixel
  fdist = pdist(samples.fv(2:end,:)');
  cdist = pdist(samples.ab');
  
  figure; scatter(fdist, cdist, '.');
  xlabel('Feature distance');
  ylabel('Color distance');
end    

end

