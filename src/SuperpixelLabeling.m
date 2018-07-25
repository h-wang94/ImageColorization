function [SPClusterLabels] = SuperpixelLabeling(pClustersIdxs, spIdxsLin, LABEL_MAJORITY, ...
  PLOT, spIdxs, imgSize)
%SuperpixelLabeling: Define class labels for superpixels
% Based on label majority of contained pixels.
% Requires at least LABEL_MAJORITY% to label, otherwise labeled as doubt (-1).

SPClusterLabels = zeros(1, max(spIdxsLin));
for i = 1:length(SPClusterLabels)
  sp_labels = pClustersIdxs(spIdxsLin == i);

  % If majority is less than kMajTol %, mark as doubt.
  if(sum(mode(sp_labels) == sp_labels) < LABEL_MAJORITY*length(sp_labels))
    SPClusterLabels(i) = -1;
  else
    SPClusterLabels(i) = mode(sp_labels);
  end
end

if (PLOT)
  src_col_labels = CreateLabeledImage(SPClusterLabels, spIdxs, imgSize);
  figure; imshow(src_col_labels, []); title('Automatic Labeling of Superpixels');
  colormap jet; drawnow;
end

end

