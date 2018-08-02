function [SPClusterLabels, SPChrom] = SuperpixelLabeling(imgLab, pClustersIdxs, SPIdxsLin, ...
  LABEL_MAJORITY, PLOT, spIdxs, imgSize)
%SuperpixelLabeling: Define class labels for superpixels
% Based on label majority of contained pixels.
% Requires at least LABEL_MAJORITY% to label, otherwise labeled as doubt (-1).

SPClusterLabels = zeros(1, max(SPIdxsLin));
SPChrom = zeros(2, max(SPIdxsLin));
for i = 1:length(SPClusterLabels)
  idxLin_i = SPIdxsLin == i;
  sp_labels = pClustersIdxs(idxLin_i);

  % Compute superpixel color median
  mask = (spIdxs == i);
  A = mask.*imgLab(:,:,2);
  B = mask.*imgLab(:,:,3);
  SPChrom(1,i) = median(A(idxLin_i));
  SPChrom(2,i) = median(B(idxLin_i));

  
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

