function newLabels = FeatureClustersRelabeling(target, labels)
%
newLabels = zeros(size(labels));
tgt_feat_labels = CreateLabeledImage(target.fv_sp_labels, target.sp, size(target.image));

for fli = 1:max(target.fv_sp_labels)
  set_i = (tgt_feat_labels == fli);
  CC = bwconncomp(set_i);

  for cci = 1:CC.NumObjects
    %Linear index of pixels of each connected region
    cci_p_idxs = CC.PixelIdxList(cci);

    %Convert pixel indexes to superpixel indexes
    cci_sp_idxs = unique(target.lin_sp(cci_p_idxs{1}));

    %TODO: caso especial quando for pequeno (absorvido pela redondeza)

    %Find set of color labels assigned to the pixels of current region
    cci_labels = labels(cci_sp_idxs);
    cci_labels = cci_labels(cci_labels ~= -1);
    if (isempty(cci_labels))
      newLabels(cci_sp_idxs) = -1;
    else
      newLabels(cci_sp_idxs) = mode(cci_labels);
    end
  end
end


end

