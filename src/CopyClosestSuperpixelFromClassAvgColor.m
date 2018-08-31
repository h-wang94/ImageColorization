function lab_out = CopyClosestSuperpixelFromClassAvgColor(source, target, ...
  neighbor_idxs, neighbor_classes, labels)
%Each superpixel receives the color of the closest superpixel from the
%majority class.

%Output image
lab_out = zeros([size(target.image) 3]);
lab_out(:,:,1) = target.luminance*100;

for c = 2:3
  for i = 1:target.nSuperpixels
    if (labels(i) == -1)
      tgt_mask = (target.sp==i);
      lab_out(:,:,1) = lab_out(:,:,1).*~tgt_mask;
      %if label is marked as doubt, assign black (for debug purposes)  
      continue
    end
    
    % Instances from chosen class
    [~, majority_instances] = find(neighbor_classes(i,:) == labels(i));

    %Matching superpixels ROI masks
    tgt_mask = (target.sp==i);
    src_mask = (source.sp==neighbor_idxs(i,majority_instances(1)));

%     if (isempty(majority_instances));
%       majority_instances = find(neighbor_classes == labels(i));
%       src_mask = (source.sp==neighbor_idxs(majority_instances(1)));
%     else
%       src_mask = (source.sp==neighbor_idxs(i, majority_instances(1)));
%     end
    
    %Prototype color transfer (Superpixel average)
    mask_c = source.lab(:,:,c).*src_mask;
    avg_sp = sum(sum(mask_c))/length(find(src_mask));
    lab_out(:,:,c) = lab_out(:,:,c) + avg_sp*tgt_mask;
  end
end

end

