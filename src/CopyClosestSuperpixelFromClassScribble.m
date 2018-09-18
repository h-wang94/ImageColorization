function [lab_out, scribbles_mask] = CopyClosestSuperpixelFromClassScribble(source, target, ...
  neighbor_idxs, neighbor_classes, labels, Kct)
%TODO

%Output image
lab_out = zeros([size(target.image) 3]);
lab_out(:,:,1) = target.luminance*100;

scribbles_mask = zeros(size(target.image));
for c = 2:3
  for i = 1:target.nSuperpixels
    if (labels(i) == -1)
      continue;
    end
    % Instances from chosen class
    [~, majority_instances] = find(neighbor_classes(i,:) == labels(i));

    %Matching superpixels ROI masks
    src_mask = zeros(size(source.sp));
    for ni = 1:min([Kct length(majority_instances)])
      src_mask = src_mask | (source.sp==neighbor_idxs(i,majority_instances(ni)));
    end
    cntrd = round(target.sp_centroids(:,i));
    
    %Prototype color transfer (Superpixel average)
    mask_c = source.lab(:,:,c).*src_mask;
    avg_sp = sum(sum(mask_c))/length(find(src_mask));
    lab_out(cntrd(1), cntrd(2), c) = avg_sp;
    scribbles_mask(cntrd(1), cntrd(2)) = 1;
  end
end

end

