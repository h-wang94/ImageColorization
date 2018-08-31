function [lab_out] = CopyClosestSuperpixelAvgColor(source, target, neighbor_idxs)
%Perform superpixel matching followed by transfer by copy

%Output image
lab_out = zeros([size(target.image) 3]);
lab_out(:,:,1) = target.luminance*100;

for i = 1:length(neighbor_idxs)
% Masks
  src_mask = (source.sp==neighbor_idxs(i));
  tgt_mask = (target.sp==i);

  %Prototype color transfer (Superpixel average)
  for c = 2:3
    mask_c = source.lab(:,:,c).*src_mask;
    avg_sp = sum(sum(mask_c))/length(find(src_mask));
    lab_out(:,:,c) = lab_out(:,:,c) + avg_sp*tgt_mask;
  end
end

end

