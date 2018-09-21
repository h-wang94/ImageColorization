function [lab_out, scribbles_mask] = CopyClosestSuperpixelAvgScribble(source, target, neighbor_idxs)
%TODO

%Output image
lab_out = zeros([size(target.image) 3]);
lab_out(:,:,1) = target.luminance*100;

scribbles_mask = zeros(size(target.image));
for i = 1:length(neighbor_idxs)
  % Masks
  src_mask = (source.sp==neighbor_idxs(i));
  cntrd = round(target.sp_centroids(:,i));
  
  %Prototype color transfer (Superpixel average)
  for c = 2:3
    mask_c = source.lab(:,:,c).*src_mask;
    avg_sp = sum(sum(mask_c))/length(find(src_mask));
    lab_out(cntrd(1), cntrd(2), c) = avg_sp;
    scribbles_mask(cntrd(1), cntrd(2)) = 1;
  end
end

end

