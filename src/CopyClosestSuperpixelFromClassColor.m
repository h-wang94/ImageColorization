function [lab_out, bsspft] = CopyClosestSuperpixelFromClassColor(source, target, K)
%Each superpixel receives the color of the closest superpixel from the
%majority class.

%best super pixel samples for each target
[bsspft, dists] = knnsearch(source.fv_sp', target.fv_sp', 'K', K);

%Output image
lab_out = zeros([size(target.image) 3]);
lab_out(:,:,1) = target.luminance*100;

for i = 1:max(target.lin_sp)
    %Find classes of each source superpixel
    class_hipts = source.sp_clusters(bsspft(i,:));
    
    %if there is no majority, returns the first which is the closest.
    class_majority = mode(class_hipts);
    [~, majority_instances] = find(class_hipts == class_majority);

    %Matching superpixels ROI masks
    tgt_mask = (target.sp==i);
    src_mask = (source.sp==bsspft(i,majority_instances(1)));
    
    %Prototype color transfer (Superpixel average)
    for c = 2:3
        mask_c = source.lab(:,:,c).*src_mask;
        avg_sp = sum(sum(mask_c))/length(find(src_mask));
        lab_out(:,:,c) = lab_out(:,:,c) + avg_sp*tgt_mask;
%         figure(c*100); imshow(lab_out(:,:,c),[]);
    end
    
end


end

