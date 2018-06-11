function [lab_out, bsspft] = CopyClosestSuperpixelFromClassColor(samples, source, target)
%Each superpixel receives the color of the closest superpixel from the
%majority class.
kK = 15;

%best super pixel samples for each target
[bsspft, dists] = knnsearch(samples.fv_sp', target.fv_sp', 'K', kK);

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
        avg_sp = sum(sum(mask_c))/length(find(mask_c));
        lab_out(:,:,c) = lab_out(:,:,c) + avg_sp*tgt_mask;
    end
    
end


end

