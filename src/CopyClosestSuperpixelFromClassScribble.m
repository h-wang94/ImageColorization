function [lab_out, bsspft] = CopyClosestSuperpixelFromClassScribble(source, target, K)
%
%best super pixel samples for each target
[bsspft, dists] = knnsearch(source.fv_sp', target.fv_sp', 'K', K);

%Output image
lab_out = zeros([size(target.image) 3]);
lab_out(:,:,1) = target.luminance*100;

for i = 1:target.nSuperpixels
    %Find classes of each source superpixel
    class_hipts = source.sp_clusters(bsspft(i,:));
    
    %if there is no majority, returns the first which is the closest.
    class_majority = mode(class_hipts);
    [~, majority_instances] = find(class_hipts == class_majority);

    %Matching superpixels ROI masks
    src_mask = (source.sp==bsspft(i,majority_instances(1)));
    
    %Prototype color transfer (Superpixel average)
    cntrd = round(target.sp_centroids(:,i));
    for c = 2:3
        mask_c = source.lab(:,:,c).*src_mask;
        avg_sp = sum(sum(mask_c))/length(find(src_mask));
        lab_out(cntrd(1), cntrd(2), c) = avg_sp;
    end
end


end

