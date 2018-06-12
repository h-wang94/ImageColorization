function [tgt_spfv, source_spfv] = SuperpixelFeatures(source, samples, target)
%Builds the feature vector of superpixels by averaging 
    tgt_nSP = max(target.lin_sp);
    src_nSP = max(source.lin_sp);
    
    tgt_spfv = zeros(size(target.fv,1), tgt_nSP);
    source_spfv = zeros(size(target.fv,1), src_nSP);

    for i = 1:tgt_nSP
        tgt_spfv(:,i) = mean(target.fv(:, target.lin_sp == i), 2);
    end
    for i = 1:src_nSP
        source_spfv(:,i) = mean(samples.fv(:, source.lin_sp == i), 2);
    end

end

