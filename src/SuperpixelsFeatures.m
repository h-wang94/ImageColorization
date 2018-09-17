function [tgt_spfv, src_spfv] = SuperpixelsFeatures(source, samples, target)
%Computes the features of superpixels based on its pixel elements.

nBins = 10;  %TODO: Parametererize!
% nFeats = length(target.fvl);
% ftsBounds = cumsum(target.fvl);
% ftsBounds = [0 ftsBounds];
% minmax = @(x) (x - min(x))/(max(x) - min(x));

%% Feature Aggregation 1:
% %- Scalars: remain the same.
% %- Vectors: apply statistic.
% agg = { @(x) x, ...
%         @(x) x, ...
%         @(x) minmax(std(x)), ...
%         @(x) minmax(median(x))};
% 
% t_fv = zeros(nFeats, size(target.fv,2));
% s_fv = zeros(nFeats, size(samples.fv,2));
% for fi = 1:nFeats
%   t_fv(fi,:) = agg{fi}(target.fv(ftsBounds(fi)+1:ftsBounds(fi+1),:));
%   s_fv(fi,:) = agg{fi}(samples.fv(ftsBounds(fi)+1:ftsBounds(fi+1),:));
% end

tgt_spfv = fillSPFV(target.fv, target.lin_sp, target.nSuperpixels, nBins, target.fvl);
src_spfv = fillSPFV(samples.fv, source.lin_sp, source.nSuperpixels, nBins, target.fvl);

end

function sp_fv = fillSPFV(fv, linSPIdxs, nSP, nBins, fvLens)
  minmax = @(x) (x - min(x))/(max(x) - min(x));
  
  %
  nFeats = length(fvLens);
  ftsBounds = [0 cumsum(fvLens)];
  
  %Mappings between pixel -> superpixel features. 
  % Each row represents the operations over each feature.
  binsEdges = 0:(1/nBins):1;
  dirEdges = -180:(360/(2*nBins)):180;
  descripts = { {@(x) mean(x), @(x) std(x), @(x) histcounts(x, binsEdges) / length(x)}, ...
                {@(x) mean(x)}, ...
                {@(x) max(histcounts(x, dirEdges))/length(x) }, ...  % [(1/nBins),1] do not use histcounts (rotation)
                {@(x) mean(x, 2)}, ...
                {@(x) mean(x, 2)}};
  descsLen = [1 1 nBins 1 1 fvLens(end-1) fvLens(end)];
  descsBounds = [0 cumsum(descsLen)];
  
    
  %Compute the Superpixel Features
  sp_fv = zeros(sum(descsLen), nSP);
  for spi = 1:nSP
    mask = (linSPIdxs == spi);
    
    desc_idx = 1;
    for fi = 1:nFeats
      %Pixel features belongin to spi
      p_feats_spi = fv((ftsBounds(fi)+1):ftsBounds(fi+1),mask);
      
%       if(fi == 3)
%         histogram(p_feats_spi,dirEdges);
%         pause(0.1);
%       end
            
      for si = 1:length(descripts{fi})
        sp_fv((descsBounds(desc_idx)+1):descsBounds(desc_idx+1), spi) = descripts{fi}{si}(p_feats_spi);
        desc_idx = desc_idx + 1;
      end
    end
  end
  
  %MinMax should not be used!
  % Normalization is performed on each image separately but then used for
  % comparison against a different image.
  
%   for di = 1:length(descsLen)
%     if (descsLen(di) == 1)
%       sp_fv(di,:) = minmax(sp_fv(di,:));
%     end
%   end

%   for di = 1:sum(descsLen)
%     sp_fv(di,:) = minmax(sp_fv(di,:));
%   end
   
end
