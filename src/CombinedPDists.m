function PD  = CombinedPDists(spfM1, spfM2, featsWeigths)
  %@anonymous minmax
%   minmax = @(x) (x - min(x))/(max(x) - min(x)); 
  
  %Sync with function SuperpixelsFeatures.
  nBins = 10;
  %mean(int) std(int) hist(int) peak(grad) gabor sift saliency
  int_mean = 1;
  int_std = 1;
  int_hist = nBins;
  grad_mean = 1;
  grad_peak = 1;
  %gabor
  sift = 128;
  sal1 = 1;
  sal2 = 1;
  gabor = (size(spfM1,1) - (int_mean + int_std + int_hist + grad_peak + grad_mean + sift + sal1 + sal2));

  descsLen = [int_mean int_std int_hist grad_mean grad_peak gabor sift sal1 sal2];
  descsBounds = [0 cumsum(descsLen)];  
  
  %% Distances computations
  src_nSP = size(spfM1,2);
  
  PDs = {};
  kIdxs = {};
  for pdi = 1:length(descsLen)
    idxs = (descsBounds(pdi)+1):descsBounds(pdi+1);
%     [kIdxs{pdi}, PDs{pdi}] = knnsearch(spfM1(idxs,:)', spfM2(idxs,:)', 'K', src_nSP);
    PDs{pdi} = pdist2(spfM1(idxs,:)', spfM2(idxs,:)');
    PDs{pdi} = PDs{pdi} / max(PDs{pdi}(:));
  end
  idxs = (descsBounds(3)+1):descsBounds(3+1);
  PDs{3} = pdist2(spfM1(idxs,:)', spfM2(idxs,:)', @match_distance);
  PDs{3} = PDs{3} / max(PDs{3}(:));
  
%   figure;
%   for i = 1:length(PDs)
%     subplot(3,3, i); histogram(PDs{i}(:));
%     title(num2str(i));
%   end
  
  PD = zeros(size(PDs{1}));
  featsWeigths = featsWeigths/sum(featsWeigths);
  for i = 1:length(PDs)
    PD = PD + featsWeigths(i)*PDs{i};
  end
  PD = PD';
end