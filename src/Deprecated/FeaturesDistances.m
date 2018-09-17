function distance = FeaturesDistances(spfv1, spfM2)
  %Function pdist2 tries to compute element-wise first.
  nSP = size(spfM2,1);
  assert(nSP ~= 1); 

  %@anonymous minmax
  minmax = @(x) (x - min(x))/(max(x) - min(x)); 
  
  %Sync with function SuperpixelsFeatures.
  nBins = 10;
  %mean(int) std(int) hist(int) peak(grad) gabor sift saliency
  int_mean = 1;
  int_std = 1;
  int_hist = nBins;
  grad_peak = 1;
  %gabor
  sift = 128;
  sal1 = 1;
  sal2 = 1;
  gabor = (length(spfv1) - (int_mean + int_std + int_hist + grad_peak + sift + sal1 + sal2));

  descsLen = [int_mean int_std int_hist grad_peak gabor sift sal1 sal2];
  descsBounds = [0 cumsum(descsLen)];  

  %% Distances computation
  spfM1 = repmat(spfv1, nSP, 1);

  %Compute distances of scalar features
  scalar_feats = descsBounds((descsLen == 1)) + 1;
  scalar_dists = abs(spfM1(:,scalar_feats) - spfM2(:,scalar_feats));
  %Compute distances of vector features
  gab_difs = spfM1(:,(descsBounds(5)+1):descsBounds(5+1)) - ...
                    spfM2(:,(descsBounds(5)+1):descsBounds(5+1));
	sift_difs = spfM1(:,(descsBounds(6)+1):descsBounds(6+1)) - ...
                    spfM2(:,(descsBounds(6)+1):descsBounds(6+1));
  vec_dists = [sqrt(sum(gab_difs.^2, 2)) sqrt(sum(sift_difs.^2, 2))];
  %Compute histogram distances
  hist_dist = match_distance(spfv1(1,(descsBounds(3)+1):descsBounds(4)), ...
    spfM2(:,(descsBounds(3)+1):descsBounds(4)));
%   hist_dist = minmax(hist_dist);
  
  scalar_weights = [1 1 1 1 1]';
  scalar_weights = scalar_weights/sum(scalar_weights);
  combined_dists = [scalar_dists*scalar_weights vec_dists hist_dist];    

%   figure(748);
%   for i = 1:5
%     subplot(5,1,i); histogram(scalar_dists(:,i), 0:(1/(2*nBins)):(2*1));
%   end
%   figure(749);
%   for i = 1:2
%     subplot(2,1,i); histogram(vec_dists(:,i));
%   end
%   figure(750);
%   for i = 1:4
%     subplot(4,1,i); histogram(combined_dists(:,i));
%     %, 0:(1/(2*nBins)):(2*1));
%   end
%   pause(0.1);
%   drawnow;
  
  total_weights = [1 1/3*ones(1,2) 1]';
  distance = combined_dists*total_weights;
end