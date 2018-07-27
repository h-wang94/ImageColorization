function iiDists = ClusterDistNormMedian(ldata, mcCost, compType)
%TODO: Change function to implement the median idea.

maxLen = length(find(ldata(:,end) ~= -1));
maxLen = (maxLen*maxLen - maxLen)/2;

% Initialization
intra = zeros(1, maxLen);
inter = zeros(1, maxLen);
i_intra = 1;
i_inter = 1;

% Computation
if (strcmp(compType,'loop'))  
  for i = 1:size(ldata,1)
    for j = (i+1):size(ldata,1)
      pwDist = norm(ldata(i,1:end-1) - ldata(j,1:end-1));

      if (ldata(i,end) ~= -1 && ldata(j,end) ~= -1)
        if ldata(i,end) == ldata(j,end)
          intra(i_intra) = pwDist;
          i_intra = i_intra + 1;
        else
          inter(i_inter) = pwDist/mcCost(ldata(i,end),ldata(j,end));
          i_inter = i_inter + 1;
        end
      end
    end
  end
else
  pwDist = squareform(pdist(ldata(:,1:end-1)));
  intra = [];
  inter = [];
  %For each class
  for i = 1:max(ldata(:,end))
    maskSame = (ldata(:,end) == i);
    maskOthers = (and(ldata(:,end) ~= i, ldata(:,end) > 0));
    
    intraDists = pwDist(maskSame, maskSame);
    %Special treatment for diagonal removal.
    intraDists_lin = sort(intraDists(:));
    intraDists_lin = intraDists_lin((length(diag(intraDists))+1):end);
    intraDists_lin = intraDists_lin(1:2:end);
    intra = [intra intraDists_lin];
    
    interDists = pwDist(maskSame, maskOthers);
    inter = [inter median(interDists(:))];

  end
end
intra = median(intra(find(intra)));
inter = median(inter(find(inter)));

iiDists = [intra;inter];
end