function iiDists = ClusterDistNormMedian(ldata, mcCost, compType)
%TODO: Change function to implement the median idea.

% Computation
if (strcmp(compType,'loop'))
  disp('This might take a while');
  
  % Initialization
  maxLen = length(find(ldata(:,end) ~= -1));
  maxLen = (maxLen*maxLen - maxLen)/2;
  intra = zeros(1, maxLen);
  inter = zeros(1, maxLen);
  i_intra = 1;
  i_inter = 1;
  
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
  
  intra = median(intra(1:i_intra));
  inter = median(inter(1:i_inter));

  iiDists = [intra; inter];
else
  pwDist = squareform(pdist(ldata(:,1:end-1)));
  %For each class
  for i = 1:max(ldata(:,end))
    maskSame = (ldata(:,end) == i);
    maskOthers = (and(ldata(:,end) ~= i, ldata(:,end) > 0));
    
    intraDists = pwDist(maskSame, maskSame);
    %Special treatment for diagonal removal.
    intraDists_lin = sort(intraDists(:));
    intraDists_lin = intraDists_lin((length(diag(intraDists))+1):end);
    intraDists_lin = intraDists_lin(1:2:end);
    intras{i} = intraDists_lin;
    
    interDists = pwDist(maskSame, maskOthers);
    inters{i} = median(interDists(:));
  end
  intra = zeros(1,max(ldata(:,end)));
  inter = intra;
  for ci = 1:max(ldata(:,end))
    intra(ci) = median(intras{ci});
    inter(ci) = median(inters{ci});
  end
  iiDists = [sum(intra); sum(inter)];
end

end