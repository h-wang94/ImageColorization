function relabels = EdgeAwareRelabeling(clusters, labels, clCosts)
  relabels = zeros(size(labels));
  if (~isempty(labels))
    for ci = 1:length(clusters)
      cluster_labels = labels(clusters{ci});
      cluster_labels = cluster_labels(cluster_labels ~= -1);
      if (~isempty(cluster_labels))
        relabels(clusters{ci}) = mode(cluster_labels);
      else
        relabels(clusters{ci}) = -1;
      end
    end
  else
    if ( sum(sum(clCosts,2) - ones(size(clCosts,1),1) < 1e-5) == size(clCosts,1) )
      %Scores
      for ci = 1:length(clusters)
        scores = clCosts(clusters{ci},:);
        [~, argmax] = max(sum(scores));
        relabels(clusters{ci}) = argmax;
      end
    else
      %Costs
      for ci = 1:length(clusters)
        costs = clCosts(clusters{ci},:);
        [~, argmin] = min(sum(costs));
        relabels(clusters{ci}) = argmin;
      end
    end
  end
    
end

