function [groups, clusters] = EdgeAwareClustering(target)
%

remainingSPs = 1:target.nSuperpixels;
groups = {};

%Compute nnList (image space)
nnList = GetListofNeighbouringSuperPixels(target.sp, 3);

%Compute image edges (and..)
edges = edge(target.luminance, 'canny');
% edges = imdilate(edges, ones(3));

i = 1;
while(~isempty(remainingSPs))
  groups{i} = [remainingSPs(1)];
  remainingSPs = remainingSPs(2:end);
  
  %Expand the group{i} until it reaches its boundaries.
  [groups{i}, remainingSPs] = ExpandGroup(groups{i}, nnList, remainingSPs, target.sp_centroids, edges);
  
  i = i + 1;
end

%Build connection map (label style)
clusters = zeros(size(target.sp));
for gi = 1:length(groups)
  for i = 1:length(groups{gi})
    clusters = clusters + gi*(target.sp==groups{gi}(i));
  end
end
figure(73);
imshow(clusters,[]);
colormap 'jet'

end

function [frontier, remSPs] = ExpandGroup(frontier, nnList, remSPs, centroidsSPs, edges)
  i = 1;
  while(i <= length(frontier))
    %Frontier gets updated by final insertion until cursor reaches the end.
    
    conNeigh = intersect(nnList{frontier(i),3}, remSPs);
    conNeigh = CheckPathEdge(frontier(i), conNeigh, centroidsSPs, edges);
    
    frontier = [frontier setdiff(conNeigh,frontier)];
    remSPs = setdiff(remSPs, conNeigh);
    
    i = i + 1;
  end
end

function [notBlocked] = CheckPathEdge(reference, candidates, centroidsSPs, edges)
%Estimate
  center = @(x,y) (max(x,y)-min(x,y))/2 + min(x,y);
  bounds = @(x,y) floor(min(x,y)):ceil(max(x,y));

  notBlocked = [];
  for i = 1:length(candidates)
    p1 = centroidsSPs(:,reference);
    p2 = centroidsSPs(:,candidates(i));
    
    center_pix = [center(p1(1),p2(1));
                  center(p1(2),p2(2))];
    
    center_pix = round(center_pix);
%     imshow(edges); hold on; scatter(center_pix(2), center_pix(1), 'y.');
%     scatter(p1(2), p1(1), 'g.')
    
    rect = edges(bounds(p1(1),p2(1)),bounds(p1(2),p2(2)));
    if (~sum(rect(:)))
      notBlocked = [notBlocked candidates(i)];
      
%       scatter(p2(2), p2(1), 'b.')
%       drawnow; pause(0.5);
%     else
%       scatter(p2(2), p2(1), 'r.')
%       drawnow; pause(0.5);
    end
  end
end


function Result = GetListofNeighbouringSuperPixels(labels, diskSize)

TotalSuperpixels = max(labels(:));
Result = cell(TotalSuperpixels, 3);

se = strel('disk', diskSize);
for I = 1:TotalSuperpixels
    sp = (labels == I);
    spd = imdilate(sp, se);
    spNN = unique(labels(spd));
    spNN = spNN(spNN ~= I);
    Result{I, 1} = I;
    Result{I, 2} = size(spNN, 1);
    Result{I, 3} = spNN;
end

end

