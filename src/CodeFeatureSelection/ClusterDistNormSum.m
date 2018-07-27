function f = ClusterDistNormSum(ldata, mcCost)
%TODO: Change function to implement the median idea.

% Initialization
intra=0;
inter=0;

% Computation
for i=1:size(ldata,1)
    for j=i+1:size(ldata,1)
        D = norm(ldata(i,1:end-1)-ldata(j,1:end-1));
        
        if (ldata(i,end) ~= -1 && ldata(j,end) ~= -1)
          if ldata(i,end) == ldata(j,end)
              intra = intra + D;
          else
              inter = inter + mcCost(ldata(i,end),ldata(j,end))*D;
          end
        end
    end
end

f = [intra;inter];