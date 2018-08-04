function f=clusterDist(ldata, mcCost)

% FUNCTION DESCRIPTION %
% This function is to compute the inter and intra class distance vector of
% the dataset without feature selection or weighting
%
% INPUT %
% 'ldata': It is the learning data. It is a matrix of size n x m+1, where n
% being the number of data points and m being the number of features per
% data point. The last column contains the labels for each data point.
%
% OUTPUTS %
% 'f' is a matrix of size 2 x m, s.t. f(1,:) contain the intra class
% distance vector and f(2,:) contain the inter class distance vector
%
% Written By: Sujoy Paul, Jadavpur University, India
% Email: paul.sujoy.ju@gmail.com

% Initialization
intra=zeros(1,size(ldata,2)-1);
inter=zeros(1,size(ldata,2)-1);

% Computation
for i=1:size(ldata,1)
    for j=i+1:size(ldata,1)
        D = ldata(i,1:end-1)-ldata(j,1:end-1);
        D = abs(D);
        if (ldata(i,end) == ldata(j,end))
            intra = intra + D;
        else
            inter = inter + mcCost(ldata(i,end), ldata(j,end))*D;
        end
    end
end

f=[intra;inter];