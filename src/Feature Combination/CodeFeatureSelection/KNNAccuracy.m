function f=KNNAccuracy(ldata,tdata,K)

% FUNCTION DESCRIPTION %
% This function is to compute the classification accuracy by KNN
%
% INPUTS %
% 'ldata': Learning data of size L x m+1, L being the number of data points
% for learning and m being the number of features per data point. The last
% column contains the labels for each data point.
% 'tdata': Test data of size T x m+1 , T being the number of test data
% points. The last column contains the labels for each data point, to
% compute the classification accuracy by comparing the predicted and actual
% labels
%
% OUTPUTS %
% 'f': It is the classification accuracy by KNN
%
% Written By: Sujoy Paul, Jadavpur University, India
% Email: paul.sujoy.ju@gmail.com


% Initialization
sl = size(ldata);
st = size(tdata);
pre_act=[];
p=0;

%More efficient using knnsearch/pdist.
v = zeros(st(1), sl(1));
for i=1:st(1)
  for j=1:sl(1)
    v(i,j)=norm(tdata(i,1:st(2)-1)-ldata(j,1:sl(2)-1));
  end
  [v(i,:), index]=sort(v(i,:));
  r=ldata(index,sl(2));
  m=mode(r(1:K));      %Predicted class
  pre_act(i,:)=[m tdata(i,st(2))];
  if m==tdata(i,st(2))
    p=p+1;
  end
end

f=p/st(1);


