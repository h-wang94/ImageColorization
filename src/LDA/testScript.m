%    This script aims to domenstrate FDA for 2-class 3-D data
%    Number of samples are 10
%
%    (c) Sultan Alzahrani, PhD Student, Arizona State University.
%    ssalzahr@asu.edu,  http://www.public.asu.edu/~ssalzahr/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
clear
x_base= [ 1,1,4
    2,0.5,.5
    3,-0.5,2
    4,2,-1
    5,4,1
    ];
f1 = figure();
% We try to create 10 sample points that linearly seperable, and
% intuitively do reduction and see results
% x_1 are points for class one (red dots)
x_1 = x_base-3;
% x_1 are points for class two (blue dots)
x_2 = x_base+3;

%we combined all datapoints to a single set where first 5 points are
%corresponding to class 1 and the rest for class 2.
x=[x_1;x_2];
% hold on
plot3(x_1(:,1),x_1(:,2),x_1(:,3),'r*')
hold on
plot3(x_2(:,1),x_2(:,2),x_2(:,3),'b*')

grid on
xlabel('X');
ylabel('Y');
zlabel('Z');


% gererating y
% first 5 points from class 1
% last 5 ponits from class 2
% we gave class label '1' to class 1
% we gave class label '2' to class 2

%y=[1,1,1,1,1,2,2,2,2,2];
y=[1*ones(1,5) 2*ones(1,5)];

X=x';
Y=y;
r=2;
% Now we do dimentionality reduction to 2D
[Z,W] = FDA(X,Y',r)
%plot figure as a one dirmentional
f2 = figure();
plot(Z(1,Y==1),Z(2,Y==1),'r*')
hold on
plot(Z(1,Y==2),Z(2,Y==2),'b*')


X=x';
Y=y;
r=1;
% Now we do dimentionality reduction to 1D
[Z,W] = FDA(X,Y',r)
%plot figure as a one dirmentional
f2 = figure();
plot(Z(Y==1),0,'r*')
hold on
plot(Z(Y==2),0,'b*')
%Now W holds the new datapoint in 1-D space.
% Again first 5 points are
%corresponding to class 1 and the rest for class 2.

% If there are more data to test, you can simply mutliply 
%W' by X to reduce dimentionality.
%Example:

%f3=figure()
%plot(W'*X,0,'*')




