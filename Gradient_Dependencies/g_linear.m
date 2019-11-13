function [ XG,BG] = g_linear( x,b,sigma )
%{
This family of functions describes the space over which we seek to maximize
L. This function is 2 dimensional, with a N=6+3K dimensional parameter
space. The function consists of a quadratic part and a gaussian part.
%}

%% Unpacking b
[N,~]=size(b); K=(N-6)/3;
M=[b(1) b(2); b(2) b(3)]; v=[b(4); b(5)]; C=b(6); %Coefficients of the quadratic terms
mag=b(7:3:end)'; %Magnitude of the gaussian terms
centre=[b(8:3:end) b(9:3:end)]; %Centres of the gaussian terms
BXG = zeros(2,N); BG = zeros(N,1);
d1 = zeros(N,1); d2 = zeros(N-1,1); d3 = zeros(N-2,1);
Gx = 2*x(1,:)* b(1)+2*x(2,:)*b(2)+b(4);
Gy = 2*x(1,:)*b(2)+2*x(2,:)*b(3)+b(5);
% Gxx = 2*b(1);
% Gxy = 2*b(2);
% Gyy = 2*b(3);

%% Construct XG,XGx,XGy
XG=[Gx; Gy];

%% Construct BXG
% BXG(:,1:6) = [2*x(1) 2*x(2) 0 1 0 0; 0 2*x(1) 2*x(2) 0 1 0];


%% Constructing BG
BG(1:6) = [x(1)^2; 2*x(1)*x(2); x(2)^2; x(1); x(2); 1];

%% Constructing BBG
% BBG=zeros(N);

end