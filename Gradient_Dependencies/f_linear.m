function [ XF,AFx,AFy ] = f_linear( x,a,sigma )
%{
This family of functions describes the space over which we seek to minimize
L. This function is 2 dimensional, with a N=5+3K dimensional parameter
space. The function consists of a quadratic part and a gaussian part.
%}

%% Unpacking a
[N,~]=size(a); K = (N-5)/3;
M=[a(1) a(2); a(2) a(3)]; v=[a(4); a(5)]; %Coefficients of the quadratic terms
mag=a(6:3:end); %Magnitude of the gaussian terms
centre=[a(7:3:end) a(8:3:end)]; %Centres of the gaussian terms

%% 
AFx = zeros(N,1); AFy = zeros(N,1);
d1 = zeros(N,1); d2 = zeros(N-1,1); d3 = zeros(N-2,1);


%% Compute XF
XF = 2*M*x + v;

%% Construct AFx
AFx(1:5) = [ 2*x(1); 2*x(2); 0; 1; 0];

%% Construct AFy
AFy(1:5) = [0; 2*x(1); 2*x(2); 0; 1];

% %% Construct AAFx
% AAFx=zeros(N);
% %% Construct AAfy
% AAFy=zeros(N);
end