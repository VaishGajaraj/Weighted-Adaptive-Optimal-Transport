function [ G,BG] = g1( x,b,sigma )
%% Unpacking b
[N,~]=size(b); K=(N-6)/3;
M=[b(1) b(2); b(2) b(3)]; v=[b(4); b(5)]; C=b(6); %Coefficients of the quadratic terms
mag=b(7:3:end); %Magnitude of the gaussian terms
centre=[b(8:3:end) b(9:3:end)]; %Centres of the gaussian terms
BG = zeros(N,1);
d1 = zeros(N,1); d2 = zeros(N-1,1); d3 = zeros(N-2,1);

%% Constructing G
G = x'*M*x + v'*x + C ;

%% Constructing BG
BG(1:6) = [x(1)^2; 2*x(1)*x(2); x(2)^2; x(1); x(2); 1];

%% Constructing BBG
% BBG=zeros(N);
end

