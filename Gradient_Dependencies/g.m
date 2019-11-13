function [ XG,XGx,XGy,BG,BXG,BBG] = g( x,b,sigma )
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

one = ones(K,1);
gauss = exp(-((x(1)*one-centre(:,1)).^2 + (x(2)*one-centre(:,2)).^2) / (2*sigma^2)); % Gaussian terms
Gx = 2*x(1)* b(1)+2*x(2)*b(2)+b(4)+(mag.*(centre(:,1)-x(1)*one)')*gauss / (sigma^2);
Gy = 2*x(1)*b(2)+2*x(2)*b(3)+b(5)+(mag.*(centre(:,2)-x(2)*one)')*gauss / (sigma ^ 2);
Gxx = 2*b(1)+(mag.*((centre(:,1)-x(1)*one)'.^2/(sigma^2) - one'))*gauss / (sigma^2);
Gxy = 2*b(2)+(mag.*(centre(:,1)-x(1)*one)'.*(centre(:,2)-x(2)*one)')*gauss / (sigma^4);
Gyy = 2*b(3)+(mag.*((centre(:,2)-x(2)*one)'.^2 / (sigma^2) - one'))*gauss / (sigma^2);

%% Construct XG,XGx,XGy
XG=[Gx; Gy]; 
XGx=[Gxx; Gxy]; 
XGy=[Gxy; Gyy];

%% Construct BXG
BXG(:,1:6) = [2*x(1) 2*x(2) 0 1 0 0; 0 2*x(1) 2*x(2) 0 1 0];
BXG(:,7:3:N) = [ ((centre(:,1)-x(1)*one).*gauss/sigma^2)'; ((centre(:,2)-x(2)*one).*gauss/sigma^2)'];
BXG(:,8:3:N) = [ (-mag'.*((centre(:,1)-x(1)*one).^2 - sigma^2*one).*gauss/sigma^4)'; (-mag'.*(centre(:,1)-x(1)*one).*(centre(:,2)-x(2)).*gauss/sigma^4)'];
BXG(:,9:3:N) = [ (-mag'.*(centre(:,1)-x(1)*one).*(centre(:,2)-x(2)).*gauss/sigma^4)'; (-mag'.*((centre(:,2)-x(2)*one).^2 - sigma^2*one).*gauss/sigma^4)'];

%% Constructing BG
BG(1:6) = [x(1)^2; 2*x(1)*x(2); x(2)^2; x(1); x(2); 1];
BG(7:3:N) = gauss;
BG(8:3:N) = -mag'.*(centre(:,1)-x(1)*one).*gauss / sigma^2;
BG(9:3:N) = -mag'.*(centre(:,2)-x(2)*one).*gauss / sigma^2;

%% Constructing BBG
d1(8:3:N) = mag'.*((centre(:,1)-x(1)*one).^2 - (sigma^2)*one).*gauss / sigma^4;
d1(9:3:N) = mag'.*((centre(:,2)-x(2)*one).^2 - (sigma^2)*one).*gauss / sigma^4;
d2(7:3:N-1) = -(centre(:,1)-x(1)*one).*gauss / sigma^2;
d2(8:3:N-1) = mag'.*(centre(:,1)-x(1)*one).*(centre(:,2)-x(2)*one).*gauss / sigma^4;
d3(7:3:N-2) = -(centre(:,2)-x(2)*one).*gauss / sigma^2;
BBG = diag(d3,-2)+diag(d2,-1)+diag(d1,0)+diag(d2,1)+diag(d3,2);
end