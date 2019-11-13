function [ XF,AFx,AFy,AAFx,AAFy ] = f( x,a,sigma )
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
one = ones(K,1);
gauss = exp(-((x(1)*one-centre(:,1)).^2 + (x(2)*one-centre(:,2)).^2)/(2*sigma^2));

%% Compute XF
XF = 2*M*x + v + ((centre' - x*one').*[gauss'; gauss'])*mag / sigma^2;

%% Construct AFx
AFx(1:5) = [ 2*x(1); 2*x(2); 0; 1; 0];
AFx(6:3:N) = (centre(:,1) - x(1)*one).*gauss / (sigma^2);
AFx(7:3:N) = -mag.*((centre(:,1)-x(1)*one).^2 - (sigma^2)*one).*gauss / sigma^4;
AFx(8:3:N) = -mag.*(centre(:,1)-x(1)*one).*(centre(:,2)-x(2)*one).*gauss / sigma^4;

%% Construct AFy
AFy(1:5) = [0; 2*x(1); 2*x(2); 0; 1];
AFy(6:3:N) = (centre(:,2)-x(2)*one).*gauss / (sigma^2);
AFy(7:3:N) = -mag.*(centre(:,1)-x(1)*one).*(centre(:,2)-x(2)*one).*gauss / sigma^4;
AFy(8:3:N) = -mag.*((centre(:,2)-x(2)*one).^2 - (sigma^2)*one).*gauss / sigma^4;

%% Construct AAFx
d1(7:3:N) = mag.*(centre(:,1)-x(1)*one).*(centre(:,1).^2 - 2*x(1)*centre(:,1) + (x(1)^2 -3*sigma^2)*one).*gauss / sigma^6;
d1(8:3:N) = mag.*(centre(:,1)-x(1)*one).*((centre(:,2)-x(2)*one).^2 - (sigma^2)*one).* gauss / sigma^6;
d2(6:3:N-1) = -((centre(:,1)-x(1)*one).^2 - (sigma^2)*one).*gauss / sigma^4;
d2(7:3:N-1) = mag.*(centre(:,2)-x(2)*one).*((centre(:,1)-x(1)*one).^2 - (sigma^2)*one).* gauss / sigma^6;
d3(6:3:N-2) = -(centre(:,1)-x(1)*one).*(centre(:,2)-x(2)*one).*gauss / sigma^4;
AAFx = diag(d3,-2)+diag(d2,-1)+diag(d1,0)+diag(d2,1)+diag(d3,2);

%% Construct AAfy
d1(7:3:N) = mag.*(centre(:,2)-x(2)*one).*((centre(:,1)-x(1)*one).^2 - (sigma^2)*one).* gauss / sigma^6;
d1(8:3:N) = mag.*(centre(:,2)-x(2)*one).*(centre(:,2).^2 - 2*x(2)*centre(:,2) + (x(2)^2 -3*sigma^2)*one).*gauss / sigma^6;
d2(6:3:N-1) = -(centre(:,1)-x(1)*one).*(centre(:,2)-x(2)*one).*gauss / sigma^4;
d2(7:3:N-1) = mag.*(centre(:,1)-x(1)*one).*((centre(:,2)-x(2)*one).^2 - (sigma^2)*one).* gauss / sigma^6;
d3(6:3:N-2) = -((centre(:,2)-x(2)*one).^2 - (sigma^2)*one).*gauss / sigma^4;
AAFy = diag(d3,-2)+diag(d2,-1)+diag(d1,0)+diag(d2,1)+diag(d3,2);
end