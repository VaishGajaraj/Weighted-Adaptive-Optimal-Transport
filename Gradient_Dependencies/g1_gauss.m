function [ G,BG,BBG ] = g1_gauss( x,b,sigma )
%% Unpacking b
[N,~]=size(b); K=(N-6)/3;
M=[b(1) b(2); b(2) b(3)]; v=[b(4); b(5)]; C=b(6); %Coefficients of the quadratic terms
mag=b(7:3:end); %Magnitude of the gaussian terms
centre=[b(8:3:end) b(9:3:end)]; %Centres of the gaussian terms
BG = zeros(N,1);
d1 = zeros(N,1); d2 = zeros(N-1,1); d3 = zeros(N-2,1);
one = ones(K,1);
gauss = exp(-((x(1)*one-centre(:,1)).^2 + (x(2)*one-centre(:,2)).^2)/(2*sigma^2));

%% Constructing G
G = x'*M*x + v'*x + C + mag'*gauss;

%% Constructing BG
BG(1:6) = [x(1)^2; 2*x(1)*x(2); x(2)^2; x(1); x(2); 1];
BG(7:3:N) = gauss;
BG(8:3:N) = -mag.*(centre(:,1)-x(1)*one).*gauss / sigma^2;
BG(9:3:N) = -mag.*(centre(:,2)-x(2)*one).*gauss / sigma^2;

%% Constructing BBG
d1(8:3:N) = mag.*((centre(:,1)-x(1)*one).^2 - (sigma^2)*one).*gauss / sigma^4;
d1(9:3:N) = mag.*((centre(:,2)-x(2)*one).^2 - (sigma^2)*one).*gauss / sigma^4;
d2(7:3:N-1) = -(centre(:,1)-x(1)*one).*gauss / sigma^2;
d2(8:3:N-1) = mag.*(centre(:,1)-x(1)*one).*(centre(:,2)-x(2)*one).*gauss / sigma^4;
d3(7:3:N-2) = -(centre(:,2)-x(2)*one).*gauss / sigma^2;
BBG = diag(d3,-2)+diag(d2,-1)+diag(d1,0)+diag(d2,1)+diag(d3,2);
end



