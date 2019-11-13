function [p] = Local_transport_function(x,y,alength,blength)
% %UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%weights
sigma = 1;
[~,xlength] = size(x);
[~,ylength] = size(y);
mu = ones(xlength,1)/xlength;
nu = ones(ylength,1)/ylength;
mu_old = 0.05; % Initial step size
growthfactor = 1.5; % Determines the growthrate of mu 
mu_max =100; % Maximum permissible value of mu
maxstep = 3; % Maximum number of itterations
count=1;
J = diag([ones(alength,1);-ones(blength,1)]);
z_old = zeros(alength+blength,1);
E=.05;

%% Initial Conditions



x_bar=mean(x,2);
y_bar=mean(y,2);
sigma_x = cov(x');
sigma_y = cov(y');
det_x = det(sigma_x);
det_y = det(sigma_y);
inv_x=inv(sigma_x);
inv_y=inv(sigma_y);

z_old(1)=0.5; 
z_old(3)=0.5; 
% z_old(7) = 3;
C= (inv_y-inv_x)/2;
z_old(alength+1)=C(1,1); 
z_old(alength+2)=C(1,2); 
z_old(alength+3)=C(2,2);
C= inv_x*x_bar - inv_y*y_bar;
z_old(alength+4:alength+5)=inv_y*y_bar - inv_x*x_bar;
z_old(alength+6)=(log(det_y/det_x) + y_bar'*sigma_y*y_bar - x_bar'*sigma_x*x_bar)/2;
number_of_clusters = (blength-6) / 3;

%% placing gaussians using K-Means
% [idx,Cx] = kmeans(x',number_of_clusters);

[idy,Cy] = kmeans(y',number_of_clusters);
for i = 0:number_of_clusters-1
    z_old(7+(3*i)) = Cy(i+1,1);
    z_old(8+(3*i)) = Cy(i+1,2);
    
%     z_old(alength+7+(3*i)) = .5;

    z_old(alength+8+(3*i)) = Cy(i+1,1);
    z_old(alength+9+(3*i)) = Cy(i+1,2);
end

[G,~] = gh( z_old(1:alength),z_old(1+alength:alength+blength),x,mu,y,nu,sigma );

B=J;
%% Iteration
while norm(G)>E 
%     && maxstep>count
    
    
    p = zeros(size(x));
    for j=1:xlength
        [ XF,~,~,~,~ ] = f( x(:,j),z_old(1:alength),sigma );
        p(:,j)=XF;
    end
%     scatter(p(1,:),p(2,:),'b*');
%     hold on
%     scatter(y(1,:),y(2,:),'r*');
%     for i = 0:number_of_clusters-1
% 
%         plot(Cx(i+1,1),Cx(i+1,2), 'bx','MarkerSize',15,'LineWidth',3);
%         plot(Cy(i+1,1),Cy(i+1,2), 'kx','MarkerSize',15,'LineWidth',3);
% 
%     end
%     drawnow
%     hold off
    mu_new = mu_old;
    mu_new = min(growthfactor*mu_old ,mu_max); % find new stepsize

    eta= mu_new/norm(G);
    
    z_new = z_old - eta*B*(G);



    % Decompose z
    a_new = z_new(1:alength); 
    a_old = z_old(1:alength); 
    b_new = z_new(1+alength:alength+blength); 
    b_old = z_old(1+alength:alength+blength);
    % Acceptance test
   while (L_weights(a_new,b_old,x,mu,y,nu,sigma) > L_weights(a_new,b_new,x,mu,y,nu,sigma) || L_weights(a_new,b_new,x,mu,y,nu,sigma) > L_weights(a_old,b_new,x,mu,y,nu,sigma)) && mu_new>0.001
        mu_new = mu_new/2;
        eta = mu_new/ norm(G);
        z_new = z_old - eta*(B*G);
   end
   
   count = count+1;
   mu_old = mu_new;
   z_old= z_new;
   G_store = G;

   [G,~] = gh( z_old(1:alength),z_old(1+alength:alength+blength),x,mu,y,nu,sigma );

   B = B_Updated(G_store, G, B, J);




end
z_old

p = zeros(size(x));
    for j=1:xlength
        [ XF,~,~,~,~ ] = f( x(:,j),z_old(1:alength),sigma );
        p(:,j)=XF;
    end
    


end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

