function [ z_new ] = localOT( x,y,alength,blength,E )
%{
The following code implements implicit twisted gradient descent; We seek to
find the saddle point of the Lagrangian L(z). The code terminates when the
magnitude of the gradient passes below a certain threshold E, or a
prescribed number of steps are preformed. The stepsize between successive
approximations is mu, which is dynamically updated to accelerate
convergence. The domain z is split into two components; x in which 
%}

%% Parameters
sigma = 0.1;
[~,xlength] = size(x);
[~,ylength] = size(y);
mu = ones(xlength,1)/xlength;
nu = ones(ylength,1)/ylength;
mu_old = 0.5; % Initial step size
growthfactor = 2; % Determines the growthrate of mu 
mu_max =10000; % Maximum permissible value of mu
maxstep = 100000; % Maximum number of itterations
count=1;
J = diag([ones(alength,1);-ones(blength,1)]);
z_old = zeros(alength+blength,1);

x_bar=mean(x,2);
y_bar=mean(y,2);
sigma_x = cov(x');
sigma_y = cov(y');
det_x = det(sigma_x);
det_y = det(sigma_y);
inv_x=inv(cov(x'));
inv_y=inv(cov(y'));

z_old(1)=0.5; z_old(3)=0.5;
C= (inv_y-inv_x)/2;
z_old(alength+1)=C(1,1); z_old(alength+2)=C(1,2); z_old(alength+3)=C(2,2);
C= inv_x*x_bar - inv_y*y_bar;
z_old(alength+4:alength+5)=inv_y*y_bar - inv_x*x_bar;
z_old(alength+6)=(log(det_y/det_x) + y_bar'*sigma_y*y_bar - x_bar'*sigma_x*x_bar)/2;

[G_old] = gh(z_old(1:alength),z_old(1+alength:alength+blength),x,mu,y,nu,sigma);

%% Iteration
while norm(G_old)>E && count <= maxstep
    
    
    p = zeros(size(x));
    for j=1:xlength
        [ XF,~,~,~,~ ] = f( x(:,j),z_old(1:alength),sigma );
        p(:,j)=XF;
    end
    scatter(p(1,:),p(2,:),'b*');
    hold on
    scatter(y(1,:),y(2,:),'r*');
    drawnow
    hold off
    
    mu_new = min(growthfactor*mu_old ,mu_max); % find new stepsize
    eta = mu_new/ norm(G);
    z_new =  z_old - eta*((J+eta*H)\G); % find new approximation
%     %% Decompose z
%     a_new = z_new(1:alength); 
%     a_old = z_old(1:alength); 
%     b_new = z_new(1+alength:alength+blength); 
%     b_old = z_old(1+alength:alength+blength);
%     %% Acceptance test
%    while (L(a_new,b_old,x,mu,y,nu,sigma) > L(a_new,b_new,x,mu,y,nu,sigma) || L(a_new,b_new,x,mu,y,nu,sigma) > L(a_old,b_new,x,mu,y,nu,sigma)) && mu_new>0.01
%         mu_new = mu_new/2;
%         eta = mu_new/ norm(G);
%         z_new = z_old - eta*((J+eta*H)\G);
%     end
    %% Update
    count=count+1;
    mu_old = mu_new;
    z_old = z_new
    [G,H] = gh( z_old(1:alength),z_old(1+alength:alength+blength),x,mu,y,nu,sigma );
    norm(G);
end
end