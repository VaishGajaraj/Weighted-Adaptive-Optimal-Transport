

clc;
close all;


mu1 = [2 3];
sigma1 = [1 1.5; 1.5 3];

mu2 = [1 4];
sigma2 = [1.5 2; 2 3.5];

xlength = 100;
ylength= 100;

%rng('default')  % For reproducibility
z_xi = mvnrnd(mu1,sigma1,ylength);

figure('Name','Original Data','NumberTitle','off');
plot(z_xi(:,1),z_xi(:,2),'b+');

hold on
%rng('default')  % For reproducibility
z_yj = mvnrnd(mu2,sigma2,xlength);
 
plot(z_yj(:,1),z_yj(:,2),'r*')

E = 0.1; % Error threshold
growthfactor = 2; % Determines the growthrate of mu 
mu_max =100; % Maximum permissible value of mu
maxstep = 10; % Maximum number of itterations
count=1;
sigma = 1;



mu_old = .005;
%gamma_initial= [.5 .5 .5 0 .5 0 1 .5 .5 1 0 .5 1 0 1 0 1 0 1 0];
gamma_initial = [0 .5 0 0 .5 0 .4 .6 .3 .5 0 .5 0 .5 0 1 .5 .5 0 0];
%gamma_initial(1), gamma_initial(2),gamma_initial(3),gamma_initial(4), gamma_initial(5),gamma_initial(6),gamma_initial(7),gamma_initial(8),gamma_initial(9),gamma_initial(10),gamma_initial(11),gamma_initial(12),gamma_initial(13),gamma_initial(14),gamma_initial(15),gamma_initial(16),gamma_initial(17)
array_of_ones = [ones(8,1)' (-1.*ones(12,1))' ];

J = diag(array_of_ones);
B_Quasi = J;
G_initial = G(z_xi,z_yj, gamma_initial(1), gamma_initial(2),gamma_initial(3),gamma_initial(4), gamma_initial(5),gamma_initial(6),gamma_initial(7),gamma_initial(8),gamma_initial(9),gamma_initial(10),gamma_initial(11),gamma_initial(12),gamma_initial(13),gamma_initial(14),gamma_initial(15),gamma_initial(16),gamma_initial(17),gamma_initial(18),gamma_initial(19),gamma_initial(20), sigma, xlength, ylength);

eta = mu_old/ norm(G_initial);

gamma_one = gamma_initial' + eta*B_Quasi*G_initial;
%gamma_next = gamma_one;
%% Iterationorm(G(z_xi,z_yj, gamma_next(1), gamma_next(2),gamma_next(3),gamma_next(4), gamma_next(5),gamma_next(6),gamma_next(7),gamma_next(8),gamma_next(9),gamma_next(10),gamma_next(11),gamma_next(12),gamma_next(13),gamma_next(14),gamma_next(15),gamma_next(16),gamma_next(17),gamma_next(18),gamma_next(19),gamma_next(20), sigma, xlength, ylength))>E ||
while  count <= maxstep
    count
    mu_new = min(growthfactor*mu_old ,mu_max); % find new stepsize
    G_next = G(z_xi,z_yj, gamma_next(1), gamma_next(2),gamma_next(3),gamma_next(4), gamma_next(5),gamma_next(6),gamma_next(7),gamma_next(8),gamma_next(9),gamma_next(10),gamma_next(11),gamma_next(12),gamma_next(13),gamma_next(14),gamma_next(15),gamma_next(16),gamma_next(17),gamma_next(18),gamma_next(19),gamma_next(20), sigma, xlength,ylength );
    
    eta = mu_new/ norm(G_next);
    
    B_Quasi = B_Updated(G_initial, G_next,B_Quasi,J)
    
    gamma_next = gamma_next + eta*B_Quasi*G_next;
    if isnan(G_next)
        fprintf("Illegal gradient value. Code is updating learning rate:")
    end
    
    
    while L(z_xi, z_yj, xlength,ylength,gamma_next(1), gamma_next(2),gamma_next(3),gamma_next(4), gamma_next(5),gamma_next(6),gamma_next(7),gamma_next(8),gamma_initial(9),gamma_initial(10),gamma_initial(11),gamma_initial(12),gamma_initial(13),gamma_initial(14),gamma_initial(15),gamma_initial(16),gamma_initial(17),gamma_initial(18),gamma_initial(19),gamma_initial(20), sigma)>= L(z_xi, z_yj, xlength,ylength,gamma_next(1), gamma_next(2),gamma_next(3),gamma_next(4), gamma_next(5),gamma_next(6),gamma_next(7),gamma_next(8),gamma_next(9),gamma_next(10),gamma_next(11),gamma_next(12),gamma_next(13),gamma_next(14),gamma_next(15),gamma_next(16),gamma_next(17),gamma_next(18),gamma_next(19),gamma_next(20), sigma) && L(z_xi, z_yj, xlength,ylength,gamma_next(1), gamma_next(2),gamma_next(3),gamma_next(4), gamma_next(5),gamma_next(6),gamma_next(7),gamma_next(8),gamma_next(9),gamma_next(10),gamma_next(11),gamma_next(12),gamma_next(13),gamma_next(14),gamma_next(15),gamma_next(16),gamma_next(17),gamma_next(18),gamma_next(19),gamma_next(20), sigma) >= L(z_xi, z_yj, xlength,ylength,gamma_initial(1), gamma_initial(2),gamma_initial(3),gamma_initial(4), gamma_initial(5),gamma_initial(6),gamma_initial(7),gamma_initial(8),gamma_next(9),gamma_next(10),gamma_next(11),gamma_next(12),gamma_next(13),gamma_next(14),gamma_next(15),gamma_next(16),gamma_next(17),gamma_next(18),gamma_next(19),gamma_next(20), sigma)
        mu_new = mu_new/2;
        eta = mu_new/ norm(G_next);
        gamma_next = gamma_next - eta*B_Quasi*G_next;
        
        B_Quasi = B_Updated(G_initial, G_next,B_Quasi,J);
                
        gamma_next = gamma_next - eta*B_Quasi*G_next;
        %G_initial = G_next;
        if isnan(G_next)
            fprintf("debug:")
        end
        
        z_xi = phi_grad(z_xi,gamma_next(1), gamma_next(2),gamma_next(3),gamma_next(4), gamma_next(5),gamma_next(6),gamma_next(7),gamma_next(8),sigma);

        figure

        plot(z_xi(:,1),z_xi(:,2),'+')
        hold on

        plot(z_yj(:,1),z_yj(:,2),'r*');
        
        if count >= maxstep
            
            break;
        end
        
    
    end
    
        
    
    

    
    z_xi = phi_grad(z_xi,gamma_next(1), gamma_next(2),gamma_next(3),gamma_next(4), gamma_next(5),gamma_next(6),gamma_next(7),gamma_next(8),sigma);
    
    figure('Name','Transported Data')
    

    plot(z_xi(:,1),z_xi(:,2),'+');
    hold on
    plot(z_yj(:,1),z_yj(:,2),'r*')
    drawnow
    
    count = count +1;
    if count >=maxstep
        break;
    end


    
end
