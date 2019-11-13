function [count, mu_new, gamma_new] = Minimax(gamma_old, mu_old,xlength, ylength)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
clc; clear all; close all;
%{
The following code implements implicit twisted gradient descent; We seek to
find the saddle point of the Lagrangian L(z). The code terminates when the
magnitude of the gradient passes below a certain threshold E, or a
prescribed number of steps are preformed. The stepsize between successive
approximations is mu, which is dynamically updated to accelerate
convergence. The domain z is split into two components; x in which 
%}


%% Parameters

E = 0.0001; % Error threshold
growthfactor = 2; % Determines the growthrate of mu 
mu_max =100; % Maximum permissible value of mu
maxstep = 50; % Maximum number of itterations
count=1;

figure(1)
hold on

%% Iteration
while norm(G(z_old,xlength,ylength))>E || count <= maxstep
    mu_new = min(growthfactor*mu_old ,mu_max); % find new stepsize
    eta = mu_new/ norm(G(z_old,xlength,ylength)); 
    z_new = zstep( G(z_old,xlength,ylength),H(z_old,xlength,ylength),z_old,xlength,ylength,eta ); % find new approximation
    %% Decompose z
    x_new = z_new(1:xlength); 
    x_old = z_old(1:xlength); 
    y_new = z_new(1+xlength:xlength+ylength); 
    y_old = z_old(1+xlength:xlength+ylength);
    %% Acceptance test
    while (L([x_new;y_old],xlength,ylength) >= L([x_new;y_new],xlength,ylength) || L([x_new;y_new],xlength,ylength) >= L([x_old;y_new],xlength,ylength)) && mu_new>0.1
        mu_new = mu_new/2;
        eta = mu_new/ norm(G(z_old,xlength,ylength));
        gamma_new = zstep( G(z_old,xlength,ylength),H(z_old,xlength,ylength),z_old,xlength,ylength,eta );
    end
    %% Update
    count=count+1;
    mu_old = mu_new;
    z_old = z_new;
end
end

