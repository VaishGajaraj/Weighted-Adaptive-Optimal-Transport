function [ l ] = L_weights( a,b,x,mu,y,nu,sigma )
%{
This function computes the Lagrangian of the minimax problem. The function
takes the arguments a,b. The parameters are (x,mu) and (y,nu) which
represent the pixel locations and their respective intensity.
%}
%% Setup
[~,xlength]=size(x);
[~,ylength]=size(y);
s1 = zeros(1,xlength); %Summands of first sum
s2 = zeros(1,ylength); %Summands of second sum

%% Compute Summands
for j=1:xlength
    [ XF,~,~] = f_linear( x(:,j),a,sigma );
    [ G,~ ] = g1( XF,b,sigma );
    s1(j) = G;
end
for j=1:ylength
    [ G,~ ] = g1( y(:,j),b,sigma );
    s2(j) = exp(G);
end

%% Compute l
l = s1*mu - s2*nu;
end

