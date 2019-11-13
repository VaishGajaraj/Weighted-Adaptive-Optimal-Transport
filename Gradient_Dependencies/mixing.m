function [z_interpolated] = mixing(xlength,ylength,z_xi, z_yj,N,k)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
a = randperm(xlength, floor(xlength*((N-k)/(N))));
b = randperm(ylength, floor(ylength*(k/N)));


for j = 1:N

    x = z_xi(a,:);
    y = z_yj(b,:);
    
end
z_interpolated = [x;y];

end

