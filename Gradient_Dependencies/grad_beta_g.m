function [grad_beta_g] = grad_beta_g(z, b1,b2,b3,b4,b5,b6,b7,b8,b9, b10,b11,b12, sigma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    y0=  z(:,1);
    y1=  z(:,2);
    
    arg1 = (-1/(sigma*2)).*((y0-b8).^2+(y1-b9).^2);
    arg2 = (-1/(sigma*2)).*((y0-b11).^2+(y1-b12).^2);
    
    b_8= b7*(-1/(sigma*2))*(-2*y0*b8+2*b8);
    b_9= b7*(-1/(sigma*2))*(-2*y1*b9+2*b9);
    b_11= b10*(-1/(sigma*2))*(-2*y0*b11+2*b11);
    b_12= b10*(-1/(sigma*2))*(-2*y1*b12+2*b12);
    
    grad_beta_g= [zeros(size(y0)) y0 y1 .5*y0.^2 y0.*y1 .5*y1.^2 exp(arg1) b_8 b_9 exp(arg2) b_11 b_12];
       
    
end