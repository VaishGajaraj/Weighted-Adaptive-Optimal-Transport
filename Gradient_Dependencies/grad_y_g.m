function [grad_y_g] = grad_y_g(z, b1,b2,b3,b4,b5,b6,b7,b8,b9 ,b10,b11,b12, sigma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    y0=  z(:,1);
    y1=  z(:,2);
    
    arg1 = (-1/(sigma*2))*((y0-b8).^2+(y1-b9).^2);
    arg2 = (-1/(sigma*2))*((y0-b11).^2+(y1-b12).^2);
    
    first_row = b1+ b4*y0+b5*y1+b7*(-1/(sigma*2))*(2*y0 - 2*b8).*exp(arg1)+b10*(2*y0 - 2*b11).*exp(arg2);
    
    second_row = b3+b5*y1+b6*y1+b7*(-1/(sigma*2))*(2*y1 - 2*b9).*exp(arg1)+b10*(2*y1 - 2*b12).*exp(arg2);
    
    grad_y_g = [first_row second_row];
    
    
    
end


