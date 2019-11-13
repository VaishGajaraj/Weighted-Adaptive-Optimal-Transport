function [g] = highlighter_g(z, b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12, sigma)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    
    y0=  z(:,1);
    y1=  z(:,2);
    
    
    arg1 = (-1/(sigma*2)).*((y0-b8).^2+(y1-b9).^2);
    arg2 = (-1/(sigma*2)).*((y0-b11).^2+(y1-b12).^2);
    g = b1+ b2.*y0+b3.*y1+0.5.*b4.*y0.^2+b5.*y0.*y1+0.5*b6*y1.^2 + b7.*exp(arg1) +b10.*exp(arg2);
    
    
  
end

