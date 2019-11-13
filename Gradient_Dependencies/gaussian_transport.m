function [ A ] = gaussian_transport( z_xi,z_yj, xlength, mean_one, ylength, mean_two)

%Closed form solution to transport between two gaussians (sigma_x^-.5(sigma_x^.5*sigma_y*sigma_x^.5)*sigma_x^-.5)
    x=  z_xi(:,1);
    y=  z_xi(:,2);
    
    x_mean = mean_one(:,1);
    y_mean = mean_one(:,2);
    
    s11_x = 1/xlength *sum( (x-x_mean).^2);
    
    s12_x = 1/xlength *sum( (x-x_mean).*(y-y_mean));
    
    s22_x = 1/xlength * sum((y-y_mean).^2);
    
    sigma_x = [s11_x s12_x; s12_x s22_x];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% first covariance matrix
    
    a=  z_yj(:,1);
    b=  z_yj(:,2);
    
    a_mean = mean_two(:,1);
    b_mean = mean_two(:,2);
    
    s11_y = 1/ylength *sum( (a-a_mean).^2);
    
    s12_y = 1/ylength *sum( (a-a_mean).*(b-b_mean));
    
    s22_y = 1/ylength * sum((b-b_mean).^2);
    
    sigma_y = [s11_y s12_y; s12_y s22_y];
    
    A =inv(sqrtm(sigma_x))*sqrtm((sqrtm(sigma_x)*sigma_y*sqrtm(sigma_x)))*(inv(sqrtm(sigma_x)));
    
    