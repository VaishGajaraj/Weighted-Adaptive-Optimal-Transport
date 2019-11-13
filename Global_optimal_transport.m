clc;
clear;
close all;

mean_x = [0 0];
cov_x = [1 0;0 1];

% mean_y = [2 4];
% cov_y = [2 1;1 2];

xlength = 200;
ylength= 300;

rng(500)  % For reproducibility
z_xi = mvnrnd(mean_x,cov_x,xlength);
z_xi= z_xi';

figure('Name','Original Data','NumberTitle','off');
plot(z_xi(1,:),z_xi(2,:),'b+');

hold on
% rng(500)  % For reproducibility
%  r1 =2; 
%  r2 = 1;
%  r = .5*rand(ylength)+2*ones(ylength); % Using square root here ensures distribution uniformity by area
%  t = rand(1,ylength) *2*pi;
% x = r.*cos(t);
% y = r.*sin(t);
% z_xi = [x; x ];
% z_yj = [x ;y ];
% mean_y = mean(z_yj');
% cov_y = cov(z_yj');
%%%%%%%
rng(500)
x=rand(1,ylength)*5;
y=rand(1,ylength)*5;
z_yj = [x ;y ];
mean_y = mean(z_yj')
cov_y = cov(z_yj');
%%%%%%%%
plot(z_yj(2,:),z_yj(1,:),'r*')
% z_yj = mvnrnd(mean_y,cov_y,ylength);
% z_yj= z_yj';
% plot(z_xi(1,:),z_xi(1,:),'b*')

alength =32;
blength =33;
A = gaussian_transport(z_xi',z_yj', xlength, mean_x, ylength, mean_y);

x_bar = ones(xlength,1)*mean_x;

T_x = ((ones(xlength,1)*mean_y) +(z_xi'-x_bar)*A)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_star = gaussian_transport(z_yj', z_xi',ylength, mean_y, xlength, mean_x);

y_bar = ones(ylength,1).*mean_y;

T_y = ((ones(ylength,1).*mean_x)+(z_yj'-y_bar)*A_star)';

color =[0.1 0.6740 0.1880];
plot(T_x(1,:),T_x(2,:),'ko');
color =[0.4660 0.1 0.1880];
% scatter(T_y(1,:),T_y(2,:),30,color,'x');
hold off;

N=4;

cloud_forward = cell(xlength, 2) ;
cloud_backward = cell(ylength, 2) ;
a = randperm(xlength);
b = randperm(ylength);
%interpolating my two sets of transported points (T_X & T_Y)

for k = 0 : N
%     size((((N-k)/(N)))*z_xi(a(1:floor(xlength*((N-k)/(N)))),:))
%     size(((k/N)) *T_x(a(1:floor(xlength*((N-k)/(N)))),:))
    z0k = ((N-k)/(N))*z_xi(:,a(1:floor(xlength*((N-k)/N))))+((k/N)) *T_x(:,a(1:floor(xlength*((N-k)/N))));

    zNk = ((k/N)) *z_yj(:,b(1:floor(ylength*(k/N)))) +(((N-k)/(N)))*T_y(:,b(1:floor(ylength*(k/N))));

    
    cloud_interpolated{k+1} = [z0k  zNk];
figure
plot(cloud_interpolated{end}(1,:),cloud_interpolated{end}(2,:),'k*');
hold on
plot(z_xi(1,:),z_xi(2,:),'bo')
plot(z_yj(1,:),z_yj(2,:),'r*')
drawnow
hold off
end
% plot(cloud_interpolated{end}(1,:),cloud_interpolated{end}(2,:),'b+');
for k =1:N


        locally_transported_cloud = Local_transport_function( cloud_interpolated{k},cloud_interpolated{k+1},alength,blength);
        
        k
%         plot(locally_transported_cloud(1,:),locally_transported_cloud(2,:),'r*');

        cloud_interpolated{k+1} = locally_transported_cloud;
        
        
        

%         hold on
%         plot(cloud_interpolated{k}(1,:),cloud_interpolated{k}(2,:),'b+');
%         plot(locally_transported_cloud(1,:),locally_transported_cloud(2,:),'o');
        
%         
%         plot(z_xi(1,:),z_xi(2,:),'^')
%          plot(z_yj(1,:),z_yj(2,:),'*')
%         drawnow
%         hold off

end
% figure('Name','Transported Data (in blue)')
% plot(cloud_interpolated{end}(1,:),cloud_interpolated{end}(2,:),'b+');
% hold on
% plot(z_xi(1,:),z_xi(2,:),'o')
% plot(z_yj(1,:),z_yj(2,:),'*')
% hold off
index = randperm(xlength);

% interpolating my two sets of transported points (T_X & T_Y)

z_yj_transported = cloud_interpolated{end};
optimal_cloud_interpolated ={z_xi};
figure
for k = 1:100
        disp('Running through McCann Interpolation scheme. Iteration:')
        k
        for k = 0:N

                    interpolated_optimal_cloud = ((N-k)/(N))*z_xi +((k)/(N))*z_yj_transported;
            %         size(z0k)
            %         size(zNk)
                    optimal_cloud_interpolated{k+1} = interpolated_optimal_cloud;

%                     figure
%                     scatter(optimal_cloud_interpolated{k+1}(:,1),optimal_cloud_interpolated{k+1}(:,2),30,color,'o');
%                     hold on
%                     plot(z_xi(:,1),z_xi(:,2),'b+')
%                     plot(z_yj(:,1),z_yj(:,2),'r*')

        end

        for i=1:N
                optimal_cloud = Local_transport_function( optimal_cloud_interpolated{i},optimal_cloud_interpolated{i+1}, alength,blength);
                optimal_cloud_interpolated{i+1} = optimal_cloud;
%                 plot(optimal_cloud(1,:),optimal_cloud(2,:),'r+')
%                 hold on
%                 plot(optimal_cloud_interpolated{1}(1,:),optimal_cloud_interpolated{1}(2,:),'+')
%                 plot(z_yj(1,:),z_yj(2,:),'ro')
%                 drawnow
%                 hold off
                
                
                
                
               
        end
        z_yj_transported =optimal_cloud_interpolated{end}; 
        final_cloud{k} = optimal_cloud_interpolated{end};

% scatter(optimal_cloud_interpolated{end}(1,:),optimal_cloud_interpolated{end}(2,:),30,color,'+');
% 
% hold on
% plot(z_xi(1,:),z_xi(2,:),'b+')
% plot(z_yj(1,:),z_yj(2,:),'r*')
% drawnow
% hold off
    %         size(optimal_cloud_interpolated{end})
% cov(optimal_cloud_interpolated{end}')
% mean(optimal_cloud_interpolated{end}')
end
%             


    
    



% figure('Name','Final Transported Data')
%  color =[0.4660 0.6740 0.1880];
% % 

% cov(cloud_interpolated{3})
% mean(cloud_interpolated{3})
% plot(cloud_interpolated{}(:,1),cloud{3}(:,2),'b+');
    