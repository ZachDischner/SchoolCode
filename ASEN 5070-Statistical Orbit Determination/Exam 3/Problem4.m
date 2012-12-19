%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   ASEN 5070-Problem 4
% 
% Zach Dischner
%   Exam 3
%       Problem 4
% 
% 
% 
% Solves and answers questions relating to problem 4 of the STATOD final
% exam
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prepare Workspace
clc;clear all;close all

%% Load Data
tmp=load('Exam3_Problem4_data.txt');

%% Form into useful datasets
for ii=1:8
    data(ii).xy = tmp(ii,2:3);
    x(ii) = tmp(ii,2);
    y(ii) = tmp(ii,3);
    data(ii).P  = reshape(tmp(ii,4:end),2,2);
    data(ii).sigma = data(ii).P(1,1)^.5;
end


%% a - Find best guess of xhat, and covariance matrix
sumP_inv=0;
sumx=0;
xhat = [0 0]';
% Perform Summing Algorithm
for ii=1:8
    sumP_inv = sumP_inv + inv(data(ii).P);
    sumx = sumx + inv(data(ii).P)*data(ii).xy';
end

xhat = inv(sumP_inv)*sumx;
Pfinal = inv(sumP_inv);

fprintf('\na:\n\n')
fprintf('The guess for xhat is: [%3.3f,%3.3f]',xhat(1),xhat(2))
fprintf('\n\nWith a covariance matrix of:  [%3.3f,%3.3f]',Pfinal(1),Pfinal(2))
fprintf('\n..............................[%3.3f,%3.3f]\n',Pfinal(3),Pfinal(4))
  


%% b - Plot points and ellipses
figure
plot(x,y,'r*','MarkerSize',12)
hold on
plot(xhat(1),xhat(2),'xk','MarkerSize',15)

% Plot Error Ellipses
cov_ellipse(Pfinal,xhat(1),xhat(2),'G');
cov_ellipse(Pfinal.*2^2,xhat(1),xhat(2),'m');
cov_ellipse(Pfinal.*3^2,xhat(1),xhat(2),'r');

legend('$\hat{x}_i$','$\hat{x}_{comb}$','$1-\sigma$','$2-\sigma$','$3-\sigma$')
xlabel('x');ylabel('y');title('Problem 4b - Combined $\hat{x}$ Estimate')



%% C - Distance of xhat6
% Move point of interest to origin
xhat2 = xhat - xhat;
xhat6 =[x(6),y(6)]' - xhat;
figure
plot(xhat6(1),xhat6(2),'r*','MarkerSize',12)
hold on
plot(xhat2(1),xhat2(2),'xk','MarkerSize',15)

cov_ellipse(Pfinal,xhat2(1),xhat2(2),'G');
cov_ellipse(Pfinal*2^2,xhat2(1),xhat2(2),'m');
cov_ellipse(Pfinal*3^2,xhat2(1),xhat2(2),'r');

a= Pfinal(1,1)^.5;
b=Pfinal(2,2)^.5;


% Get U, the rotation vector
[vec,val]=eigs(Pfinal);
U=vec;
%Rotate point into principle axis
xhat6prime = U'*xhat6;

plot(xhat6prime(1) ,xhat6prime(2),'b*','Markersize',12)

% Scale
a= Pfinal(1,1)^.5;
b=Pfinal(2,2)^.5;

% Angle of PA rel to xy
phi = atan2(U(2,1),U(1,1));

% Angle of xhat6 rel to XY
th=atan2(xhat6prime(2),xhat6prime(1));
% Angle of xhat6 rel to PA
rel=th-phi;
% Now get distance
d=norm(xhat6prime);
delx=d*cos(rel);
dely=d*sin(rel);
delxstd=delx/a;
delystd=dely/b;

delstd = sqrt(delxstd^2 + delystd^2);

cov_ellipse(Pfinal.*delstd^2,xhat2(1),xhat2(2),'y')

legend('$\hat{x}_6$','$\hat{x}_{comb}$','$1-\sigma$','$2-\sigma$','$3-\sigma$','$\hat{x}_6''$ Rel to $\hat{x}_{comb}$','15.708-$\sigma$','location','best')
xlabel('x');ylabel('y');title('Problem 4c - Distance from $\hat{x}_{comb}$ to $\hat{x}_6$')

fprintf('\nc:\n\n')
fprintf('\n\nxhat6 is %3.3f STD_comb from xhat in X',delxstd);
fprintf('\nxhat6 is %3.3f STD_comb from xhat in Y',delystd);
fprintf('\nFor a total sigma_comb distance of %3.3f\n\n',delstd)



%% d - Find dist from x6 to xhat
% Move point of interest to Origin
xhat6 =[x(6),y(6)]' - [x(6),y(6)]';
xhat2 = xhat - [x(6);y(6)];

P6 = data(6).P;
figure
plot(xhat6(1),xhat6(2),'r*','MarkerSize',12)
hold on
plot(xhat2(1),xhat2(2),'xk','MarkerSize',15)

cov_ellipse(P6,xhat6(1),xhat6(2),'g')

% STD distances
a= P6(1,1)^.5;
b= P6(2,2)^.5;

% Get rotation vector
[vec,val]=eigs(P6);
U=vec;

% Angle of PA rel to XY
phi = atan2(U(2,1),U(1,1));

% Rotate point about orgin
xhat2prime = U'*(xhat2);
% Angle of xhat rel to XY
th=atan2(xhat2prime(2),xhat2prime(1));
% Angle of xhat rel to PA
rel=th-phi;
% Find distances
d=norm(xhat2prime);
delx=d*cos(rel);
dely=d*sin(rel);
delxstd=delx/a;
delystd=dely/b;

delstd = sqrt(delxstd^2 + delystd^2);

cov_ellipse(P6*delstd^2,xhat6(1),xhat6(2),'y');

fprintf('\nd:\n\n')
fprintf('\n\nxhat is %3.3f STD_6 from xhat6 in X',delxstd);
fprintf('\nxhat6 is %3.3f STD_6 from xhat6 in Y',delystd);
fprintf('\nFor a total sigma_6 distance of %3.3f\n\n',delstd)

legend('$\hat{x}_6$','$\hat{x}_{comb}$','$1-\sigma$','1.336-$\sigma$','location','best')
xlabel('x');ylabel('y');title('Problem 4d - Distance from $\hat{x}_{6}$ to $\hat{x}_{comb}$')

%% e - The plot for all covariance matrices
figure
hold on
cmap = hsv(9);
for ii = 1:8
    plot(x(ii),y(ii),'.','Color',cmap(ii,:))
    cov_ellipse(data(ii).P,x(ii),y(ii),cmap(ii,:));
end

plot(xhat(1),xhat(2),'xk','MarkerSize',15)
cov_ellipse(Pfinal,xhat(1),xhat(2),cmap(9,:));
lines = sort(findobj(gca,'Type','line'));
set(lines(end),'linewidth',3)

title('Problem 4e - All Probability Ellipsoids')

legend('$\hat{x}_1$','$P_1$','$\hat{x}_2$','$P_2$','$\hat{x}_3$','$P_3$',...
        '$\hat{x}_4$','$P_4$','$\hat{x}_5$','$P_5$','$\hat{x}_6$','$P_6$','$\hat{x}_7$',...
        '$P_7$','$\hat{x}_8$','$P_8$','$\hat{x}_{comb}$','$P_{\hat{x}_{comb}}$','location','best')
    
xlabel('x');ylabel('y');title('Problem 4e - All Covariance Ellipsoids')    
    

%% g - Monte Carlo
figure
hold on
plot(xhat(1),xhat(2),'xk','MarkerSize',15)
cov_ellipse(Pfinal,xhat(1),xhat(2),'b');
cov_ellipse(Pfinal*2^2,xhat(1),xhat(2),'y');
cov_ellipse(Pfinal*3^2,xhat(1),xhat(2),'r');


S = (chol(Pfinal))';
x_montbar = [0;0];
P_mont = [0 0;0 0];
for ii = 1:1000
    e=randn(2,1);
    x_mont(:,ii)=(S'*e)+xhat;
    x_montbar = (x_montbar + x_mont(:,ii))/2;
    P_mont = (P_mont + (x_mont(:,ii) -x_montbar)*(x_mont(:,ii)-x_montbar)')/2; 
end
plot(x_mont(1,:),x_mont(2,:),'k.','markersize',1)

plot(x_montbar(1),x_montbar(2),'og','Markersize',15)


cov_ellipse(P_mont,x_montbar(1),x_montbar(2),'G');


legend('$\hat{x}_{comb}$','1-$\sigma$','2-$\sigma$','3-$\sigma$','Monte Carlo Iter','$\bar{x}$','$P_{MonteCarlo}$','location','best')
xlabel('x');ylabel('y');title('Problem 4g - Monte Carlo Simulation')   
    
fprintf('\nh:\n\n')
fprintf('The guess for xhat is: [%3.3f,%3.3f]',x_montbar(1),x_montbar(2))
fprintf('\n\nWith a covariance matrix of:  [%3.3f,%3.3f]',P_mont(1),P_mont(2))
fprintf('\n..............................[%3.3f,%3.3f]\n',P_mont(3),P_mont(4))
  

