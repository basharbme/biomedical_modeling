clc; clear all;

fid=fopen('dwi.Bfloat', 'r', 'b');
dwis = fread(fid, 'float');
fclose(fid);

dwis = reshape(dwis, 33, 112, 112, 50);

% Middle slice of the first image volume, which has b=0.
imshow(squeeze(dwis(1,:,:,25)), []);
% Middle slice of the second volume, which has b=1000.
figure;imshow(squeeze(dwis(2,:,:,25)), []);


fid=fopen('grad_dirs.txt','r','b');
qhat=fscanf(fid,'%f',[3,inf]);
fclose(fid);

bvals=1000*sum(qhat.*qhat);

Avox = dwis(:,52,62,25);

% Define a starting point for the non-linear fit
startx = [250000 1E-3 0.5 0 0];



% Define various options for the non-linear fitting
% algorithm.
% Note: matlab version 2011b replaces the second line
% below with 'Algorithm', 'levenberg-marquardt',...
h=optimset('MaxFunEvals',20000,...
   'LevenbergMarquardt','on',...
   'TolX',1e-10,...
   'TolFun',1e-10);
 
A=[];
b=[];
Aeq=[];
beq=[];
lb=[0 0 0 -3.14 -3.14];
ub=[0.25*startx(1) 0.25*startx(2) 1 3.14 3.14] 
nonlcon = [];

re = 100; %repeating with 100 different starting points (you can change it if you want to)

val = zeros(1,re); %for RESNORM for all realization
S0=zeros(1,re);
d=zeros(1,re);
f=zeros(1,re);
theta=zeros(1,re);
phi=zeros(1,re);

final_S0=zeros(1,re);
final_d=zeros(1,re);
final_f=zeros(1,re);
final_theta=zeros(1,re);
final_phi=zeros(1,re);

for i=1:re
    
    %generating random starting points
    startx_diff(1) = startx(1) + startx(1)*random('normal',0,100); 
    startx_diff(2) = startx(2) + startx(2)*random('normal',0,100);
    startx_diff(3) = startx(3) + startx(3)*random('normal',0,100);
    startx_diff(4) = random('unif',-3.14,3.14);
    startx_diff(5) = random('unif',-3.14,3.14);
    
    %function minimization with the constraints
    [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fmincon('BallStickSSD',startx,A,b,Aeq,beq,lb,ub,nonlcon,h,Avox,bvals,qhat);
    
    val(1,i)=RESNORM;
    
    S0(1,i)=startx_diff(1);
    d(1,i)=startx_diff(2);
    f(1,i)=startx_diff(3);
    theta(1,i)=startx_diff(4);
    phi(1,i)=startx_diff(5);
    
    final_S0(1,i)=parameter_hat(1);
    final_d(1,i)=parameter_hat(2);
    final_f(1,i)=parameter_hat(3);
    final_theta(1,i)=parameter_hat(4);
    final_phi(1,i)=parameter_hat(5);
end

figure;
subplot(3,2,1);plot(val,'o-');title('RESNORM for different initial parameters values');
subplot(3,2,2);plot(S0,'o-');title('Different "S0" values initial (blue) final (red)');hold on;
plot(final_S0,'r.-');
subplot(3,2,3);plot(d,'o-');title('Different "d" values initial (blue) final (red)');hold on;
plot(final_d,'r.-');
subplot(3,2,4);plot(f,'o-');title('Different "f" values initial (blue) final (red)');hold on;
plot(final_f,'r.-');
subplot(3,2,5);plot(theta,'o-');title('Different "theta" values initial (blue) final (red)');hold on;
plot(final_theta,'r.-');
subplot(3,2,6);plot(phi,'o-');title('Different "phi" values initial (blue) final (red)');hold on;
plot(final_phi,'r.-');