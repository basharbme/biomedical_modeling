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

%BOOTSTRAP PROCEDURE TO GENERATE SAMPLES
Mvox = dwis(:,42-2:42+2,62,25); %five samples were used to resample to create 1000 different data sets

no_of_samples=1000; %you can change the number of samples
stats=zeros(33,no_of_samples);
for k=1:33
    stats(k,:) = bootstrp(no_of_samples,@mean,Mvox(k,:));
end

parameters=zeros(no_of_samples,5);
for k=1:no_of_samples
% Define a starting point for the non-linear fit
startx = [250000 1E-3 0.5 0 0]; %



% Define various options for the non-linear fitting
% algorithm.
% Note: matlab version 2011b replaces the second line
% below with 'Algorithm', 'levenberg-marquardt',...
h=optimset('MaxFunEvals',20000,...
   'LevenbergMarquardt','on',...
   'TolX',1e-10,...
   'TolFun',1e-10);

%Define paramters for constraint nonlinear optimization
A=[];
b=[];
Aeq=[];
beq=[];
lb=[0 0 0 -3.14 -3.14];
ub=[0.25*startx(1) 0.25*startx(2) 1 3.14 3.14] 
nonlcon = [];

% Now run the fitting
%constrained nonlinear optimzation
            [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fmincon('BallStickSSD',startx,A,b,Aeq,beq,lb,ub,nonlcon,h,stats(:,k),bvals,qhat);
parameters(k,:)=parameter_hat;
end


%%
%DISPLAYING 2-SIGMA RANGE
disp('2-sigma range for S0:')
mean(parameters(:,1))-2*std(parameters(:,1))
mean(parameters(:,1))+2*std(parameters(:,1))

disp('2-sigma range for d:')
mean(parameters(:,2))-2*std(parameters(:,2))
mean(parameters(:,2))+2*std(parameters(:,2))

disp('2-sigma range for f:')
mean(parameters(:,3))-2*std(parameters(:,3))
mean(parameters(:,3))+2*std(parameters(:,3))

%DISPLAYING 95% RANGE
disp('2-sigma range for S0:')
mean(parameters(:,1))-1.96*std(parameters(:,1))
mean(parameters(:,1))+1.96*std(parameters(:,1))

disp('95 percent range for d:')
mean(parameters(:,2))-1.96*std(parameters(:,2))
mean(parameters(:,2))+1.96*std(parameters(:,2))

disp('95 percent range for f:')
mean(parameters(:,3))-1.96*std(parameters(:,3))
mean(parameters(:,3))+1.96*std(parameters(:,3))



