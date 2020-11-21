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

Avox = dwis(:,42,62,25);

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
 
% Now run the fitting
[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD',startx,h,Avox,bvals,qhat);
parameter_hat(1)
parameter_hat(2)
parameter_hat(3)
parameter_hat(4)
parameter_hat(5)

[sumRes] = BallStickSSD1(parameter_hat, Avox, bvals, qhat)
RESNORM

%function minimization with the constraints
[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSDC',startx,h,Avox,bvals,qhat);
parameter_hat(1)
parameter_hat(2)
parameter_hat(3)
parameter_hat(4)
parameter_hat(5)
RESNORM
[sumRes] = BallStickSSD1(parameter_hat, Avox, bvals, qhat)

