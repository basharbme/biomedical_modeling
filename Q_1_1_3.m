clc; clear all;

fid=fopen('dwi.Bfloat', 'r', 'b');
dwis = fread(fid, 'float');
fclose(fid);

% B = reshape(A,m,n,p,...) or B = reshape(A,[m n p ...]) 
% returns an n-dimensional array with the same elements as A but reshaped to have the size m-by-n-by-p-by-....
% The product of the specified dimensions, m*n*p*..., must be the same as prod(size(A)).
prod(size(dwis));
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
    'LargeScale','on',...
   'LevenbergMarquardt','on',...
   'TolX',1e-10,...
   'TolFun',1e-10);

re = 100; %repeating with 100 different starting points (you can change it if you want to)

val = zeros(1,re); %for RESNORM for all realization

for i=1:re
    startx_diff(1) = random('normal',startx(1),1);
    startx_diff(2) = startx(2);
    startx_diff(3) = startx(3);
    startx_diff(4) = 0;
    startx_diff(5) = 0;
    
    %function minimization with the constraints
    [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSDC',startx_diff,h,Avox,bvals,qhat);
    val(1,i)=RESNORM;
end

figure;plot(val,'o-'); title('Residuals');