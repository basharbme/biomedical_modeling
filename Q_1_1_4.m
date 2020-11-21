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

%Define paramters for constraint nonlinear optimization
A=[];
b=[];
Aeq=[];
beq=[];
lb=[0 0 0 -3.14 -3.14];
ub=[0.25*startx(1) 0.25*startx(2) 1 3.14 3.14] 
nonlcon = [];

%Define matrix to store the values of S0, d, f, theta and phi
S0_map=zeros(112,112);
d_map=zeros(112,112);
f_map=zeros(112,112);
theta_map=zeros(112,112);
phi_map=zeros(112,112);


im=squeeze(dwis(1,:,:,25)); %choose a slice form the DWI images
im=mat2gray(im); %convert that slice to gray level
img=zeros(112,112);
l=find(im>0.1); %find the values of the voxel that are greater than 0.1
img(l)=1; %produce a binary images with those voxels

for i=1:112
    i
    for j=1:112
        if (img(i,j)~=0) %execute the procedure only for those voxels that are defined as 1 in the binary image to exclude outside object voxels
            Avox = dwis(:,i,j,25);
            
            %constrained nonlinear optimzation
            [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fmincon('BallStickSSD',startx,A,b,Aeq,beq,lb,ub,nonlcon,h,Avox,bvals,qhat);
            
            %save the values in the matrix
            S0_map(i,j)=parameter_hat(1);
            d_map(i,j)=parameter_hat(2);
            f_map(i,j)=parameter_hat(3);
            theta_map(i,j)=parameter_hat(4);
            phi_map(i,j)=parameter_hat(5);
            
        end
        
        
    end
end

%save all teh parameters maps
save S0_map S0_map;
save d_map d_map;
save f_map f_map;
save theta_map theta_map;
save phi_map phi_map;

%%
[X,Y]=meshgrid(1:112,1:112);
subplot(2,2,1);imshow(S0_map,[]);title('S0 Map for slice 25');
subplot(2,2,2);imshow(d_map,[]);title('d Map for slice 25');
subplot(2,2,3);imshow(f_map,[]);title('f Map for slice 25');
subplot(2,2,4);quiver(X,Y,theta_map.*f_map,phi_map.*f_map);title('Fibre direction map for slice 25'); %fibre direction weighted by volume fraction f