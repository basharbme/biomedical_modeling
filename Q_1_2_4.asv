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

            %MCMC PROCEDURE TO GENERATE SAMPLES
            Avox = dwis(:,52,62,25);

            no_of_samples=10; %you can change the number of samples
            stats=zeros(33,no_of_samples);
            
            for k=1:33
                stats(k,:) = Avox(k,1)+0.1*Avox(k,1)*random('norm',1,1,no_of_samples,1);
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
            
            theta_map(i,j)=std(parameters(:,4));
            phi_map(i,j)=std(parameters(:,5));
        end
    end
end


%%
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

disp('95 percent range for phi:')
mean(parameters(:,4))-1.96*std(parameters(:,4))
mean(parameters(:,4))+1.96*std(parameters(:,4))

disp('95 percent range for theta:')
mean(parameters(:,5))-1.96*std(parameters(:,5))
mean(parameters(:,5))+1.96*std(parameters(:,5))

%%
angle_uncertainity=zeros(112,112);
for i=1:112
    for j=1:112
        angle_uncertainity(i,j)=atan2(theta_map(i,j),phi_map(i,j))*180/pi;
    end
end

imshow(angle_uncertainity,[]);


