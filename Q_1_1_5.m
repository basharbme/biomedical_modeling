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


S0=dwis(1,:,:,:); %baseline image
S=dwis([2:length(bvals)],:,:,:); %other diffusion images

GradientOrientations=qhat';
GradientOrientations=GradientOrientations([2:length(bvals)],:); % the first orientation is neglected since b=0

b_value = bvals';
b_value = b_value([2:length(b_value)]); %the b values where b~=0

G=constructMatrixOfMonomials(GradientOrientations, 2);
C=constructSetOf81Polynomials(2)';

P=G*C;
%%
P=-diag(b_value)*P;

%%
DTI=zeros(3,3,size(S,2),size(S,3));

UniqueTensorCoefficients = zeros(6,size(S,2),size(S,3));
for i=1:size(S,2)
   for j=1:size(S,3)
      
      if (S0(:,i,j,25)~=0) %execute the procedure only for those voxels that are defined as 1 in the binary image to exclude outside object voxels
        
        arr = S(:,i,j,25);
        zero_detect=find(arr==0);
        
        if (isempty(zero_detect)==1)
        
            y=log(squeeze(S(:,i,j,25)/S0(:,i,j,25))); %log (S/S0) for each voxel
      
            x=lsqnonneg(P, y); %solving nonlinear least square problem
            
            T = C * x;
            
            UniqueTensorCoefficients(:,i,j)=T; %tensor coefficient
            
            DTI(:,:,i,j)=[T(6) T(5)/2 T(4)/2
            T(5)/2 T(3) T(2)/2
            T(4)/2 T(2)/2 T(1)]; %estimation of diffusion tensor
        
        end
    
      end
   end
end

%%
figure;
for i=1:6
    subplot(3,2,i);imshow(squeeze(UniqueTensorCoefficients(i,:,:)),[]); title(strcat('Tensor Coefficient-',num2str(i)));
end