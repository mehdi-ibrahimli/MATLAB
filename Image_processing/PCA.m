%Mehdi Ibrahimli
%ID 2040467
clear;
clc;
%read the hdf
data=hdfread('myfile.L1R','EO1H1660322004097110PX.L1R');

%256x256 area
clipdata=data(1000:1255,:,:);
clipdata=permute(clipdata,[1,3,2]);

%% plotting
plot(clipdata(:,100,40),clipdata(:,100,39),'b.');
hold on
title 'Band 39 vs Band 40'
xlabel 'Band 39'
ylabel 'Band 40'
grid on
hold off

%% rescaling
B_39 = reshape(clipdata(:,:,39),[],1);
B_39 = int16(rescale(B_39,0,100));

B_40 = reshape(clipdata(:,:,40),[],1);
B_40 = int16(rescale(B_40,0,100));
all=[B_39,B_40];
%% histogram
hist = zeros(101,101);
for i = 1:101
    for j=1:101
        hist(i,j)=sum(all(:, 1) == i & all(:, 2) == j);
    end
end
%%
mesh(hist);
xlabel('Band 39');
ylabel('Band 40');
%% contour plot
contour(hist);
xlabel('Band 39');
ylabel('Band 40');
%% covariance

C = double(reshape(clipdata,[],242,1));
Substracted_means = bsxfun(@minus, C, mean(C));
COvar=(Substracted_means.'*Substracted_means)/65535;
% correlation
sigma = sqrt(diag(COvar));
COR = bsxfun(@rdivide,COvar,sigma);
COR = bsxfun(@rdivide,COR,sigma');
imagesc(COR)
colormap jet

%% PCA
[V,D] = eig(COvar);
[d,ind] = sort(diag(D),'descend');
Ds = D(ind,ind);
Trmatrix = V(:,ind);
PCA=Trmatrix'*(C');
PCA = reshape(PCA',256,256,242);
PC1=PCA(:,:,1);
PC2=PCA(:,:,2);
%% image 1
imagesc(PC1);
%% image 2
imagesc(PC2);
%% scatter plot

plot (PC1(:,100),PC2(:,100),'b.');

