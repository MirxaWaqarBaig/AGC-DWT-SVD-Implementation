%% 
% un- official implementation of 
% Kallel, Fathi, and Ahmed Ben Hamida. "A new adaptive gamma correction 
% based algorithm using DWT-SVD for non-contrast CT image enhancement.
% " IEEE transactions on nanobioscience 16.8 (2017): 666-675.

% this (main) file is used to read dicom file and apply all the required functions. 
% you can change the path or file according to the need  
clear all 
clc
%%
img = dicomread('Subject_13.dcm');
img = im2double(img);
img(isinf(img)|isnan(img)) = 0;
%%
%img = mat2gray(img);
%% converting image into uint8
into=linspace(min(img(:)),max(img(:)),256)
into=uint8(arrayfun(@(x) find(abs(into(:)-x)==min(abs(into(:)-x))),img))
%% First part (GHE on the original image) Left side first, of the proposed technique
%% GHE
GH_eq = GHE_gray(into);
%%
% imshow(GH_eq,[]);  %Result of GHE on the original image
%% DWT on this equalized image
i= GH_eq;
sX=size(i);
%Splits image into 4 layers according to paper
[LL_GH,LH_GH,HL_GH,HH_GH]=dwt2(i,'db1');
figure(1)
subplot(2,2,1);imshow(LL_GH,[]);title('LL band of image');
subplot(2,2,2);imshow(LH_GH,[]);title('LH band of image');
subplot(2,2,3);imshow(HL_GH,[]);title('HL band of image');
subplot(2,2,4);imshow(HH_GH,[]);title('HH band of image');
%% SVD on the LL component  of enhanced 
[u_gh,s_gh,v_gh] = svd(LL_GH);
%% Computing Maximum of U_gh and V_gh 
max_u_gh = max(u_gh(:));
max_v_gh = max(v_gh(:));
%% Right hand side approach of the technique
%% DWT on the original image 
i_org= img;
sX=size(i_org);
%Splits image into 4 layers
[LL,LH,HL,HH]=dwt2(i_org,'db1');
figure(1)
subplot(2,2,1);imshow(LL,[]);title('LL band of original image');
subplot(2,2,2);imshow(LH,[]);title('LH band of original image');
subplot(2,2,3);imshow(HL,[]);title('HL band of original image');
subplot(2,2,4);imshow(HH,[]);title('HH band of original image');
%% SVD on the LL component original
[u_org,s_org,v_org] = svd(LL);
%% Computing Maximum of U_org and V_org
max_u_org = max(u_org(:));
max_v_org = max(v_org(:));
%% Now merging left and right side top down appraoch 
%% Correction Factor 
cor_f =  (max_u_gh + max_v_gh / max_u_org + max_v_org) ;
%% Enhanced Singular Value
en_sv = cor_f*s_gh;
%% LL SVD (inverse SVD) 
ll_isvd = u_gh*s_gh*v_gh';

%% Mean and Standard Deviation of (LL_SVD) image
mean_isvd = mean2(ll_isvd);
std_isvd = std2(ll_isvd);
%% Difference Formula
D = (mean_isvd + 2*std_isvd)-(mean_isvd - 2*std_isvd);  % if it is required
%% Checking the class of the image 
%% First Formula 
taw = 3;
f_std = 4*std_isvd;
%% Classifying usin If statement 
if f_std <= taw/1;
    display('Low contrast image');
else 
    display('Moderate level contrast image');
end
%% now as it is a moderate level contrast image 
%%
gamma = exp(1 - (mean_isvd + std_isvd) ./ 2);
%%
[rows columns] = size(ll_isvd);
%%
onesMat = ones(rows, columns, 'double');
%%
k_mat = onesMat - ll_isvd;
%%
k_mat = k_mat .* power(mean_isvd, gamma);
%%
k_mat = ll_isvd + k_mat;
%%
denominator = k_mat - onesMat;
%%
denominator = denominator .* heaviside(0.5 - mean_isvd);
%%
denominator = onesMat + denominator;
%%
c_param = onesMat ./ denominator;
%%
SVD_LL_img = c_param .* ll_isvd;
%%

%%
imshow(SVD_LL_img,[]);
%%
ints=linspace(min(SVD_LL_img(:)),max(SVD_LL_img(:)),256)
ints=uint8(arrayfun(@(x) find(abs(ints(:)-x)==min(abs(ints(:)-x))),SVD_LL_img))
%%
% imshow(ints,[]);


%% Inverse WT

I_wt = idwt2(SVD_LL_img,LH,HL,HH,'db1');
%%
% I_wt = idwt2(LL_GH,LH,HL,HH,'db1');
% %%
% I_wt_2 = mat2gray(I_wt);
% I_wt_2 = I_wt_2 / 255;
%%
imshow(I_wt,[]);
% imshow(I_wt_2,[]);