%%
% mainGenerateVideo
% This script loads a segmented 3D-image of a brain, then places a
% synthetic tumour and lets it grow according to the PDE-model
% Two video files are generated that visualize 1) the tumour placement in the
% brain and 2) the shape of the tumour
%
% This script solves the mathematical model with p=2, rho=const, a=b=c=1,
% D(x)=d(x)*eye. Synthetic tumour is planted in.
%%
clear 
%Gunzip files
files = gunzip('Dataset\per_subject_MRI_volumes\MTP_2023_0004','Unziped_Dataset\MTP_2023_0004');
files_seg = gunzip('BraTS_version\BraTS2021_01516_seg.nii.gz','BraTS_version\');
segvol1 = niftiread("T1_half_seg.nii");
vol1 = niftiread('Dataset_unzip\MTP_2023_0004\T1.nii');
segvol2 = niftiread('T2_half_seg.nii');
seg_tumor = niftiread('BraTS_version\BraTS2021_01516_seg.nii');
%%


% Parameters used in the model
% High-Grade Gliomas (HGG)
rho = 0.012; % Proliferation rate
dw  = 0.65;  % Diffusion coefficient, white matter
dg  = 0.13;  % Diffusion coefficient, grey matter
alpha = 1;
beta = 1;
gamma = 1;
% PDE-solver parameters
tStep  = 0.05;    % Discretization step for variable t
noIter = 50;    % Number of iterations in t-variable
freqImgSave = 10; % Every 10th image is to be saved (adjust FrameRate for myVideo1 and myVideo2 below accordingly)

% Load brain-data
%load('coreg_7940_mni152_2009bet.mat', 'segvol', 'grayMatterSegmentationValues', 'whiteMatterSegmentationValues');
% segvol contains a segmented 3D-scan of a brain (256x256x256 voxels)
% ...MatterSegmentationValues contain a list of values in segvol that
%                    correspond to grey and whitematter, respectively

% Generate the D-function that appears in the PDE
D =   dg * ismember(segvol1,2) + dw * ismember(segvol1,3);
% Crop D so that it is not surrounded by unnecessary/dummy zero voxels
idx1 = squeeze(any(any(D, 2), 3));
idx2 = squeeze(any(any(D, 1), 3));
idx3 = squeeze(any(any(D, 1), 2));
D = D(idx1, idx2, idx3);

% D however needs to be surrounded by at least one dummy layer of zeroes
% due to boundary conditions of the model
D = padarray(D, [1, 1, 1], 0, 'both');

% Generate Omega, the domain in which the PDE is being solved
Om = +(D > 0);
OmSize = size(Om);

% Plant a seed of a tumour
if false  % Change to true if the seed is to be placed randomly. If false place the tumour "in the middle"
    idxList = find(Om);
    seedDiam   = 19;
    [seedPos(1), seedPos(2), seedPos(3)] = ind2sub(OmSize, idxList(randi([1, length(idxList)])));
else    
    % Position and diameter of the seed
    seedPos    = [60,100,100];
    seedDiam   = 4;
    % make sure that the seed is planted in Om
    seedTmp = find(D(seedPos(1), seedPos(2), :));
    seedTmp(seedTmp < seedPos(3)) = [];
    seedPos(3) = seedTmp(3); 
end

% Generate the seed
u1tmp = reshape(((1:OmSize(1)) - seedPos(1)).^2, [OmSize(1), 1, 1]);
u2tmp = reshape(((1:OmSize(2)) - seedPos(2)).^2, [1, OmSize(2), 1]);
u3tmp = reshape(((1:OmSize(3)) - seedPos(3)).^2, [1, 1, OmSize(3)]);
u = min(1, 10*exp(-(u1tmp+u2tmp+u3tmp)*3/seedDiam));
u(u<1e-3) = 0;
u = u .* Om;
phi = u;

% Show a 3D-image of the brain and the planted seed
fig1 = figure(1);
%fig1.Position = [10, 200, 1280, 720];

myVideo1 = VideoWriter('tumGrowth.mp4', 'MPEG-4');
myVideo1.FrameRate = 10;
myVideo1.Quality = 98;
open(myVideo1);

tumorPatch = visTum3d(D, u, seedPos);
drawnow;

frame1 = getframe(fig1);
writeVideo(myVideo1, frame1);

% Show a 3D-image of the tumour only
fig2 = figure(2);
%fig2.Position = [10, 200, 1280, 720];

myVideo2 = VideoWriter('tumOnlyGrowth.mp4', 'MPEG-4');
myVideo2.FrameRate = 10;
myVideo2.Quality = 98;
open(myVideo2);

tumorOnlyPatch = visTumOnly3d(u);
drawnow;
frame2 = getframe(fig2);
writeVideo(myVideo2, frame2);

% Prepare auxilliary variables for Finite Difference Method
% Edges of Omega (needed for the Neumann condition)
OmEdges.Pos1 = 1 - circshift(1-Om, -1, 1) .* Om;
OmEdges.Neg1 = 1 - circshift(1-Om,  1, 1) .* Om;
OmEdges.Pos2 = 1 - circshift(1-Om, -1, 2) .* Om;
OmEdges.Neg2 = 1 - circshift(1-Om,  1, 2) .* Om;
OmEdges.Pos3 = 1 - circshift(1-Om, -1, 3) .* Om;
OmEdges.Neg3 = 1 - circshift(1-Om,  1, 3) .* Om;

% Averages of D in each direction (needed in the elliptic diff. operator)
% Davg.Pos1 = 0.5 * (D(x, y+1, z) + D(x, y, z))
% Davg.Neg1 = 0.5 * (D(x, y-1, z) + D(x, y, z))
%      and analogously for Pos2 (x), Neg2 (x), Pos3 (z), and Neg3 (z)
Davg.Pos1 = (circshift(D, -1, 1) + D)/2;
Davg.Neg1 = (circshift(D,  1, 1) + D)/2;
Davg.Pos2 = (circshift(D, -1, 2) + D)/2;
Davg.Neg2 = (circshift(D,  1, 2) + D)/2;
Davg.Pos3 = (circshift(D, -1, 3) + D)/2;
Davg.Neg3 = (circshift(D,  1, 3) + D)/2;

% run the model, i.e., solve the PDE
for currIter=1:noIter
  fprintf('Iteration %d\n', currIter);
  u = u + tStep * (divDdu(u, Davg, OmEdges) + rho * (u.^alpha)*beta .* ...
      (1-u.^(1/beta)).^gamma);
   if mod(currIter, freqImgSave) == 0
    figure(fig1);
      tumorPatch = visTum3dUpdate(u, tumorPatch);
      drawnow;
      frame1 = getframe(fig1);
      writeVideo(myVideo1, frame1);

     figure(fig2);
      tumorOnlyPatch = visTum3dUpdate(u, tumorOnlyPatch);
      drawnow;
      frame2 = getframe(fig2);
      writeVideo(myVideo2, frame2);
  end
end

writeVideo(myVideo1, frame1);
close(myVideo1);

writeVideo(myVideo2, frame2);
close(myVideo2);

save('reactDiffModelSolution.mat', 'phi', 'u', 'seedPos', 'seedDiam', ...
    'D', 'Om', 'tStep', 'noIter');

%%
%Add tumor in Omega
Om_wt = ceil(Om + u);


%%
% Find the volume of synthetic tumor
% 1 voxel = 1mm^3
nonzero_indices = find(u>0.05); % find indices of non-zero elements
[num_nonzero, ~] = size(nonzero_indices); % count non-zero elements

% convert linear indices to subscripts
[sub1, sub2, sub3] = ind2sub(size(u),nonzero_indices);

% print the subscripts of non-zero elements
disp("Subscripts of non-zero elements in synthetic tumor:");
disp([sub1, sub2, sub3]);

% print the number of non-zero elements
disp("Number of non-zero elements in synthetic tumor:");
disp(num_nonzero);

%%
% Find the volume of tumor
% 1 voxel = 1mm^3
nonzero_indices2 = find(seg_tumor); % find indices of non-zero elements
[num_nonzero2, ~] = size(nonzero_indices2); % count non-zero elements

% convert linear indices to subscripts
[sub12, sub22, sub32] = ind2sub(size(seg_tumor),nonzero_indices2);

% print the subscripts of non-zero elements
disp("Subscripts of non-zero elements in synthetic tumor:");
disp([sub12, sub22, sub32]);

% print the number of non-zero elements
disp("Number of non-zero elements in synthetic tumor:");
disp(num_nonzero)
disp(num_nonzero2);
