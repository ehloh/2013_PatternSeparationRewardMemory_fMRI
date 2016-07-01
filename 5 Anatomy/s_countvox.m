% Count voxels of of overlap (and other things) 
clear all, close all hidden, clc


where='D:\1 [Context-Memory] fMRI Data\3a Anatomical Ants\1 A priori ROI\';
rois={
    [where 'HPC2_Left.nii'];
    [where 'HPC2_Right.nii'];
    [where 'HPC2_anterior_Left.nii'];
    [where 'HPC2_anterior_Right.nii'];
    [where 'SNVTA3_Left.nii'];
    [where 'SNVTA3_Right.nii'];
}; 


%% 

rvol= cell(size(rois,1),1); 
for r=1:size(rois,1)
    rvol{r} = spm_read_vols(spm_vol(rois{r} ));
    rvol{r,2} =sum( rvol{r}(:)~=0); 
end

rvol


% rvol

% sum( (rvol{1}(:)~=0))
% sum( (rvol{2}(:)~=0))

% disp(['No. overlapping voxels = '  num2str(sum( (rvol{1}(:)~=0).* (rvol{2}(:)~=0)))])


