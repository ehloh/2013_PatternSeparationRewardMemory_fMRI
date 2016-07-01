% Re-slice anatomical ROIs to the fxnal data spex
clear all


where='D:\1 [Context-Memory] fMRI Data\3a Anatomical Ants\1 A priori ROI\';

matlabbatch{1}.spm.spatial.coreg.write.ref = {'D:\1 [Context-Memory] fMRI Data\2 Second level results s4FullCardiacWithDerivAnts\m_ci13_ContextStickItemNomotor_Hit_Landmarks5\cm_m1_2x2\OrthogMask05.img,1'};
matlabbatch{1}.spm.spatial.coreg.write.source = {
    [where 'HPC2_Left.nii,1'];  
    [where 'HPC2_Right.nii,1'];
    [where 'HPC2_anterior_Left.nii,1'];  
    [where 'HPC2_anterior_Right.nii,1'];
    [where 'HPC2_anterior_bilat.nii,1'];
    [where 'HPC2_bilat.nii,1'];  
    [where 'HPC2_posterior_Left.nii,1'];
    [where 'HPC2_posterior_Right.nii,1'];
    [where 'HPC2_posterior_bilat.nii,1'];
    [where 'SNVTA3_Left.nii,1'];  
    [where 'SNVTA3_Right.nii,1']; 
    [where 'SNVTA3_bilat.nii,1'];   
                                                 };
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';


spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);







matlabbatch=[];




