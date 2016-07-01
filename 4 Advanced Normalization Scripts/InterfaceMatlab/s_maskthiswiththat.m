%  maskthiswiththat
% clear all, close all hidden, clc

cd('/Users/EleanorL/Desktop/1 CONTEXT fmri data/2 Second level results s4FullCardiacWithDerivAnts')
% cd('D:\1 [Context-Memory] fMRI Data\2 Second level results s4FullCardiacWithDerivAnts')


% FL ='m_ci13_ContextStickItemNomotor_Hit_Landmarks5'; 
FL ='m_ci3_ContextItemNomotor_Hit_Landmarks5';
% FL ='m_ci2_ContextonlyItem_Hit_Landmarks5';
% FL ='m_i3_ItemContextpresonly_Hit_Landmarks5'; 
% FL ='m_i4_ItemContextevent_Hit_Landmarks5';
% FL ='m_i6_ItemMempmodContextpresonly_Hit_Landmarks5'; 
% FL ='m_c4_ContextallItempresonly_Hit_Landmarks5'; 
% wherethis=[FL filesep 'cm_m1_2x2' filesep]; 
wherethis=[ FL filesep 'co_m1_2x2' filesep]; 
% wherethis=[ FL filesep 'iom_m1_2x2x2' filesep]; 
% wherethis=[ FL filesep 'im_m1_2x2' filesep]; 

this={ 
    'OrthogMask05' 
};

% wherethat='D:\1 [Context-Memory] fMRI Data\3a Anatomical Ants\1 A priori ROI\Resliced s4WDFC Ants\';    % FIXED anatomical masks 
wherethat= '/Users/EleanorL/Desktop/1 CONTEXT fmri data/3a Anatomical Ants/1 A priori ROI/Resliced s4Ants/';    % FIXED anatomical masks 
that={ 
    'rHPC2_bilat'        '_HPC2'
    'rSNVTA3_bilat'      '_SNVTA3'
};




%% Create Orthog Masks x Anatomy!
putwhere= wherethis; 

for i= 1:size(that,1)
    % BATCH
    matlabbatch{1}.spm.util.imcalc.input = {
        [wherethis this{1} '.img,1']
        [wherethat that{i,1} '.nii,1']   };
    matlabbatch{1}.spm.util.imcalc.output = [this{1}  that{i,2}  '.img'];
    matlabbatch{1}.spm.util.imcalc.outdir = {putwhere};
    matlabbatch{1}.spm.util.imcalc.expression = 'i1.*(i2>0)';
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];
    
end
disp('DONE ! :)') 
