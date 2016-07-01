% s3_CheckregInverseROIs
clear all;close all hidden; clc

where.where='F:\1 fMRI analysis'; where.where='D:\1 [Context-Memory] fMRI Data\1 fMRI analysis';
where.exp_folder='D:\1 [Context-Memory] fMRI Data'; where.data_beh=[where.where '\2 Behavioural data'];where.data_brain=[where.exp_folder filesep '1 MRI data'];
where.antsfolder='D:\host';


%
request.checkreg=1;
request.display=0;
%
request.roiname='_LM5cm_SRvSN_001_HPC_aLpeak';
request.AntsType='Landmarks5';
%
request.checkreg_n_subsper=1;

% Which subjects
log.specificsubjects={}; 
request.WhichFLmodel='m_ci3_ContextItemNomotor_Hit';

for o1=1:1 % General settings and specifications
    
    % Load subjects
    addpath(where.where)
    log.w=load([where.data_brain filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    log.allsubjects=log.subjects; 
    [w.s w.s1 log.koshertable]=xlsread([where.where filesep 'i8_subjectdataok_secondlevel.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable, request.WhichFLmodel);
    
    
end


%% Check reg

if request.checkreg
% Assemble files
k=1; allfiles=cell(log.n_subjs*2,1);
for s=1:log.n_subjs
    allfiles{k}=[where.antsfolder '\2_Data_AdjustCons\' log.subjects{s} filesep request.AntsType '\InverseTransform\' log.subjects{s} request.roiname '.nii,1'];
    allfiles{k+1}=[where.antsfolder '\2_Data_AdjustCons\' log.subjects{s} filesep log.subjects{s} '_T1w_coreg.nii,1'];
    k=k+2;
end

% Put into figures
n_figs=size(allfiles,1)/(request.checkreg_n_subsper*2); files=cell(n_figs,1);
for f=1:n_figs
    files{f}=allfiles((f-1)*(request.checkreg_n_subsper*2)+1:f*(request.checkreg_n_subsper*2));
end

% Display
for f=1:n_figs
    matlabbatch{1}.spm.util.checkreg.data=files{f};
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    input('Next batch?  ');
    matlabbatch=[];
end


end

%% Display (single subject)


