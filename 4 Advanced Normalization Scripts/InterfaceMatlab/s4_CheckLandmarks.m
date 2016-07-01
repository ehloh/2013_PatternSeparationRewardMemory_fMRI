% s3_CheckregInverseROIs
clear all;close all hidden; clc

where.where='F:\1 fMRI analysis'; where.where='D:\1 [Context-Memory] fMRI Data\1 fMRI analysis';
where.exp_folder='D:\1 [Context-Memory] fMRI Data'; where.data_beh=[where.where '\2 Behavioural data'];where.data_brain=[where.exp_folder filesep '1 MRI data'];
where.antsfolder='D:\host';

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


%%


