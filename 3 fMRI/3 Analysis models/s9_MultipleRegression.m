% s9_MultipleRegression

% Second level analysis
clear all;close all hidden; clc

% where.where='I:\1 fMRI analysis'; 
where.where='D:\1 [Context-Memory] fMRI Data\1 fMRI analysis';
where.data_brain='D:\1 [Context-Memory] fMRI Data\1 MRI data'; where.data_beh=[where.where filesep '2 Behavioural data'];  where.data_behfolder=where.where;
% where.where='/Volumes/PENNYDISK/1 fMRI analysis';  where.data_brain='/Volumes/SANDISK/1 CONTEXT Brain data/1 MRI data'; where.data_behfolder='/Volumes/PENNYDISK/2 fMRI behaviour analysis';

% Requested analysis
log.specificsubjects={};

% BASIC MODEL TYPES ############################################
log.AnalysisType=' s4FullCardiacWithDerivAnts';
log.onsetsmodel='m_ci3_ContextItemNomotor';                  % FIRST LEVEL MODEL###############
for  o1=1:1 % Unused onsets models  
    % log.onsetsmodel='m_ci1_ContextItem';
% log.onsetsmodel='m_c4_ContextallItempresonly';
% log.onsetsmodel='m_ci10_ContextEventItemNomotor';
% log.onsetsmodel='m_ci11_ContextEvent';
% log.onsetsmodel='m_i3_ItemContextpresonly';
    % log.onsetsmodel='m_c1_Contextall';
% log.onsetsmodel='m_c2_Contextonly';
% log.onsetsmodel='m_c3_ContextallNooutcome';
% log.onsetsmodel='m_c4_ContextallItempresonly';
% log.onsetsmodel='m_c5_ContextonlyItempresonly';
% log.onsetsmodel='m_ci1_ContextItem';
% log.onsetsmodel='m_ci2_ContextonlyItem';
% log.onsetsmodel='m_ci4_ContextItemNomotorNooutcome';
% log.onsetsmodel='m_ci5_ContextItempmod';
% log.onsetsmodel='m_ci6_ContextItempmodWithmotor';
% log.onsetsmodel='m_ci7_ContextMemcatItemNomotor';
% log.onsetsmodel='m_ci8_ContextMemcatItempresonly';
% log.onsetsmodel='m_ci9_ContextMemPcatItemNomotor';
% log.onsetsmodel='m_ci12_ContextSimboxDisstickItemNomotor';
% log.onsetsmodel='m_i1_Item';
% log.onsetsmodel='m_i2_ItemNooutcome';
end
log.memtype='Hit';

for o1=1:1 % General settings and specifications
    
    % Subjects
    addpath(where.where);addpath(where.data_behfolder);
    log.w=load([where.data_brain filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    w.onsets=load([where.data_brain filesep 'Onsetslog_' log.firstlevelmodel]); % This filter is only applied first to catch subjects with bad onsets files, where this has yet to be marked in excel sheet
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog,log.specificsubjects,vertcat({'Subject' 'onsets'}, w.onsets.subjectsok),'onsets');
    [w.s w.s1 log.koshertable]=xlsread([where.where filesep 'i8_subjectdataok_secondlevel.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.firstlevelmodel);
    
    % Log
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' (' date ')  ' log.AnalysisType ' ' log.firstlevelmodel])
    errorlog=cell(1,1); e=1;
    log.firstlevelmodel=[log.onsetsmodel '_' log.memtype];
    
    % Specify memory type
    memdata=load([where.data_behfolder filesep '(31-May-2013) Memory scores.mat']); % Load memdata for covariates
    switch log.memtype
        case 'Hit'
            w.memtype{1}='Hit';
            w.memtype{2}='Miss';
        case 'Surehit'
            w.memtype{1}='Surehit';
            w.memtype{2}='NotSurehit';
        case 'Rem'
            w.memtype{1}='Rem';
            w.memtype{2}='NotRem';
        case 'Surerem'
            w.memtype{1}='Surerem';
            w.memtype{2}='NotSurerem';
        case 'Roc'
            w.memtype='Roc';
    end
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0;  disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location (brain): ' where.data_brain])
    disp('CHECK HERE: 1st and 2nd level models ---------');  disp(' ')
    disp(['             Model type:               ' log.AnalysisType])
    disp(['             First level model:      ' log.firstlevelmodel])
    disp(' '); input('Hit Enter to start      ')
    disp('=======================================================')
end

%%

log.secondlevelbranch=[where.data_brain filesep '2 Second level results' log.AnalysisType filesep];
log.secondlevelfolder=[where.data_brain filesep '2 Second level results' log.AnalysisType filesep log.firstlevelmodel];





%% Which contrasts, which covariates?





%%

log.scans_b4=[where.data_brain filesep];
log.scans_after=


\2 First level s4FullCardiacWithDeriv\m_ci3_ContextItemNomotor_Hit Contrasted\con_0005.img,1'

matlabbatch{1}.spm.stats.factorial_design.dir = {log.secondlevelfolder};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'RT';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;



matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = {
                                                            'D:\1 [Context-Memory] fMRI Data\1 MRI data\p01_CW\2 First level s4FullCardiacWithDeriv\m_ci3_ContextItemNomotor_Hit Contrasted\con_0005.img,1'
                                                            'D:\1 [Context-Memory] fMRI Data\1 MRI data\p03_EA\2 First level s4FullCardiacWithDeriv\m_ci3_ContextItemNomotor_Hit Contrasted\con_0005.img,1'
                                                            };
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = [1
                                                             2];
