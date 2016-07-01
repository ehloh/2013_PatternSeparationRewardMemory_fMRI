% Second level analysis
clear all;close all hidden; clc; where.where='D:\Dropbox\SCRIPPS\1 ContextMem fMRI'; where.data_brain='D:\1 [Context-Memory] fMRI Data\1 MRI data'; where.data_beh=[where.where filesep '2 Behavioural data'];   where.data_behfolder='D:\Dropbox\SCRIPPS\2a ContextMem behaviour'; 

% Requested analysis
log.specificsubjects={};
%
log.AntsType='_Landmarks5';
log.AnalysisType=' s4FullCardiacWithDerivAnts';                              % MODEL TYPES ###################
log.covariates=[];   % 'subjects', 'memeffect'
% log.onsetsmodel='m_ci3_ContextItemNomotor';                  % FIRST LEVEL MODEL###############
% log.onsetsmodel='m_ci2_ContextonlyItem'; 
% log.onsetsmodel='m_ci1_ContextItem';
% log.onsetsmodel='m_c4_ContextallItempresonly';
% log.onsetsmodel='m_c6_ContextallItemNomem'; 
% log.onsetsmodel='m_c7_ContexteventItempresonly';
% log.onsetsmodel='m_ci10_ContextEventItemNomotor';
% log.onsetsmodel='m_ci11_ContextEvent';
% log.onsetsmodel='m_i3_ItemContextpresonly';
% log.onsetsmodel='m_i4_ItemContextevent'; 
% log.onsetsmodel='m_i5_ItemContextonly'; 
log.onsetsmodel='m_i6_ItemMempmodContextpresonly'; % log.secondlevelmodel='im_m1_2x2';
% log.onsetsmodel='m_ci13_ContextStickItemNomotor'; 
for  o1=1:1 % Unused onsets models  
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
log.firstlevelmodel=[log.onsetsmodel '_' log.memtype];
% log.secondlevelmodel='cm_m1_2x2';                 % SECOND LEVEL MODEL###############
% log.secondlevelmodel='co_m1_2x2';
% log.secondlevelmodel='com_m1_2x2x2';
% log.secondlevelmodel='em_m1_2x2';
% log.secondlevelmodel='im_m1_2x2';
log.secondlevelmodel='io_m1_2x2';
% log.secondlevelmodel='iom_m1_2x2x2';
% log.secondlevelmodel='o_outcomecontrasts';

for o1=1:1 % General settings and specifications
    
    % Subjects
    addpath(where.where);addpath(where.data_behfolder);
    log.w=load([where.data_brain filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    w.onsets=load([where.data_brain filesep 'Onsetslog_' log.firstlevelmodel]); % This filter is only applied first to catch subjects with bad onsets files, where this has yet to be marked in excel sheet
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog,log.specificsubjects,vertcat({'Subject' 'onsets'}, w.onsets.subjectsok),'onsets');
    
%     disp('Good mem subjects: DISABLED')
    [w.s w.s1 log.koshertable]=xlsread([where.where filesep 'i8_subjectdataok_secondlevel.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.firstlevelmodel);
    
    
    log.firstlevelmodel=[log.firstlevelmodel log.AntsType];

    
    % Log
    addpath([where.where filesep '5b Scripts - Set up model' filesep '1b Specify Second level models']) 
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' (' date ')  ' log.firstlevelmodel])
    errorlog=cell(1,1); e=1;
    
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
    disp(['             Second level model:  ' log.secondlevelmodel])
    disp(['             Covariates included:   ' log.covariates]);
    disp(' '); input('Hit Enter to start      ')
    disp('=======================================================')
end


%% (1) Second-level model specification + Estimation

disp('STEP 1: Specify 2nd-level model  ####################')

% Set up folders (2nd level overall, Specific 1st level. 2nd level folders made later)
cd(where.data_brain); cd ..; w.here=pwd; cd(where.where); w.mainresultsfolder=[w.here filesep '2 Second level results' log.AnalysisType filesep];
if isdir(w.mainresultsfolder)==0; mkdir(w.mainresultsfolder); end
if isdir([w.mainresultsfolder log.firstlevelmodel])==0; mkdir([w.mainresultsfolder log.firstlevelmodel]); end
log.secondlevelfolder=[w.mainresultsfolder log.firstlevelmodel filesep log.secondlevelmodel];

input('Continue to specify model?    ')

% Specify requested model
eval([' [ matlabbatch contrastfiles  log.secondlevelfolder] = model_' log.secondlevelmodel '(where.data_brain, log.secondlevelfolder, log.subjects, log.AnalysisType, log.firstlevelmodel, log.covariates, w.memtype, memdata);'])

% Run SPM batch (Model specification)
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];


% Estimate 2nd-level model
disp('STEP 2: Estimating model  ####################')
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[log.secondlevelfolder filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];
mkdir([log.secondlevelfolder filesep 'ROI'])

%% (2) Add additional sensible contrasts, if Factorial design

% error('pause')

if strcmp(log.secondlevelmodel(length(log.secondlevelmodel)-3:end), '_2x2')==1 || strcmp(log.secondlevelmodel(length(log.secondlevelmodel)-6:end), '_2x2_ad')==1 
    disp('STEP 3: Adding sensible contrasts for the 2x2 model  ####################')
%     ws=load([log.secondlevelfolder filesep 'SPM.mat']);

    % Compile instructions
    for o1=1:1 
        cons={'Negative effect of Similarity' [-1 -1 1 1];
            'Negative effect of Valence' [-1 1 -1 1];
            'Negative Interaction: Similarity x Valence' [-1 1 1 -1];
            'ME Sim_effect'          [1 1 0 0]; % Main effects (against baseline)
            'ME Dis_effect'          [0 0 1 1];
            'ME Val_effect'          [1 0 1 0];
            'ME Neu_effect'         [0 1 0 1];
            'SR_effect'                 [1 0 0 0];  % SR
            'SR_vs_others'          [3 -1 -1 -1];
            'SR_vs_SN'               [1 -1 0 0];
            'SR_vs_DR'               [1 0 -1 0];
            'SN_effect'                 [0 1 0 0];  % SN
            'SN_vs_others'          [-1 3 -1 -1];
            'SN_vs_SR'               [-1 1 0 0];
            'SN_vs_DN'               [0 1 0 -1];
            'DR_effect'                 [0 0 1 0]; % DR
            'DR_vs_others'          [-1 -1 3 -1];
            'DR_vs_DN'               [0 0 1 -1];
            'DN_effect'                 [0 0 0 1]; % DN
            'DN_vs_others'          [-1 -1 -1 3];
            'DN_vs_DR'               [0 0 -1 1]};
    end
    
    % Compile batch
    matlabbatch{1}.spm.stats.con.spmmat =  {[log.secondlevelfolder filesep 'SPM.mat']};
    matlabbatch{1}.spm.stats.con.delete = 0;
    for c=1: size(cons,1)
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = cons{c,1};
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec=cons{c,2};
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
%         matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = zeros(size(ws.SPM.xX.X,2),1);
%         matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec(1:4)=cons{c,2}';
    end
    
    
    % Identity contrast
    matlabbatch{1}.spm.stats.con.consess{c}.fcon.name = 'identity';
    matlabbatch{1}.spm.stats.con.consess{c}.fcon.convec = eye(4);
    matlabbatch{1}.spm.stats.con.consess{c}.fcon.sessrep = 'none';


    % Run SPM
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];
    
end


%% END

disp('======================================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
disp('Analysis completed (s8_Secondlevel)'); disp(['No. of subjects: ' num2str(log.n_subjs)])
if isempty(log.specificsubjects)==0;  disp('   Subset of subjects only:'); disp(log.subjects); end
disp(' ')
disp(['Data location (brain): ' where.data_brain])
disp(' ')
disp(['             Model type:               ' log.AnalysisType])
disp(['             First level model:  ' log.firstlevelmodel])
disp(['             Second level model:  ' log.secondlevelmodel])
disp(' ')
disp('Errors:')
disp(errorlog')
disp(' ')
disp('=======================================================')
cd(log.secondlevelfolder)
diary off
save([log.secondlevelfolder filesep 'details_2ndlevel'], 'log', 'contrastfiles'); % Save details in 2nd level folder
try % Notify researcher
    f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s8_Secondlevel)'), ' ',1);
end
