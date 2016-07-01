% First level checks: Does each subject show the results that we are
% looking for?
clear all;close all hidden; clc

where.where='I:\1 fMRI analysis';
% where.data_brain='C:\Users\eloh\Desktop\0 TEST Context fMRI';
where.data_brain='C:\Users\eloh\Desktop\1 [Context-Memory] fMRI Data\1 MRI data';              %  #### REAL DATA ###
% where.where='/Volumes/PENNYDISK/1 fMRI analysis';
% where.data_brain='/Volumes/SANDISK/1 CONTEXT Brain data/1 MRI data';
where.data_beh=[where.where filesep '2 Behavioural data'];
log.w=load([where.data_brain filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
scan=load([where.where filesep 'i_scanningdetails.mat']);

% Requested analysis
log.specificsubjects={'p01_CW'}; % BLANK to process all subjects

% ############# Which first-level Model #######################
% log.onsetsmodel='m_c1_Contextall';
% log.onsetsmodel='m_c2_Contextonly';
% log.onsetsmodel='m_c3_ContextallNooutcome';
% log.onsetsmodel='m_ci1_ContextItem';
% log.onsetsmodel='m_ci2_ContextonlyItem';
log.onsetsmodel='m_ci3_ContextItemNomotor';
% log.onsetsmodel='m_ci4_ContextItemNomotorNooutcome';
% log.onsetsmodel='m_ci5_ContextItempmod';
% log.onsetsmodel='m_i1_Item';
% log.onsetsmodel='m_i2_ItemNooutcome';
log.memtype='Hit';   % Options: Hit\Surehit\Rem\Surerem
log.firstlevelmodel=[log.onsetsmodel '_' log.memtype];
% log.secondlevelmodel='c_m1_2x2x2';

% ############# Which checks to perform? #######################
log.checkingscript='c7_checkfactorial_factorresults';
log.ContextorItem=1;

for o1=1:1 % General settings and specifications
        
    % Subjects
    allsubjects=0;
    if allsubjects==0
        [w.s w.s1 log.koshertable]=xlsread([where.where filesep 'i8_subjectdataok_secondlevel.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.firstlevelmodel);
    else
        w.onsets=load([where.data_brain filesep 'Onsetslog_' log.firstlevelmodel]); % Apply subject filter first according to validity of onsets file
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog,[],vertcat({'Subject' 'onsets'}, w.onsets.subjectsok),'onsets');
    end
    
    %
    errorlog=cell(1,1); e=1;
    addpath([where.where filesep '5 Checks'])
    
    % Interface
    disp('=======================================================')
    w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' ')
    disp('Check fMRI results after first-level analysis')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0
        disp('   Subset of subjects only:')
        disp(log.specificsubjects)
    end
    disp(' ')
    disp(['Data location (brain): ' where.data_brain])
    disp(' ')
    disp(['Model: ' log.firstlevelmodel])
    disp(' ')
    input('Hit Enter to start      ')
    disp('=======================================================')
end

% Checks to 1st level: coverage (p<.99), visual/motor effects? subsequent
% mem effects, valence effects, similarity effects?



%% Settings for viewing results

% General settings
settings.pthresh=0.05;
settings.multiplecompar='none'; % 'none' or 'FWE'


%% Load details from all subject (contrasts to run, in checks)

error('stop here - bug in compiling contrasts. need to read the cntrast #, not the contrast weights')

checks=cell(log.n_subjs,1); disp('Compiling contrasts from all subjects')
for s=1:log.n_subjs
    ws=load([where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.firstlevelmodel ' Contrasted' filesep 'SPM.mat']);
    eval(['[subjchecks{s}] =' log.checkingscript ' (ws.SPM, log.ContextorItem, log.memtype);'])
end

n_checks=size(subjchecks{1},1);
checks=cell(n_checks,1);
for c=1:n_checks
    for s=1:log.n_subjs
        checks{c}{s,1}=subjchecks{s}{c,1};
        checks{c}{s,2}=subjchecks{s}{c,2};
    end
end

error('stop here')

%% View results (Execution) + Record assessment

% Standard settings
spm('defaults','FMRI'); spm_jobman('initcfg');
for o1=1:1 % Batch results report 
    defaultreport.units = 1;
    defaultreport.print = true;
    defaultreport.conspec.threshdesc =settings.multiplecompar; % 'FWE';
    defaultreport.conspec.thresh = settings.pthresh;
    defaultreport.conspec.extent = 0;
    defaultreport.conspec.mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
end

% Display
res=cell(log.n_subjs, n_checks);
for c=1:n_checks
    s=1;
    
    while s~=999
        % Load up requested result
        matlabbatch{1}.spm.stats.results=defaultreport;
        matlabbatch{1}.spm.stats.results.spmmat = {[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.firstlevelmodel ' Contrasted' filesep 'SPM.mat']};
        matlabbatch{1}.spm.stats.results.conspec.titlestr = checks{c}{s,1}; % Load title
        matlabbatch{1}.spm.stats.results.conspec.contrasts = checks{c}{s,2}; % Load contrasts
        spm_jobman('run' , matlabbatch);
        
        % Request & record inputs
        disp('----------------------------------------------------------------------------------')
        disp(['CURRENT COMPARISON:          ' checks{c}{s}  '     -     ' log.subjects{s} '   (' num2str(s) ')'])
        % disp('[SCROLL THRU] Enter to go to next subject; # to request specific subject; 999 to go to next contrast/comparison')
        disp('[Progress] Enter to go to next subject; 999 to go to next contrast/comparison')
        disp('[Log res] Log results: 1 if OK, else type reason')
        f=input('Next request:   ','s');
        
        % Progress
        if isnumeric(str2num(f)) && str2num(f)==999
            s=999;
        elseif isempty(f)==1
            s=s+1;
        else
            res{s,c}=f;
            s=s+1;
        end
        
        matlabbatch=[];
    end
    
end

%% END

disp('##########################################')
disp('DONE')
disp('##########################################')

disp(res) % Display results - copy to excel file