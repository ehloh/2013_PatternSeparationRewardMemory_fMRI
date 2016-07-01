% First level: Specify & Estimate model
clear all;close all hidden; clc

% where.where='I:\1 fMRI analysis'; where.data_brain='D:\1 [Context-Memory] fMRI Data\1 MRI data'; where.data_beh=[where.where filesep '2 Behavioural data'];
% where.where='/Volumes/PENNYDISK/1 fMRI analysis';  where.data_brain='/Volumes/SANDISK/1 CONTEXT Brain data/1 MRI data'; where.data_beh=[where.where filesep '2 Behavioural data'];
% where.where='D:\1 [Context-Memory] fMRI Data\1 fMRI analysis';
where.where='D:\Dropbox\SCRIPPS\1 ContextMem fMRI';  where.data_brain='D:\1 [Context-Memory] fMRI Data\1 MRI data';  where.data_beh=[where.where filesep '2 Behavioural data'];


% Requested analysis
log.specificsubjects={}; 
% log.specificsubjects={'p01_CW'};
request.specify=1;
request.estimate=1;
%  ###### Analysis details ################
% request.AnalysisType=' s4WithDeriv';   log.prefix='s4wubf';   % SPM 
request.AnalysisType=' s4FullCardiacWithDerivAnts';   log.prefix='s4ubf';
% log.onsetsmodel='m_c1_Contextall';
% log.onsetsmodel='m_c2_Contextonly';
% log.onsetsmodel='m_c3_ContextallNooutcome';
% log.onsetsmodel='m_c4_ContextallItempresonly';
% log.onsetsmodel='m_c5_ContextonlyItempresonly';
% log.onsetsmodel='m_c6_ContextallItemNomem'; 
log.onsetsmodel='m_c7_ContexteventItempresonly';
% log.onsetsmodel='m_ci1_ContextItem';
% log.onsetsmodel='m_ci2_ContextonlyItem';
% log.onsetsmodel='m_ci3_ContextItemNomotor';
% log.onsetsmodel='m_ci4_ContextItemNomotorNooutcome';
% log.onsetsmodel='m_ci5_ContextItempmod';
% log.onsetsmodel='m_ci6_ContextItempmodWithmotor';
% log.onsetsmodel='m_ci7_ContextMemcatItemNomotor';
% log.onsetsmodel='m_ci8_ContextMemcatItempresonly';
% log.onsetsmodel='m_ci9_ContextMemPcatItemNomotor';
% log.onsetsmodel='m_ci10_ContextEventItemNomotor';
% log.onsetsmodel='m_ci11_ContextEvent';
% log.onsetsmodel='m_ci12_ContextSimboxDisstickItemNomotor';
% log.onsetsmodel='m_ci13_ContextStickItemNomotor'; 
% log.onsetsmodel='m_i1_Item';
% log.onsetsmodel='m_i3_ItemContextpresonly';
% log.onsetsmodel='m_i2_ItemNooutcome';
% log.onsetsmodel='m_i4_ItemContextevent';
% log.onsetsmodel='m_i5_ItemContextonly'; 
% log.onsetsmodel='m_i6_ItemMempmodContextpresonly';
% request.AnalysisType=' s4FullCardiacAnts';   log.prefix='s4ubf';  % FIR models #######
% log.onsetsmodel='m_io1_ItemMem';  
% log.onsetsmodel='m_co1_Contextall';    % Context doenst have a mem regressor 
log.memtype='Hit';   % Options: Hit\Surehit\Rem\Surerem
log.model=[log.onsetsmodel '_' log.memtype];

for o1=1:1 % General settings and specifications
        
    % Subjects
    addpath(where.where)
    log.w=load([where.data_brain filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    w.onsets=load([where.data_brain filesep 'Onsetslog_' log.model]); % Apply subject filter first according to validity of onsets file
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog,log.specificsubjects,vertcat({'Subject' 'onsets'}, w.onsets.subjectsok),'onsets');    
    
    % Settings that don't really change
    scan=load('i_scanningdetails.mat');
    log.smoothingsize=log.prefix(2);
    
    % Ants
    if isempty(strfind(log.prefix, 'w'))==1; disp('Non-normalized functionals chosen (for ANTS)'); log.scansuffix='ants'; else log.scansuffix=[]; end
    
    % Log
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' (' date ')  ' log.model])
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' ');disp('Requested analysis:')
    disp(request); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0;  disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location (brain): ' where.data_brain])
    disp(['Data location (behaviour): ' where.data_beh]); disp(' ')
    disp(['Analysis thread: ' request.AnalysisType])
    disp(['Smoothing size: ' num2str(log.smoothingsize)]); disp(' ')
    disp(['Model: ' log.model]); disp(' ')
    input('Hit Enter to start      ')
    disp('=======================================================')
    
    
    % Set up for 
    if sum(strcmp(log.onsetsmodel(3:4), {'io' 'co'}))>0 
        if ~isempty(strfind(request.AnalysisType, 'WithDeriv')), error('Wrong FL thread! Should not have derivs!'); end
        log.FIR_basis=1; 
        
        % Residuals are written using spm_write_residuals (spm12), which is called by spm_run_fmri_est. 
        %   Before you run estimal where you WANT the residuals written, make sure command is on in spm_run_fmri_est.
        %   Afterwards, MAKE SURE YOU TURN IT OFF. 
        % Files copied over from SPM12: spm_write_residuals,spm_data_hdr_write, spm_data_read, spm_data_write
        edit spm_run_fmri_est, input('spm_write residuals must be ON in spm_run_fmri_est! Continue?  '); 
        
    else log.FIR_basis=0;
    end
end

%% STEP 1: Specify model

if request.specify==1; % Specify: Format trial onsets & other regressors
    disp('STEP 1: Specifying model #########################')
    for o1=1:1 % Settings
        settings.firstlevelmodelspec.timing.units = 'secs';
        settings.firstlevelmodelspec.timing.RT = scan.TRms/1000;
        settings.firstlevelmodelspec.timing.fmri_t = 16;
        settings.firstlevelmodelspec.timing.fmri_t0 = 1;
        settings.firstlevelmodelspec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
        settings.firstlevelmodelspec.sess.regress = struct('name', {}, 'val', {});
        %         settings.firstlevelmodelspec.sess.multi_reg = {''};
        settings.firstlevelmodelspec.sess.hpf = 128;
        settings.firstlevelmodelspec.fact = struct('name', {}, 'levels', {});
        if log.FIR_basis  % 
%             log.FIRbasis={'Item' 2 1}; % Item regressors: 2s worth of candlesticks, 1 candlestick only 
            log.FIRbasis={'Context' 12 6}; % Context regressors: 2s worth of candlesticks, 6 candlestick only 
            disp(['FIR basis fxn: Type = ' log.FIRbasis{1} ', ' num2str(log.FIRbasis{2})  's epoch, ' num2str(log.FIRbasis{3}) ' candlesticks']), input('Check if stim type (item/context ) correct. Continue?'); 
            %
            settings.firstlevelmodelspec.bases.fir.length = log.FIRbasis{2};  % Length of window, in SECONDS 
            settings.firstlevelmodelspec.bases.fir.order = log.FIRbasis{3};   % % No. of FIR fxns
        else settings.firstlevelmodelspec.bases.hrf.derivs =[1 1]; % Standard HRF basis fxns
        end
        settings.firstlevelmodelspec.volt = 1;
        settings.firstlevelmodelspec.global = 'None';
        settings.firstlevelmodelspec.mask = {''};
        settings.firstlevelmodelspec.cvi = 'AR(1)';
%                 settings.firstlevelmodelspec.cond.sess.pmod.poly = 1; % Polynomials?
    end
    for s=1: log.n_subjs
%         try
            disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
%             wb.whereproc=[where.data_brain filesep log.subjects{s} filesep '1 Preprocessed' filesep 'Preproc functionals s' num2str(log.smoothingsize) log.scansuffix filesep];
            wb.whereproc=['F:\1 Context Mem study\1 All data\' log.subjects{s} '\1 Preprocessed\Preproc functionals s'  num2str(log.smoothingsize) log.scansuffix filesep];
            %
            wb.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' request.AnalysisType filesep];
            wb.wheremodel=[wb.where log.model ' Estimated' filesep];
            if isdir(wb.wheremodel)==0; mkdir(wb.wheremodel); end
            
            % SPECIFY MODEL -------------------
            matlabbatch{1}.spm.stats.fmri_spec= settings.firstlevelmodelspec;
            matlabbatch{1}.spm.stats.fmri_spec.dir = {wb.wheremodel}; % Specification files
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi ={[wb.where log.subjects{s} '_onsets_' log.model '.mat']};
            onsetsmatrix{s}=load([wb.where log.subjects{s} '_onsets_' log.model '.mat']);
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[wb.where log.subjects{s} '_reg_physiomovement.txt']};
            
            % Choose functional files
            f=spm_select('List', wb.whereproc, ['^' log.prefix '.*img$']); % Block 1
            wb.func=cellstr([repmat(wb.whereproc , size(f,1),1) f]);
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans =wb.func;
            files=wb.func;
            %
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            % Rename SPM file (Unidentified SPM.mat file exists only during Specification process (though still within the model's folder). After, it is immediately  renamed.
            eval('java.io.File([wb.wheremodel ''SPM.mat'']).renameTo(java.io.File([wb.wheremodel ''SPM_'' log.model  ''.mat'']));')
            matlabbatch=[];wb=[];
%         catch
%             errorlog{e}=['ERROR: Could not specify model for subject   ' log.subjects{s} ]; disp(errorlog{e}); e=e+1;
%         end
    end
end

%% STEP 2: Estimate model

% Estimate
if request.estimate==1
    disp('STEP 2: Estimate model ################################')
    for s=1:log.n_subjs
%         try
            disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
            wb.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' request.AnalysisType filesep];
            wb.wheremodel=[wb.where log.model ' Estimated' filesep];
            %
            f   = spm_select('List', wb.wheremodel, ['SPM_' log.model '.mat']);
            if isempty(f); error('Could not find Specified SPM.mat file to Estimate model. Has model been specified?'); end
            matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr([wb.wheremodel filesep f]);
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            %
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            matlabbatch=[];wb=[];
%         catch
%             errorlog{e}=['ERROR: Could not estimate model for subject   ' log.subjects{s} ]; disp(errorlog{e}); e=e+1;
%         end
    end
end

%% END

disp('======================================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
disp('Analysis completed (s7_Firstlevel_1SpecifyEstimate):'); disp(request)
disp(['No. of subjects: ' num2str(log.n_subjs)]);disp(' '); disp(log);
disp('Errors:');disp(errorlog');disp(' ')
disp('=======================================================')

diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s7_Firstlevel_1SpecifyEstimate)'), ' ',1);
end



% Reminder: Turn off residuals wrtite  
if  sum(strcmp(log.onsetsmodel(3:4), {'io' 'co'}))>0
    disp('############ WARNING ##############');  disp('############ WARNING ##############');  disp('############ WARNING ##############'); 
    disp('Turn OFF spm_write residuals in spm_run_fmri_est!!!  '), edit spm_run_fmri_est
end

