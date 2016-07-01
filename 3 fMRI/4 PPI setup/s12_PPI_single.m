% PPI - single (in SPM)
clear all;clc; % close all hidden; clc

where.where='J:\1 fMRI analysis'; 
% where.where='D:\1 [Context-Memory] fMRI Data\1 fMRI analysis';
where.exp_folder='D:\1 [Context-Memory] fMRI Data'; where.data_beh=[where.where '\2 Behavioural data'];where.data_brain=[where.exp_folder filesep '1 MRI data'];
% where.where='/Volumes/PENNYDISK/1 fMRI analysis'; where.exp_folder='/Users/EleanorL/Desktop/1 CONTEXT fmri data'; where.data_beh='/Volumes/PENNYDISK/1 fMRI analysis/2 Behavioural data';

% Requested analysis ###############
log.specificsubjects={}; % 
%     'p14_SJ'; 'p27_EW'; 'p04_JL'; 'p24_LL'; 'p16_TH'; 'p25_BS'; 'p08_AM'; 'p10_AB'; 'p26_MC'; 'p11_SS'; 'p21_SH'; 'p01_CW'}; % median split top, acc to dprime sim valfx

request.construct_PPIterm=1;
request.PPImodel_FirstLevel=0;
request.PPImodel_SecondLevel=0;
%
log.PPImodel_ComparisonType='ContextCompareSimxVal';
% log.PPImodel_ComparisonType='ContextCells';
%
log.regressor_famtype=1; % 1=Context onsets, 2=Context pmod, 3= ContextOnsetMem

% Second level setup ###############
% request.secondlevel.VOI='sph_SNVTA_R1_cmSRvSN'; 
% request.secondlevel.VOI='sph_SNVTA_R2_cmSRvSN'; 
request.secondlevel.VOI='sph_HPC_aL_cmSRvSN';
request.secondlevel.PsychVariable='DisN';
request.secondlevel.model={'onesampleT' {'PPI Pos'};};

% First level PPI setup ###############
instruc.VOIname={
% 'SNVTA_R1_cmSRvSN'; 'sph_SNVTA_R2_cmSRvSN'; 'sph_HPC_aL_cmSRvSN';
% 'SNVTA_R_cmSRvSN'; 
'HPC_aL_cmSRvSN';
    };

% Univariate model details ####################
log.AnalysisType=' s4FullCardiacWithDerivAnts'; log.memtype='Hit'; 
log.AntsType='_Landmarks5';
% log.onsetsmodel='m_c4_ContextallItempresonly';
log.onsetsmodel='m_ci3_ContextItemNomotor'; 
for o1=1:1 % Other models 
    % log.onsetsmodel='m_c4_ContextallItempresonly'; % Hit/Roc                 % FIRST LEVEL MODEL ---------------
    % log.onsetsmodel='m_ci3_ContextItemNomotor'; % Hit
    % log.onsetsmodel='m_ci5_ContextItempmod'; % Hit/Roc
% log.onsetsmodel='m_ci8_ContextMemcatItempresonly'; 
% log.onsetsmodel='m_ci7_ContextMemcatItemNomotor';
    % log.memtype='Hit';
    for o2=1:1 % Not so used onsets models
        % log.onsetsmodel='m_ci4_ContextItemNomotorNooutcome';
        % log.onsetsmodel='m_i1_Item';
        % log.onsetsmodel='m_i2_ItemNooutcome';
        % log.onsetsmodel='m_c2_Contextonly';
        % log.onsetsmodel='m_c3_ContextallNooutcome';
        % log.onsetsmodel='m_ci1_ContextItem';
        % log.onsetsmodel='m_ci2_ContextonlyItem';
    end
    % log.secondlevelmodel='cm_m1_2x2';                 % SECOND LEVEL MODEL  ---------------
    % log.secondlevelmodel='co_m1_2x2';
    % log.secondlevelmodel='im_m1_2x2';
    % log.secondlevelmodel='io_m1_2x2';
    % log.secondlevelmodel='iom_m1_2x2x2';
end

for o1=1:1 % General settings and specifications 
   
    % Settings that don't change much
    where.data_brain=[where.exp_folder filesep '1 MRI data']; addpath(where.where); addpath([where.where filesep '6 PPI setup' filesep 'Second level models']);
    log.w=load([where.exp_folder filesep '1 MRI data' filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    log.firstlevelmodel=[log.onsetsmodel '_' log.memtype log.AntsType];
    log.smoothingsize=num2str(log.AnalysisType(3));
    if isempty(strfind(log.AnalysisType, 'Ants'))==1; log.func_prefix=['s' num2str(log.AnalysisType(3)) 'wubf'];      log.funcfol=[];
    else log.func_prefix=['s' num2str(log.AnalysisType(3)) 'ubf'];   log.funcfol='ants';  end
    
    % Subjects
    w.onsets=load([where.exp_folder filesep '1 MRI data' filesep 'Onsetslog_' log.onsetsmodel '_' log.memtype]); % This filter is only applied first to catch subjects with bad onsets files, where this has yet to be marked in excel sheet
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog,[],vertcat({'Subject' 'onsets'}, w.onsets.subjectsok),'onsets');
    [w.s w.s1 log.koshertable]=xlsread([where.where filesep 'i8_subjectdataok_secondlevel.xlsx'],'PPIsingle'); 
    w.specsubs=log.specificsubjects;
    
    
    if request.PPImodel_FirstLevel || request.construct_PPIterm
        for v=1:size(instruc.VOIname,1) % Include subjects that are ok for ALL VOIs. To avoid unecessary exclusions, run 1 VOI at a time.
            [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,[log.firstlevelmodel ' ' instruc.VOIname{v}]);
            log.specificsubjects=log.subjects;
        end
    elseif request.PPImodel_SecondLevel
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,[log.firstlevelmodel ' ' request.secondlevel.VOI]);
        log.specificsubjects=log.subjects;
    end
    disp(['Subject exclusion:   No. of subjects included=' num2str(log.n_subjs)]);
    disp('Only subjects that are ok for ALL VOIs are included.'); disp('To avoid unecessary exclusions, run 1 VOI at a time'); input('OK?    ');
    instruc.log.subjects=log.subjects; instruc.log.n_subjs=log.n_subjs; log.specificsubjects=w.specsubs;
    
    % Log
    diary([where.data_brain filesep 'SPM logs' filesep 'Log ' mfilename '_' log.firstlevelmodel '- (' date ')  '])
    errorlog=cell(1,1); e=1;

    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location (brain): ' where.exp_folder]);
    disp(' '); disp('             -------------------  CHECK HERE  -------------------'); disp(' ')
    disp(['Univariate model: ' log.AnalysisType '      ' log.firstlevelmodel ]); disp(' ')
    disp('Requested analysis:'); disp(request); disp(' ');
%     disp(['Regressor type:     ' log.reg_famtype]); disp(' ');
    if request.construct_PPIterm==1 ||  request.PPImodel_FirstLevel==1
        disp('Constructing first-level models. VOIs:'); disp(instruc.VOIname);
    elseif request.PPImodel_SecondLevel==1
        disp(['VOI seed:              ' request.secondlevel.VOI]);
        disp(['Psych variable:      ' request.secondlevel.PsychVariable]);
        disp(['Second level stat:  ' request.secondlevel.model{1}  ' (first level: '  char(request.secondlevel.model{2}) ')']);
    end
    disp(' '); input('Hit Enter to start      ')
    disp('=======================================================')
end

%% (1) Set up 

for o1=1:1 % Set up PPI details  
    
    % Contrast type, VOI type, etc ####################
    switch log.regressor_famtype
        case 1;  % Context onsets
            log.reg_famtype='ContextOnset'; 
            log.reg_famprefix='co'; log.reg_famclass=1; 
            log.cellnames={'sr' 'SimR'; 'sn' 'SimN'; 'dr' 'DisR'; 'dn' 'DisN'};
        case 2; % Context memory pmod
            log.reg_famtype='ContextMem'; 
            log.reg_famprefix='cm'; 
            log.reg_famclass=2; 
            log.cellnames={'sr' 'SimRxCMem_Hit^1'; 'sn' 'SimNxCMem_Hit^1'; 'dr' 'DisRxCMem_Hit^1'; 'dn' 'DisNxCMem_Hit^1'};
        case 3; % Context Onset x Memory event
            log.reg_famtype='ContextOnsetMem'; 
            log.reg_famprefix='com'; 
            log.reg_famclass=1; 
            log.cellnames={'sr_L' 'SimR_L'; 'sr_H' 'SimR_H'; 'sn_L' 'SimN_L'; 'sn_H' 'SimN_H'; 'dr_L' 'DisR_L';  'dr_H' 'DisR_H'; 'dn_L' 'DisN_L'; 'dn_H' 'DisN_H'};
        otherwise; error('Invalid Contrast type selected. Context onsets, Context pmods, or Item?')
    end
    
    % PPI Model specification ####################
    switch log.PPImodel_ComparisonType
        case 'ContextCompareSimxVal'
            instruc.PPIweights={
%                      'FSimxVal'   {'sr'      log.reg_famclass     1;
%                                          'sn'      log.reg_famclass     -1;
%                                          'dr'      log.reg_famclass     -1;
%                                          'dn'      log.reg_famclass     1};
%                      'Sim-Dis'     {'sr'      log.reg_famclass     1;    % Factor effects
%                                          'sn'      log.reg_famclass     1;
%                                          'dr'      log.reg_famclass     -1;
%                                          'dn'      log.reg_famclass     -1};
%                      'Rew-Neu'   {'sr'      log.reg_famclass     1;    
%                                          'sn'      log.reg_famclass     -1;
%                                          'dr'      log.reg_famclass     1;
%                                          'dn'      log.reg_famclass     -1};
                     'SR-SN'       {'sr'      log.reg_famclass     1;      % Simple Valence effects
                                        'sn'      log.reg_famclass     -1};
%                      'DR-DN'       {'dr'      log.reg_famclass     1;
%                                         'dn'      log.reg_famclass     -1};
%                      'SR-DR'       {'sr'      log.reg_famclass     1;      % Simple Similarity effects
%                                          'dr'      log.reg_famclass     -1;};
%                      'SN-DN'       {'sn'      log.reg_famclass     1;
%                                          'dn'      log.reg_famclass     -1};
%                      'SR-others'   {'sr'      log.reg_famclass     1;
%                                          'sn'      log.reg_famclass     -1;
%                                          'dr'      log.reg_famclass     -1;
%                                          'dn'      log.reg_famclass     -1};
                     

                % [POSTERIOR HPC]  ---------------------------------------------
%                      'SN-SR'       {'sr'      log.reg_famclass     -1;      % Simple Valence effects
%                                         'sn'      log.reg_famclass     1};
                                         };
             %
             if sum(strcmp(log.reg_famtype, {'ContextOnset';'ContextMem'}))==0
                 error('Wrong type of regressor. Only context regressors (onset or mem) allowed!')
             end

        case 'ContextCells'
            instruc.PPIweights={
                     'SimR'         {'sr'      log.reg_famclass     1;
                                         'sn'      log.reg_famclass     0;
                                         'dr'      log.reg_famclass     0;
                                         'dn'      log.reg_famclass     0};
                     'SimN'         {'sr'      log.reg_famclass     0;
                                         'sn'      log.reg_famclass     1;
                                         'dr'      log.reg_famclass     0;
                                         'dn'      log.reg_famclass     0};
                     'DisR'         {'sr'      log.reg_famclass     0;
                                         'sn'      log.reg_famclass     0;
                                         'dr'      log.reg_famclass     1;
                                         'dn'      log.reg_famclass     0};
                     'DisN'         {'sr'      log.reg_famclass     0;
                                         'sn'      log.reg_famclass     0;
                                         'dr'      log.reg_famclass     0;
                                         'dn'      log.reg_famclass     1};
                                    };
            %
            if sum(strcmp(log.reg_famtype, {'ContextOnset';'ContextMem'}))==0
                 error('Wrong type of regressor. Only context regressors (onset or mem) allowed!')
            end

        % DISUSED PPI models (non-productive) ###########
        case 'ContextMemLevs';
           instruc.PPIweights={
%                'SimR_mfx'      {'sr_H'     log.reg_famclass     1;
%                                         'sr_L'     log.reg_famclass     -1;};
%                'SimN_mfx'      {'sn_H'     log.reg_famclass     1;
%                                         'sn_L'     log.reg_famclass     -1;};
%                'DisR_mfx'       {'dr_H'     log.reg_famclass     1;
%                                         'dr_L'     log.reg_famclass     -1;};
%                'DisN_mfx'       {'dn_H'     log.reg_famclass     1;
%                                         'dn_L'     log.reg_famclass     -1;};
                };
            %
            if sum(strcmp(log.reg_famtype, {'ContextOnsetMem';}))==0
                 error('Wrong type of regressor. Only contextxmem regressors (com) allowed!')
            end
        otherwise; error('Invalid Context model (comparison type) specified');
    end
    
    
     % -------------------------
     % Set up PPI names
    log.ppi_names=cell(size(instruc.VOIname,1)*size(instruc.PPIweights,1),3);
    for v=1:size(instruc.VOIname,1)
        for w=1:size(instruc.PPIweights,1);
            k=(v-1)*size(instruc.PPIweights,1)+w;
            log.ppi_names{k,1}=[log.reg_famprefix 'Fam_' instruc.VOIname{v} '_psy_' instruc.PPIweights{w,1}];
            log.ppi_names{k,2}=instruc.VOIname{v};
            log.ppi_names{k,3}=instruc.PPIweights{w,1};
        end
    end
    
    % Misc setup
    if request.PPImodel_SecondLevel==1
        log.PPImodelname=[log.reg_famprefix 'Fam_' request.secondlevel.VOI '_psy_' request.secondlevel.PsychVariable];
        where.secondlevelres=[where.exp_folder filesep '2 Second level results' log.AnalysisType filesep log.firstlevelmodel filesep log.PPImodelname filesep];
        if isdir(where.secondlevelres)==0; mkdir(where.secondlevelres); mkdir([where.secondlevelres filesep 'ROI']); mkdir([where.secondlevelres filesep 'ROI' filesep 'ROI analysis MarsBar']);         end
    end
end

% Displace interface ####################################################
if request.construct_PPIterm==1 || request.PPImodel_FirstLevel==1
    disp('---------------------------- REQUESTED PPI MODEL -------------------------------------'); disp(' ')
    disp(['PPI model name: ' log.PPImodel_ComparisonType])
    disp(['Regressor type:  ' log.reg_famtype]); disp(' ')
    disp('Contrast weights for all requested psychological comparisons:'); disp(' '); for i=1:size(instruc.PPIweights,1); disp([num2str(i) ')  ' instruc.PPIweights{i,1}]); disp(instruc.PPIweights{i,2});    end; disp(' ')
    disp('PPI names: ');  for i=1:size(log.ppi_names,1); disp(['  ' log.ppi_names{i}]); end; disp(' ')
    input('OK?    ');
end
% ###############################################################

%% (2) Construct PPI term

% error

if request.construct_PPIterm
    disp('########### Generating PPI terms for each VOI x Psychological comparison ###############')

    for o1=1:1 % Delete existing VOI terms? 
        log.delete_Old_PPIterms=0;
        if log.delete_Old_PPIterms
            input('User requested deleting ALL existing old PPI terms from 1st level folders. Proceed?')
            disp('Deleting old PPI terms --------')
            for s=1:log.n_subjs 
                f=spm_select('List',[where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted'],'^PPI.*.mat$');
                for i=1:size(f,1)
                    delete([where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep f(i,1:strfind(f(i,:), '.mat')+3)]);
                end
            end      
        end
    end
    for o1=1:1 % Format instructions for PPI contrast weights 
        
        % Identify which condition no. (indexed Sess.U) for each design-cell/condition (Sample 1st subject)
        p=load([where.data_brain filesep log.subjects{1} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep 'SPM.mat']);
        
%                 p=cd([where.data_brain filesep log.subjects{1} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep])
%                 
%                 p=pwd
%                 f=spm_select('List', p, 'SPM.mat')
%                 'SPM.mat']);
%         
        
        if isempty(strfind(log.reg_famtype, 'Context'))==0 % Context type regressors
            p.condnames=cell(length(p.SPM.Sess.U),2);
            for i=1:length(p.SPM.Sess.U);
                p.condnames{i,1}=p.SPM.Sess.U(i).name{1};
                if length(p.SPM.Sess.U(i).name)>1
                    p.condnames{i,2}=p.SPM.Sess.U(i).name{2};
                end
            end
        else
            error('Error - haven''t said what to do with Item regressos')
        end
        
        % Write condition #s to the instructions
        for j=1:size(instruc.PPIweights,1)
            for i=1:size(instruc.PPIweights{j,2},1)
                if sum(    strcmp(    p.condnames(:,log.reg_famclass) ,char(log.cellnames{find(strcmp(log.cellnames(:,1),instruc.PPIweights{j,2}{i,1})),2})    )  )==1
                    instruc.PPIweights{j,2}{i,1}=find(strcmp(    p.condnames(:,log.reg_famclass) ,char(log.cellnames{find(strcmp(log.cellnames(:,1),instruc.PPIweights{j,2}{i,1})),2})    ) );  % Condition number. 
                    % Condition name: instruc.PPIweights{j,2}{i,1}=log.cellnames(strcmp(    p.condnames(:,log.reg_famclass) ,char(log.cellnames{find(strcmp(log.cellnames(:,1),instruc.PPIweights{j,2}{i,1})),2})    ) ,2);
                    % OLD code: instruc.PPIweights{j,2}{i,1}=find(strcmp(p.condnames(:,log.reg_famclass),char(log.cellnames{find(strcmp(log.cellnames(:,log.reg_famclass),instruc.PPIweights{j,2}{i,1})),2})));
                else
                    error('Error formatting instructions for PPI weights. Instructions don''t match Sess.U condition names in 1st level SPM.mat')
                end
            end
            instruc.PPIweights{j,2}=cell2mat(instruc.PPIweights{j,2});
        end
        
        % Notes re contrast weights for PPI term generation
        %         Matrix  of  input  variables  and  contrast weights. This is an [nconditions x 3] matrix. The first column
        %         indexes SPM.Sess.U(i). The second column (j) indexes the name of the input or cause, see SPM.Sess.U(i).name{j}.
        %         The third column is the contrast weight.  Unless there are parametric effects the second column will generally be a 1.
        %         SPM.Sess.U(i) is regressor i. pmods are subsumed in the
        %         regressors, such that regressor 1, 2, ... i are all main regressors, & pmods are specified by the j
        disp('Contrast weights for all requested psychological comparisons:'); for i=1:size(instruc.PPIweights,1); disp(['[' num2str(i) ']  ' instruc.PPIweights{i,1}]); disp(instruc.PPIweights{i,2}); disp(' '); end
    end
    
    % Generate PPI terms according to instructions
    for v=1:size(instruc.VOIname,1) % For each VOI
        disp(['VOI #' num2str(v) ':  ' instruc.VOIname{v}  '  -------------------'])
        for w=1:size(instruc.PPIweights,1); % For each Psychological comparison
            disp(['Psych Comparison #' num2str(w) ':  ' instruc.PPIweights{w,1}  '  --------'])
            k=(v-1)*size(instruc.PPIweights,1)+w;
            
            for s=1:log.n_subjs
                disp(['Subject no. ' num2str(s) ' (' log.subjects{s} ') ----------------'])
                ws.wheremod=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep] ;
                try
                    matlabbatch{1}.spm.stats.ppi.spmmat = {[ws.wheremod 'SPM.mat']};
                    if s==1;matlabbatch{1}.spm.stats.ppi.disp = 1;
                    else matlabbatch{1}.spm.stats.ppi.disp = 0;
                    end
                    %
                    matlabbatch{1}.spm.stats.ppi.name =log.ppi_names{k,1};
                    matlabbatch{1}.spm.stats.ppi.type.ppi.voi = {[ws.wheremod 'VOI_' instruc.VOIname{v} '_1.mat']};
                    matlabbatch{1}.spm.stats.ppi.type.ppi.u =instruc.PPIweights{w,2};
                    %
                    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                    matlabbatch=[];
                    
                    % Append PPI term to add details
                    PPI=[];
                    load([ws.wheremod 'PPI_' log.ppi_names{k,1} '.mat']);
                    PPI.details.SeedVOI= instruc.VOIname{v};
                    PPI.details.PsychCond_name=instruc.PPIweights{w,1};
                    PPI.details.PsychCond_weights=instruc.PPIweights{w,2};
                    PPI.details.PsychCond_AssumedRegs=p.condnames; %
                    save([ws.wheremod 'PPI_' log.ppi_names{k,1} '.mat'], 'PPI');
                    
                    % Create directory for this PPI & move the PPI term to it
                    mkdir([ws.wheremod 'PPI ' log.ppi_names{k,1}])
                    movefile([ws.wheremod 'PPI_' log.ppi_names{k,1} '.mat'], [ws.wheremod 'PPI ' log.ppi_names{k,1} filesep 'PPI_' log.ppi_names{k,1} '.mat'])
                catch
                    errorlog{e,1}=['Couldn''t generate PPI term subject  ' log.subjects{s} '  - ' instruc.VOIname{v} '  - ' instruc.PPIweights{w,1}]; disp(errorlog{e,1}); e=e+1;
                end
                ws=[];
            end
        end
    end
end

%% (3) Set up first-level for PPI models (Specify, Estimate, Contrast)

if request.PPImodel_FirstLevel
    request.Spec=1;
    request.Est=1; 
    request.Contrast=1;
    
    if request.Spec
        disp('Specifying PPI models ################################')
        for o1=1:1 % Specification - Settings for first-level models
            scan=load([where.where filesep '5b Scripts - Set up model' filesep 'i_scanningdetails.mat']);
            settings.firstlevelmodelspec.timing.units = 'secs';
            settings.firstlevelmodelspec.timing.RT = scan.TRms/1000;
            settings.firstlevelmodelspec.timing.fmri_t = 16;
            settings.firstlevelmodelspec.timing.fmri_t0 = 1;
            settings.firstlevelmodelspec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
            settings.firstlevelmodelspec.sess.regress = struct('name', {}, 'val', {});
            %         settings.firstlevelmodelspec.sess.multi_reg = {''};
            settings.firstlevelmodelspec.sess.hpf = 128;
            settings.firstlevelmodelspec.fact = struct('name', {}, 'levels', {});
            settings.firstlevelmodelspec.bases.hrf.derivs =[1 1];
            settings.firstlevelmodelspec.volt = 1;
            settings.firstlevelmodelspec.global = 'None';
            settings.firstlevelmodelspec.mask = {''};
            settings.firstlevelmodelspec.cvi = 'AR(1)';
            %         settings.firstlevelmodelspec.cond.sess.pmod.poly = 1; % Polynomials?
        end
        for v=1: size(instruc.VOIname,1) % For each VOI
            for w=1: size(instruc.PPIweights,1); % For each Psychological comparison
                k=(v-1)*size(instruc.PPIweights,1)+w;
                disp(['Specifying models for PPI #' num2str(k) ': ' log.ppi_names{k,1} '------------------------------------------------------'])
                
                for s=1:log.n_subjs
%                     try
                        disp(['Subject ' num2str(s) '   (' log.subjects{s} ') ------------'])
                        ws.whereproc=[where.data_brain filesep log.subjects{s} filesep '1 Preprocessed' filesep 'Preproc functionals s' num2str(log.smoothingsize) log.funcfol filesep];
                        ws.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep];
                        ws.wherePPImodel=[ws.where log.firstlevelmodel ' Contrasted' filesep 'PPI ' log.ppi_names{k,1} filesep];
                        
                        % Typical model specification
                        matlabbatch{1}.spm.stats.fmri_spec= settings.firstlevelmodelspec;
                        matlabbatch{1}.spm.stats.fmri_spec.dir = {ws.wherePPImodel};
                        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[ws.where log.subjects{s} '_reg_physiomovement.txt']};
                        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
                        f=spm_select('List', ws.whereproc, ['^' log.func_prefix '.*img$']); % Functional files
                        matlabbatch{1}.spm.stats.fmri_spec.sess.scans =cellstr([repmat(ws.whereproc , size(f,1),1) f]);
                        
                        % PPI model specification
                        ws.p=load([ws.wherePPImodel 'PPI_' log.ppi_names{k,1} '.mat']);
                        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = 'PPI';    % Interaction
                        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = [ws.p.PPI.ppi];
                        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = 'P';      % Psychological
                        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = [ws.p.PPI.P];
                        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name = 'Y';      % Extracted signal
                        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val = [ws.p.PPI.Y];
                        
                        %
                        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                        matlabbatch=[];ws=[];
%                     catch
%                         errorlog{e}=['ERROR: Could not specify model for subject   ' log.subjects{s} ' - ' log.ppi_names{k,1}]; disp(errorlog{e}); e=e+1;
%                     end
                end
            end
        end
    end
    
    if request.Est
        disp('Estimating PPI models ################################')
        for v=1: size(instruc.VOIname,1)
            for w=1: size(instruc.PPIweights,1);
                k=(v-1)*size(instruc.PPIweights,1)+w;
                disp(['Estimating model for PPI #' num2str(k) ': ' log.ppi_names{k,1} '------------------------------------------------------'])
                
                for s=1:log.n_subjs
                    try
                        disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
                        ws.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep];
                        ws.wheremodel=[ws.where log.firstlevelmodel ' Contrasted' filesep 'PPI ' log.ppi_names{k,1}];
                        %
                        f   = spm_select('List', ws.wheremodel, 'SPM.mat');
                        if isempty(f); error('Could not find Specified SPM.mat file to Estimate model. Has model been specified?'); end
                        matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr([ws.wheremodel filesep f]);
                        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
                        %
                        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                        matlabbatch=[];wb=[];
                    catch
                        errorlog{e}=['ERROR: Could not estimate model for subject   ' log.subjects{s} ' - ' log.ppi_names{k,1}]; disp(errorlog{e}); e=e+1;
                    end
                end
            end
        end
    end
    
    if request.Contrast
        disp('Running Contrasts PPI models ################################')
        for v=1:size(instruc.VOIname,1)
            for w=1:size(instruc.PPIweights,1);
                k=(v-1)*size(instruc.PPIweights,1)+w;
                disp(['Running contrasts for PPI #' num2str(k) ': ' log.ppi_names{k,1} '------------------------------------------------------'])
                
                for s=1:log.n_subjs
                    try
                        disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
                        ws.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep];
                        ws.wherePPImodel=[ws.where log.firstlevelmodel ' Contrasted' filesep 'PPI ' log.ppi_names{k,1} filesep];
                        %
                        matlabbatch{1}.spm.stats.con.spmmat = {[ws.wherePPImodel 'SPM.mat']};
                        matlabbatch{1}.spm.stats.con.delete = 0;
                        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'PPI Pos';
                        matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [1];
                        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
                        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'PPI Neg';
                        matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [-1];
                        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
                        %
                        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                        matlabbatch=[]; ws=[];
                    catch
                        errorlog{e}=['ERROR: Could not run contrasts for subject   ' log.subjects{s} ' - ' log.ppi_names{k,1}]; disp(errorlog{e}); e=e+1;
                    end
                end
            end
        end
    end
end            
            

%% (4) Second level (Specific single VOI & Psych variable only)

if request.PPImodel_SecondLevel
    disp(['Running 2nd level model: ' request.secondlevel.model{1} ' ################################'])
    eval(['[ batch] = '  request.secondlevel.model{1} '(where, log, request.secondlevel.model{2});'])
    
    % Mark details if VOI is a sphere
%     if strcmp(request.secondlevel.VOI(1:4), 'sph_')
%         load([where.secondlevelres filesep 'details_2ndlevel.mat'])
%         load([where.data_brain filesep log.subjects{1} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep 'VOI_sphere_details.mat']);
%         if sum(strcmp(VOI_spherespex(:,2),request.secondlevel.VOI(5:end)))~=0
%             log.VOIspherespex= {VOI_spherespex{find(strcmp(VOI_spherespex(:,2),request.secondlevel.VOI(5:end))),3} VOI_spherespex{find(strcmp(VOI_spherespex(:,2),request.secondlevel.VOI(5:end))),4}(1) VOI_spherespex{find(strcmp(VOI_spherespex(:,2),request.secondlevel.VOI(5:end))),4}(2) VOI_spherespex{find(strcmp(VOI_spherespex(:,2),request.secondlevel.VOI(5:end))),4}(3)};
%             save([where.secondlevelres filesep 'details_2ndlevel.mat'], 'log')
%         else
%             error('VOI appears to be a sphere, but could not find details of sphere in 1st level VOI_sphere_details');
%         end        
%     end 
end

%% END

disp('======================================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
disp(['GLM model: ' log.AnalysisType '      ' log.firstlevelmodel]); disp(' ')
disp('Analysis completed:'); disp(request);
disp('Errors:'); disp(errorlog)
disp('=======================================================')


diary off
try % Notify researcher
    f_sendemail('kurzlich', ['Analysis batchscript is complete (' mfilename ')'], ' ',1);
end
    




