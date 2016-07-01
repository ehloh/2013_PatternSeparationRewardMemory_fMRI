% PPI prep: Set up the necessary contrasts
clear all;close all hidden; clc

where.where='J:\1 fMRI analysis';

where.where='D:\1 [Context-Memory] fMRI Data\1 fMRI analysis';
where.exp_folder='D:\1 [Context-Memory] fMRI Data'; where.data_beh=[where.where '\2 Behavioural data'];where.data_brain=[where.exp_folder filesep '1 MRI data'];
% where.where='/Volumes/PENNYDISK/1 fMRI analysis'; where.exp_folder='/Users/EleanorL/Desktop/1 CONTEXT fmri data'; where.data_beh='/Volumes/PENNYDISK/1 fMRI analysis/2 Behavioural data';

% Requested analysis
log.specificsubjects={'p22_EK'}
    
% 'p23_IS';'p24_LL';'p25_BS';'p26_MC';'p27_EW';'p28_CN';'p29_ET';}
%
request.extract_firstlevelVOI=1;
request.VOI_individuals=1; % 0= Group, 1 = Individual voi
request.VOI_maskwFLmask=0;
request.change_fxnscan_address=0;
%
log.regressor_famtype=1; % 1=Context onsets, 2=Context pmod, 3= ContextOnsetMem
%
request.contrast_Fomnibus=0; % If Ants, must adjust specially!
request.contrast_SpecificComparisons=0; % NOT for ants



% VOI details (seed region from which to read neural signal) ####################
request.extract_FLvoi_ImgMask={
%     'sph_SNVTA_R1_cmSRvSN'; 'sph_SNVTA_R2_cmSRvSN'; 'sph_HPC_aL_cmSRvSN';
% 'SNVTA_R_cmSRvSN'; 
'HPC_aL_cmSRvSN';
% 'SNVTA_Lpeak_coSRvSN'; 'HPC_pLpeak_cmSNvSR';
    }; 
request.extract_FLvoi_PeakSphere_Individual={
%     'HPC_aL_cmSRvSN';
    };
request.extract_FLvoi_PeakSphere_Group={
%     'HPC_aL_cmSRvSN'            [-27.4   5.5    9]       2;
    }; 

% Model details --------------
log.AnalysisType=' s4FullCardiacWithDerivAnts';
log.AntsType='_Landmarks5';
 log.onsetsmodel='m_ci3_ContextItemNomotor';log.memtype='Hit';  
 log.secondlevelmodel='cm_m1_2x2';                 
for o1=1:1 % Other models
    % log.onsetsmodel='m_c4_ContextallItempresonly'; % Hit/Roc                 % FIRST LEVEL MODEL ---------------
    % log.onsetsmodel='m_ci3_ContextItemNomotor'; % Hit
    % log.onsetsmodel='m_ci5_ContextItempmod'; % Hit/Roc
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
    addpath(where.where)
    where.data_brain=[where.exp_folder filesep '1 MRI data']; addpath(where.where)
    log.firstlevelmodel=[log.onsetsmodel '_' log.memtype log.AntsType];
    log.w=load([where.exp_folder filesep '1 MRI data' filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    log.func_prefix=['s' num2str(log.AnalysisType(3)) 'wubf'];
   
    % Subjects
    w.onsets=load([where.exp_folder filesep '1 MRI data' filesep 'Onsetslog_' log.onsetsmodel '_' log.memtype]); % This filter is only applied first to catch subjects with bad onsets files, where this has yet to be marked in excel sheet
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog,[],vertcat({'Subject' 'onsets'}, w.onsets.subjectsok),'onsets');
    [w.s w.s1 log.koshertable]=xlsread([where.where filesep 'i8_subjectdataok_secondlevel.xlsx']); 
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable, log.firstlevelmodel);
    instruc.log.subjects=log.subjects; instruc.log.n_subjs=log.n_subjs; 
    
    % Model details (Full model, 1st & 2nd levels)
    where.model=[where.exp_folder filesep '2 Second level results' log.AnalysisType filesep log.firstlevelmodel filesep log.secondlevelmodel filesep];
%     log.modeldetails=load([where.model 'details_2ndlevel.mat']);
    p=load([where.data_brain filesep log.subjects{1} filesep '2 First level' log.AnalysisType filesep log.onsetsmodel '_' log.memtype ' Contrasted' filesep 'SPM.mat']);
    
    % Log
    diary([where.data_brain filesep 'SPM logs' filesep 'Log ' mfilename '_' log.firstlevelmodel '-' log.secondlevelmodel ' (' date ')  '])
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location (brain): ' where.exp_folder])
    disp(' '); disp('Requested analysis:'); disp(' '); disp(request);
    disp(' '); disp(['GLM model: ' log.AnalysisType '      ' log.firstlevelmodel '     ' log.secondlevelmodel]); disp(' ')
    switch request.VOI_individuals; case 1; disp(' '); disp('VOIs extracted for different imgs/coordinates for each subject'); case 0;  disp(' '); disp('Extracting GROUP vois - same for all subjects'); end
    input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%% (1) Prep: Change address, implement omnibus & comparison contrasts

for o1=1:1 % Set up instructions depending on which type of model/regressor
    
    % Contrast/regressor type, VOI type, etc
    switch log.regressor_famtype
        case 1; log.reg_famtype='ContextOnset'; log.reg_famprefix='co';
        case 2; log.reg_famtype='ContextMem'; log.reg_famprefix='cm';
        case 3; log.reg_famtype='ContextOnsetMem'; log.reg_famprefix='com';
        otherwise; error('Invalid Contrast type selected. Context onsets, Context pmods, or Item?')
    end
    
    % (A) Change address of functional scans ###########################
%     request.change_fxnscan_address_newaddress='[where.data_brain filesep log.subjects{s} filesep ''1 Preprocessed'' filesep ''Preproc functionals '' log.func_prefix(1:2) filesep]';
request.change_fxnscan_address_newaddress='[where.data_brain filesep log.subjects{s} filesep ''1 Preprocessed'' filesep ''Preproc functionals '' log.func_prefix(1:2) filesep]';
    
    % (B) Implement Omnibus F-tests ######################################################
    %       Instruc: New contrast name, Regressor names, Regressor Types
    switch log.reg_famtype
        case 'ContextOnset'
            instruc.Fomni={'ContextOnset'       {'SimR';'SimN';'DisR';'DisN'}       1;};
        case 'ContextMem'
            instruc.Fomni={'ContextMem'       {['SimRxCMem_' log.memtype];['SimNxCMem_' log.memtype];['DisRxCMem_' log.memtype];['DisNxCMem_' log.memtype]}       2;};
        case 'ContextOnsetMem' % [SimxValxMemcat models]
            instruc.Fomni={'ContextOnsetMem'       {'SimR_L'; 'SimR_H';'SimN_L';'SimN_H';'DisR_L';'DisR_H';'DisN_L';'DisN_H'}       1;};
    end
    
    % (C) Specific comparisons ######################################################
    %           1=Context onset, 2=Regressor-conditions & their weights, 3=t(1) comparison or F(2)
    if sum(strcmp(log.reg_famtype, {'ContextOnset';'ContextMem'}))==1
        instruc.spec_context={
            'SR-SN'       {'SimR' 1; 'SimN' -1}       1;
            'SN-SR'       {'SimN' 1; 'SimR' -1}       1;
            'DR-DN'       {'DisR' 1; 'DisN' -1}       1;
            'DN-DR'       {'DisN' 1; 'DisR' -1}       1;
            'SR-DR'       {'SimR' 1; 'DisR' -1}       1;
            'DR-SR'       {'DisR' 1; 'SimR' -1}       1;
            'SN-DN'       {'SimN' 1; 'DisN' -1}       1;
            'DN-SN'       {'DisN' 1; 'SimN' -1}       1;
            'SimxVal'     {'SimR' 1; 'SimN' -1;'DisR' -1; 'DisN' 1}       2;
            };
        instruc.n_comparison_runs=1;
    elseif sum(strcmp(log.reg_famtype, {'ContextOnsetMem'}))==1
        instruc.spec_context={
            %                 'SR-SN'       {'SimR' 1; 'SimN' -1}       1;
            %                 'SN-SR'       {'SimN' 1; 'SimR' -1}       1;
            %                 'DR-DN'       {'DisR' 1; 'DisN' -1}       1;
            %                 'DN-DR'       {'DisN' 1; 'DisR' -1}       1;
            %                 'SR-DR'       {'SimR' 1; 'DisR' -1}       1;
            %                 'DR-SR'       {'DisR' 1; 'SimR' -1}       1;
            %                 'SN-DN'       {'SimN' 1; 'DisN' -1}       1;
            %                 'DN-SN'       {'DisN' 1; 'SimN' -1}       1;
            %                 'SimxVal'     {'SimR' 1; 'SimN' -1;'DisR' -1; 'DisN' 1}       2;
            };
        instruc.n_comparison_runs=1;
    end
end

for o2=1:1 % Execute prep  
    % (A) Change address of functional scans if necessary
    if request.change_fxnscan_address
        input('User requested that address of functional scans needs to be changed. Proceed?');
        for s=1:log.n_subjs
            try
                ws=load([where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep log.onsetsmodel '_' log.memtype ' Contrasted' filesep 'SPM.mat']);
                ws.SPM.old=ws.SPM;
                ws.preprocfolder_address=eval(request.change_fxnscan_address_newaddress); if isdir(ws.preprocfolder_address)==0; error(['Could not find location of preproc fxnals: ' ws.preprocfolder_address]); end
                
                % Identify fxnal scans (in order), re-write to SPM.xY.P & SPM.xY.VY.fname
                ws.fxnfiles=cell(size(ws.SPM.xY.P,1),2);
                for i=1: size(ws.SPM.xY.P,1) % identify
                    if isempty(strfind(ws.SPM.xY.P(i,:),log.func_prefix))==0
                        ws.fxnfiles{i,1}=[ws.preprocfolder_address ws.SPM.xY.P(i,strfind(ws.SPM.xY.P(i,:),log.func_prefix):end)];
                        ws.fxnfiles{i,2}=[ws.preprocfolder_address ws.SPM.xY.VY(i).fname(strfind(ws.SPM.xY.VY(i).fname,log.func_prefix):end)];
                        ws.SPM.xY.VY(i).fname=ws.fxnfiles{i,2};
                    else
                        error('Changing address of func scans - could not identify func scan name')
                    end
                end
                ws.SPM.xY.P=vertcat(char(ws.fxnfiles(:,1)));
                
                % Save
                disp(['Subject ' num2str(s) ' -  ' log.subjects{s} ' :  ' ws.SPM.xY.P(1,:)])
                SPM=ws.SPM;
                save([where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep log.onsetsmodel '_' log.memtype ' Contrasted' filesep 'SPM.mat'], 'SPM');
                ws=[];  SPM=[];
            catch
                errorlog{e,1}=['Couldn''t change address of fxn scans for subject  ' log.subjects{s}]; disp(errorlog{e,1}); e=e+1;
            end
        end
    end
    
    % (B) Implement omnibus F-test (all involved conditions)
    if request.contrast_Fomnibus
        disp('############# Implementing Omnibus contrasts ###############')
        disp(['fMRI model:   ' log.firstlevelmodel]); disp(instruc.Fomni); disp('Included regressors (adjusting data for effects of interest)'); disp(instruc.Fomni{:,2});
        input('Instructions OK for this fMRI model?   ');
        for f=1:size(instruc.Fomni,1)
            disp(['F omnibus for ' instruc.Fomni{f,1} ' ######'])
            
            % Identify target regressors and generate their weights
            wf.regnums=f_GetRegressorNums(p.SPM.xX.name',[instruc.Fomni{f,2} num2cell(instruc.Fomni{f,3}*ones(size(instruc.Fomni{f,2})))]);
            
            for s=1:log.n_subjs
                try
                    disp(['Subject ' num2str(s) ' - ' log.subjects{s}])
                    
                    % Compile weights
                    ws.s=load([where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep 'SPM.mat']);
                    ws.weights=zeros(size(wf.regnums,1), length(ws.s.SPM.xX.name));
                    for i=1:size(wf.regnums,1)
                        if length(wf.regnums{i,2})==1
                            ws.weights(i, wf.regnums{i,2})=1;
                        end
                    end
                    
                    % Batch in SPM
                    matlabbatch{1}.spm.stats.con.spmmat={[where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep 'SPM.mat']};
                    matlabbatch{1}.spm.stats.con.consess{1}.fcon.name=['Fomni_' instruc.Fomni{f,1}];
                    matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec = {ws.weights};
                    matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
                    matlabbatch{1}.spm.stats.con.delete = 0;
                    %
                    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                    matlabbatch=[]; ws=[];
                catch
                    errorlog{e,1}=['Couldn''t implement Fomni contrasts (' instruc.Fomni{f,1} ') for subject  ' log.subjects{s}]; disp(errorlog{e,1}); e=e+1;
                end
            end
            
            wf=[];
        end
    end
    
    % (C-1) [Context] Implement specific comparisons at the first-level (subject-specific contrasts)
    if request.contrast_SpecificComparisons
        disp('############# Implementing Specific contrast comparisons###############')
        
        for j=1:instruc.n_comparison_runs % Execute for both context onsets & context memory pmods if necessary
            instruc.spec_context_RegType=j;
            
            % --------------------------------------------------------------------------------
            % Format
            switch log.reg_famtype
                case 'ContextOnset'
                    instruc.spec_context(:,4)=num2cell(ones(size(instruc.spec_context,1),1)); 
                    instruc.spec_context_typename=log.reg_famprefix;
                case 'ContextMem'
                    instruc.spec_context(:,4)=num2cell(2*ones(size(instruc.spec_context,1),1)); 
                    for i=1:size(instruc.spec_context,1) % Alter names
                        instruc.spec_context{i,2}(:,1)=cellstr([char(instruc.spec_context{i,2}(:,1))    repmat(['xCMem_' log.memtype], size(instruc.spec_context{i,2},1),1)]);
                    end
                    instruc.spec_context_typename=log.reg_famprefix;
                otherwise
                    instruc.spec_context(:,4)=num2cell(ones(size(instruc.spec_context,1),1)); 
                    instruc.spec_context_typename=[];
            end
            disp(['Regressor type:   ' instruc.spec_context_typename '  --------------'])
            
            % Construct instructions (for contrast weights) for each comparisons
            for c=1:size(instruc.spec_context,1)
                
                % Identify regressors to be weighted for each comparison
                wc.req=[instruc.spec_context{c,2}(:,1) num2cell(instruc.spec_context_RegType*ones(size(instruc.spec_context{c,2},1),1))];
                wc.regnums=f_GetRegressorNums(p.SPM.xX.name, wc.req);
                for i=1:size(wc.regnums,1); if length(wc.regnums{i,2})~=1; error('Weighting contrasts for specific comparisons - more than 1 regressor found to weight'); end; end
                instruc.spec_context{c,2}(:,3)=wc.regnums(:,2);
                
                wc=[];
            end
            
            % Execute contrasts
            for s=1:log.n_subjs
                try
                    disp(['Subject ' num2str(s) ' - ' log.subjects{s} ' ----------------' ])
                    ws.s=load([where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep 'SPM.mat']);
                    %
                    matlabbatch{1}.spm.stats.con.spmmat={[where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep 'SPM.mat']};
                    matlabbatch{1}.spm.stats.con.delete = 0;
                    
                    for c=1:size(instruc.spec_context,1)
                        wc.weights=zeros(1,length(ws.s.SPM.xX.name));
                        for i=1:size(instruc.spec_context{c,2},1)
                            wc.weights(instruc.spec_context{c,2}{i,3})=instruc.spec_context{c,2}{i,2};
                        end
                        
                        switch instruc.spec_context{c,3}
                            case 1
                                matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = [instruc.spec_context_typename '_' instruc.spec_context{c,1}];
                                matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec=wc.weights;
                                matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
                            case 2
                                matlabbatch{1}.spm.stats.con.consess{c}.fcon.name = [instruc.spec_context_typename '_' instruc.spec_context{c,1}];
                                matlabbatch{1}.spm.stats.con.consess{c}.fcon.convec={wc.weights};
                                matlabbatch{1}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
                        end
                        wc=[];
                    end
                    
                    %
                    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                    matlabbatch=[]; ws=[];
                catch
                    errorlog{e,1}=['Couldn''t implement specific contrasts for subject  ' log.subjects{s}]; disp(errorlog{e,1}); e=e+1;
                end
            end
        end
    end
end

%% (3) Extract VOIs from all subjects


log.VoiImgs=cell(length(request.extract_FLvoi_ImgMask), log.n_subjs);  % Row = VOI, Col =Subject
if request.extract_firstlevelVOI==1 && request.VOI_individuals==0 % Group VOIs 
    disp('Group VOIs ####################################')
    
    % (A) Extract Group VOI on the basis of an image mask ####################
    for v=1:size(request.extract_FLvoi_ImgMask,1)
        disp(['VOI # ' num2str(v) ' - ' request.extract_FLvoi_ImgMask{v}])
        wv.found=0;
        
        % (1) Results image?
        wv.where=[where.exp_folder filesep '2 Second level results' log.AnalysisType filesep log.firstlevelmodel filesep log.secondlevelmodel filesep 'ROI' filesep 'PPI vois'];
        f=spm_select('List',wv.where,['^' request.extract_FLvoi_ImgMask{v} '.*.img$']);
        if isempty(f)==0
            wv.roi_ad=[wv.where filesep f]; wv.found=1;
        end
        
        % (2)  Cross-model results images (in Anatomical folder)
        if isempty(strfind(log.AnalysisType,'Ants'))==1; wv.where=[where.exp_folder filesep '3a Anatomical' filesep '3 VOI seeds' filesep 'All models'];
        else wv.where=[where.exp_folder filesep '3b Anatomical Ants' filesep '3 VOI seeds' filesep 'All models'];
        end
        f=spm_select('List', wv.where, [request.extract_FLvoi_ImgMask{v} '.*.img$']);
        if isempty(f)==0 && wv.found==0
            wv.roi_ad=[wv.where filesep f]; wv.found=1;
        end
        
        % (3) Anatomical ROI
        if isempty(strfind(log.AnalysisType,'Ants'))==1; wv.where=[where.exp_folder filesep '3a Anatomical' filesep '2 A priori ROIs' filesep 'All models'];
        else wv.where=[where.exp_folder filesep '3b Anatomical Ants' filesep '2 A priori ROIs' filesep 'All models'];
        end
        f=spm_select('List', wv.where, [request.extract_FLvoi_ImgMask{v} '.*.img$']);
        if isempty(f)==0 && wv.found==0
            wv.roi_ad=[wv.where filesep f]; wv.found=1;
        end
        
        if wv.found==0; error(['Could not find the ROI image file   (' request.extract_FLvoi_ImgMask{v} ')']); end
        log.VoiImgs(v, 1:end)=repmat(cellstr(wv.roi_ad), 1, log.n_subjs);
        wv=[];
    end
    disp('New procedure should work but not tested out yet!')
    
elseif request.extract_firstlevelVOI==1 && request.VOI_individuals==1 % Individual VOIs
    disp('Individual VOIs ####################################')
    

    % (A) Extract Individual VOI on the basis of an image mask ####################
    for v=1:size(request.extract_FLvoi_ImgMask,1)
        disp(['VOI # ' num2str(v) ' - ' request.extract_FLvoi_ImgMask{v}])
        
        % Images are assumed to be in FL subfolder 'VOI imgs'
        for s=1:log.n_subjs
            ws.wherevoi=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep  log.firstlevelmodel ' Contrasted' filesep 'VOI imgs' filesep];
            f=spm_select('List', ws.wherevoi, ['^' request.extract_FLvoi_ImgMask{v}]);
            if size(f,1)~=1; error(['No. of (individual subject) VOI images found ~=1   (' log.subjects{s} '  -  ' request.extract_FLvoi_ImgMask{v} ')']); end
            %
            log.VoiImgs{v, s}=[ws.wherevoi f(1,:)];
        end
    end
end
log.VoiSph=cell(length(request.extract_FLvoi_ImgMask), log.n_subjs);  % Row = VOI details (x y z mm), Col =Subject


%% (3) Extract VOIs from all subjects

if request.extract_firstlevelVOI
    disp('########### Extracting VOI ###############################')
    for s=1:log.n_subjs
%         try
            disp(['Subject ' num2str(s) ' (' log.subjects{s} ') ---------------'])
            ws.subfol=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep];
            ws.s=load([ws.subfol 'SPM.mat']);
            ws.contrastlist={ws.s.SPM.xCon(:).name}';
            
            % (A) Extract VOI on the basis of an image mask (Group-same or individual) ####################
            if isempty(request.extract_FLvoi_ImgMask)==0
                for v=1:size(request.extract_FLvoi_ImgMask,1)
                    disp(['VOI # ' num2str(v) ' - ' request.extract_FLvoi_ImgMask{v}])
                    matlabbatch{1}.spm.util.voi.spmmat = {[where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep 'SPM.mat']};
                    matlabbatch{1}.spm.util.voi.session = 1;
                    matlabbatch{1}.spm.util.voi.adjust = find(strcmp(ws.contrastlist,['Fomni_' log.reg_famtype]), 1, 'last'); % Omni-F contrast for THIS type of contrast only
                    wr.roi_ad=log.VoiImgs{v,s}; if isempty(wr.roi_ad)==1; error('VOI image not specified!'); end
                    disp(['VOI extracted on the basis of: ' wr.roi_ad]);
                    %
                    matlabbatch{1}.spm.util.voi.name = request.extract_FLvoi_ImgMask{v};
                    matlabbatch{1}.spm.util.voi.roi{1}.mask.image={[wr.roi_ad ',1']};
                    matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
                    matlabbatch{1}.spm.util.voi.expression = 'i1';
                    if request.VOI_maskwFLmask==1
                        matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {[where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep 'mask.img,1']};
                        matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
                        matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';
                    end
                    %
                    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                    matlabbatch=[]; wr=[];
                end
            end

            % (B-1) Extract Group VOI from peak + sphere of x radius ####################
            if isempty(request.extract_FLvoi_PeakSphere_Group)==0
                for v=1:size(request.extract_FLvoi_PeakSphere_Group,1)
                    disp(['VOI #' num2str(v) '  -  '  request.extract_FLvoi_PeakSphere_Group{v,1} ' (' num2str(request.extract_FLvoi_PeakSphere_Group{v,3}) ' mm sphere at ' num2str(request.extract_FLvoi_PeakSphere_Group{v,2}) ')'])
                    matlabbatch{1}.spm.util.voi.spmmat = {[where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep 'SPM.mat']};
                    matlabbatch{1}.spm.util.voi.session = 1;
                    matlabbatch{1}.spm.util.voi.adjust = find(strcmp(ws.contrastlist,['Fomni_' log.reg_famtype]), 1, 'last'); % Omni-F contrast for THIS type of contrast only
                    %
                    matlabbatch{1}.spm.util.voi.name = ['sph_' request.extract_flvoi_peaksphere_group{v,1}];
                    matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = request.extract_flvoi_peaksphere_group{v,2};
                    matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = request.extract_flvoi_peaksphere_group{v,3};
                    matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
                    matlabbatch{1}.spm.util.voi.expression = 'i1';
                    if request.voi_maskwflmask==1
                        matlabbatch{1}.spm.util.voi.roi{2}.mask.image={[where.data_brain filesep log.subjects{s} filesep '2 first level' log.analysistype filesep log.firstlevelmodel ' contrasted' filesep 'mask.img,1']};
                        matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
                        matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';
                    end
                    %
                    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                    matlabbatch=[]; wr=[];
                    
                    % Append to VOI mat file: Details of sphere (size & peak)
                    VOI_spherespex={};
                    if exist([ws.subfol 'VOI_sphere_details.mat'], 'file')~=0
                        wr.voispex=load([ws.subfol 'VOI_sphere_details.mat']);
                        VOI_spherespex=wr.voispex.VOI_spherespex;
                        wr=[];
                    end
                    vv=size(VOI_spherespex,1)+1;
                    VOI_spherespex{vv,1}=date;
                    VOI_spherespex{vv,2}=request.extract_FLvoi_PeakSphere_Group{v,1};
                    VOI_spherespex{vv,3}=[num2str(request.extract_FLvoi_PeakSphere_Group{v,3}) 'mm'];
                    VOI_spherespex{vv,4}=request.extract_FLvoi_PeakSphere_Group{v,2};
                    save([ws.subfol 'VOI_sphere_details.mat'], 'VOI_spherespex');
                end
            end 

            
            % (B-2) Extract Individual VOIs from peak + sphere of x radius (different for each subject) ####################
            if isempty(request.extract_FLvoi_PeakSphere_Individual)==0
                for v=1:size(request.extract_FLvoi_PeakSphere_Group,1)
%                     disp(['VOI #' num2str(v) '  -  '  request.extract_FLvoi_PeakSphere_Group{v,1} ' (' num2str(request.extract_FLvoi_PeakSphere_Group{v,3}) ' mm sphere at ' num2str(request.extract_FLvoi_PeakSphere_Group{v,2}) ')'])
%                     matlabbatch{1}.spm.util.voi.spmmat = {[where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep 'SPM.mat']};
%                     matlabbatch{1}.spm.util.voi.session = 1;
%                     matlabbatch{1}.spm.util.voi.adjust = find(strcmp(ws.contrastlist,['Fomni_' log.reg_famtype]), 1, 'last'); % Omni-F contrast for THIS type of contrast only
%                     %
%                     matlabbatch{1}.spm.util.voi.name = ['sph_' request.extract_flvoi_peaksphere_group{v,1}];
%                     matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = request.extract_flvoi_peaksphere_group{v,2};
%                     matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = request.extract_flvoi_peaksphere_group{v,3};
%                     matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
%                     matlabbatch{1}.spm.util.voi.expression = 'i1';
%                     if request.voi_maskwflmask==1
%                         matlabbatch{1}.spm.util.voi.roi{2}.mask.image={[where.data_brain filesep log.subjects{s} filesep '2 first level' log.analysistype filesep log.firstlevelmodel ' contrasted' filesep 'mask.img,1']};
%                         matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
%                         matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';
%                     end
%                     %
%                     spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
%                     matlabbatch=[]; wr=[];
%                     
%                     % Append to VOI mat file: Details of sphere (size & peak)
%                     VOI_spherespex={};
%                     if exist([ws.subfol 'VOI_sphere_details.mat'], 'file')~=0
%                         wr.voispex=load([ws.subfol 'VOI_sphere_details.mat']);
%                         VOI_spherespex=wr.voispex.VOI_spherespex;
%                         wr=[];
%                     end
%                     vv=size(VOI_spherespex,1)+1;
%                     VOI_spherespex{vv,1}=date;
%                     VOI_spherespex{vv,2}=request.extract_FLvoi_PeakSphere_Group{v,1};
%                     VOI_spherespex{vv,3}=[num2str(request.extract_FLvoi_PeakSphere_Group{v,3}) 'mm'];
%                     VOI_spherespex{vv,4}=request.extract_FLvoi_PeakSphere_Group{v,2};
%                     save([ws.subfol 'VOI_sphere_details.mat'], 'VOI_spherespex');
                end
            end 
            
%         catch
%             errorlog{e,1}=['Couldn''t extract VOI for subject  ' log.subjects{s}]; disp(errorlog{e,1}); e=e+1;
%         end
    end
end

%% END

disp('======================================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
disp(['GLM model: ' log.AnalysisType '      ' log.firstlevelmodel '     ' log.secondlevelmodel]); disp(' ')
disp('Analysis completed:'); disp(request);
disp(' '); disp('Errors:'); disp(errorlog)
disp('=======================================================')

diary off
try % Notify researcher
    f_sendemail('kurzlich', ['Analysis batchscript is complete (' mfilename ')'], ' ',1);
end
