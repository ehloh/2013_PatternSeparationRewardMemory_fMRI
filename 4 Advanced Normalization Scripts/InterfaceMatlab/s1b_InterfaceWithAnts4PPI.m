% Interface with ANTs, transferring files between SPM data folder and Ants (Neurodebian shared) folder
clear all;close all hidden; clc

where.where='J:\1 fMRI analysis'; where.where='D:\1 [Context-Memory] fMRI Data\1 fMRI analysis'; where.exp_folder='D:\1 [Context-Memory] fMRI Data'; where.data_beh=[where.where '\2 Behavioural data'];where.data_brain=[where.exp_folder filesep '1 MRI data']; where.antsfolder='D:\host\2_Data_AdjustCons';

% Steps?
log.specificsubjects={}; 
request.TransformedSeedImgs_ToPC=1;
request.ConsToAnts=0;
request.AdjustedConsToPC=0;

% PPI details
% request.seed='sph_SNVTA_R1_cmSRvSN';
% request.seed='sph_SNVTA_R2_cmSRvSN';
request.seed='sph_HPC_aL_cmSRvSN';
request.psych='DisN';
request.ppi_conditiontype=1; % 1=Context onsets, 2=Context pmod, 3= ContextOnsetMem

%
request.WhichAnalysisThread='s4FullCardiacWithDerivAnts';
request.WhichFLmodel='m_ci3_ContextItemNomotor_Hit';
request.FL_shortname='s4FullCardiacWithDeriv_ci3Hit';
request.AdjustConMethod='Landmarks5';


for o1=1:1 % General settings and specifications
    
    % Load subjects
    addpath(where.where)
    w.onsets=load([where.exp_folder filesep '1 MRI data' filesep 'Onsetslog_' request.WhichFLmodel]); % This filter is only applied first to catch subjects with bad onsets files, where this has yet to be marked in excel sheet
    log.w=load([where.exp_folder filesep '1 MRI data' filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog,[],vertcat({'Subject' 'onsets'}, w.onsets.subjectsok),'onsets');
%     w.specsubs=log.specificsubjects;
    if request.TransformedSeedImgs_ToPC==0;
        [w.s w.s1 log.koshertable]=xlsread([where.where filesep 'i8_subjectdataok_secondlevel.xlsx'],'PPIsingle'); 
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable, [request.WhichFLmodel '_' request.AdjustConMethod ' ' request.seed]);
        log.specificsubjects=log.subjects;
    else
        [w.s w.s1 log.koshertable]=xlsread([where.where filesep 'i8_subjectdataok_secondlevel.xlsx']); 
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable, [request.WhichFLmodel '_' request.AdjustConMethod]); % Select for this FL model
        log.specificsubjects=log.subjects;
    end
    
    % PPI details
    switch request.ppi_conditiontype
        case 1; request.ppi_condition_prefix='co';
        otherwise; error('not specified yet')
    end
    log.ppiname=['PPI ' request.ppi_condition_prefix 'Fam_' request.seed '_psy_' request.psych];
    log.pc_FLfol=['2 First level ' request.WhichAnalysisThread  filesep request.WhichFLmodel '_' request.AdjustConMethod ' Contrasted' filesep log.ppiname filesep];
    
    % Ants
    log.AntsTypes={'Basic';'Landmarks';'Template';'Landmarks2';'Landmarks4';'Landmarks5';'Landmarks6'};
    if isempty(strfind(request.FL_shortname, request.WhichAnalysisThread(1:end-4)))==1; error('Analysis Thread not in FL short name (for neurodebian folder)'); end
    if isempty(strfind(request.FL_shortname, request.WhichFLmodel(3:5)))==1;     error('Requested FL folder not in FL short name (for neurodebian folder)'); end
    if isempty(strfind(request.FL_shortname, request.WhichFLmodel(max(strfind(request.WhichFLmodel, '_'))+1:end)))==1;     error('Requested FL memtype not in FL short name (for neurodebian folder)'); end
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' ');     
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location: ' where.data_brain]); disp(' ');  disp(['Data location (behaviour): ' where.data_beh])
    disp(' ');  input('Hit Enter to start      ')
    disp('=======================================================')

end

%% (1)

if request.TransformedSeedImgs_ToPC
    
    
    request.Imgs2Bring={    % Col 1=Ants name,  Col 2=Name of Voi
%         'sph_SNVTA_R1_cmSRvSN2'         'sph_SNVTA_R1_cmSRvSN';       % Spheres (as images)
%         'sph_SNVTA_R2_cmSRvSN2'         'sph_SNVTA_R2_cmSRvSN';
%         'sph_HPC_aL_cmSRvSN4'               'sph_HPC_aL_cmSRvSN';


        %
%         'LM5cm_SRvSN_001_HPC_aLpeak'        'HPC_aLpeak_cmSRvSN';       % Results peaks
%         'LM5cm_SRvSN_001_SNVTA_Rpeak'      'SNVTA_Rpeak_cmSRvSN';
%         'LM5co_SRvSN_001_SNVTA_Lpeak'       'SNVTA_Lpeak_coSRvSN';
%         'LM5co_SRvO_001_SNVTA_midpeak'       'SNVTA_Mpeak_coSRvO';
%         'LM5cm_SNvSR_001_HPC_pLpeak'            'HPC_pLpeak_cmSNvSR'; 
%
'LM5_SNVTA_R_cmSRvSN'        'SNVTA_R_cmSRvSN';
'LM5_HPC_aL_cmSRvSN'        'HPC_aL_cmSRvSN';
'LM5_SNVTA_L_coSRvSN'        'SNVTA_L_coSRvSN';
        };
    
    
    % -----------------------------------------------------------------
    disp('VOIs to move over (file name, new voi name):'); disp(request.Imgs2Bring); input('OK?  ');
    disp('Bring VOIs (in subject space) to PC ################')
    for s=1:log.n_subjs
        disp(log.subjects{s})
        ws.from=[where.antsfolder filesep log.subjects{s} filesep request.AdjustConMethod filesep 'InverseTransform' filesep];
        ws.to=[where.data_brain filesep log.subjects{s} filesep '2 First level ' request.WhichAnalysisThread filesep request.WhichFLmodel '_' request.AdjustConMethod ' Contrasted' filesep 'VOI imgs' filesep];
        
        % Transfer all vois for this subject
        for v=1:size(request.Imgs2Bring,1)
            copyfile([ws.from log.subjects{s} '_' request.Imgs2Bring{v,1} '.nii'], [ws.to request.Imgs2Bring{v,2} '.nii'])
        end
        
        % Save details
        if exist([ws.to 'voi_details.mat'],'file')==0; vois=[];
        else;  ws.v=load([ws.to 'voi_details.mat']); vois=ws.v.vois;
        end
        vois=[vois; [{date} {request.FL_shortname} {request.AdjustConMethod} {request.Imgs2Bring}]];
        save([ws.to 'voi_details.mat'],'vois')
        
        
        ws=[]; vois=[];
    end
    
end

%%  (2) Unadjusted cons --> Ants 

if request.ConsToAnts
    request.saveoriginal=1;
    disp('Bring UnAdjusted contrasts from PC to Ants ################')
    for s=1:log.n_subjs
        disp(log.subjects{s})
        ws.ants_fol=[where.antsfolder filesep log.subjects{s} filesep request.AdjustConMethod filesep];
        ws.ants_ppifol=[ws.ants_fol request.FL_shortname filesep log.ppiname(5:end) filesep]; mkdir(ws.ants_ppifol);
        ws.pc_ppifol=[where.data_brain filesep log.subjects{s} filesep '2 First level ' request.WhichAnalysisThread filesep request.WhichFLmodel '_' request.AdjustConMethod ' Contrasted' filesep log.ppiname filesep];
        ws.pc_ppifol_cleancons=[ws.pc_ppifol 'PreAnts cons' filesep];
        
        
        % Save original contrasts imgs
        if request.saveoriginal
            if isdir([ws.pc_ppifol 'PreAnts cons'])==0
                mkdir([ws.pc_ppifol 'PreAnts cons'])
                movefile([ws.pc_ppifol 'con_0001.hdr'], [ws.pc_ppifol_cleancons 'con_0001.hdr'])
                movefile([ws.pc_ppifol 'con_0001.img'], [ws.pc_ppifol_cleancons 'con_0001.img'])
                movefile([ws.pc_ppifol 'con_0002.hdr'], [ws.pc_ppifol_cleancons 'con_0002.hdr'])
                movefile([ws.pc_ppifol 'con_0002.img'], [ws.pc_ppifol_cleancons 'con_0002.img'])
                % 
            else
                disp('PreAnts folder found. Not saving orignal cons, assumed already saved. OK?')
                if s==1; input('OK?   '); end
            end
        end
        
        % Convert
        v=spm_vol([ws.pc_ppifol_cleancons 'con_0001.img']) ;ima=spm_read_vols(v);
        v.fname=[ws.pc_ppifol_cleancons 'con_0001.nii']; spm_write_vol(v,ima);
%         v=spm_vol([ws.pc_ppifol_cleancons 'con_0002.img']) ;ima=spm_read_vols(v);
%         v.fname=[ws.pc_ppifol_cleancons 'con_0002.nii']; spm_write_vol(v,ima);
        
        % Copy over
        copyfile([ws.pc_ppifol_cleancons  'con_0001.nii'],  [ws.ants_ppifol 'con_0001.nii'])
%         copyfile([ws.pc_ppifol_cleancons  'con_0002.nii'],  [ws.ants_ppifol 'con_0002.nii'])
        
        %
        ws=[];
    end
end

%% (2) Adjust cons, from Ants --> SPM first-level models on PC

if request.AdjustedConsToPC
    
    disp('Bring Adjusted contrasts from Ants to PC ################')
    request.ConvertFirst=1;
    
    for s=1:log.n_subjs
        disp(log.subjects{s})
        ws.ants_fol=[where.antsfolder filesep log.subjects{s} filesep request.AdjustConMethod filesep];
        ws.ants_ppifol=[ws.ants_fol request.FL_shortname filesep log.ppiname(5:end) filesep]; 
        ws.pc_ppifol=[where.data_brain filesep log.subjects{s} filesep '2 First level ' request.WhichAnalysisThread filesep request.WhichFLmodel '_' request.AdjustConMethod ' Contrasted' filesep log.ppiname filesep];
        ws.pc_ppifol_cleancons=[ws.pc_ppifol 'PreAnts cons' filesep];
        
        % Convert
            if request.ConvertFirst
                matlabbatch{1}.spm.util.imcalc.input = {[ws.ants_ppifol 'r_con_0001.nii,1']};
                matlabbatch{1}.spm.util.imcalc.output = 'r_con_0001.img';
                matlabbatch{1}.spm.util.imcalc.outdir = {ws.ants_ppifol};
                matlabbatch{1}.spm.util.imcalc.expression = 'i1';
                matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
                matlabbatch{1}.spm.util.imcalc.options.mask = 0;
                matlabbatch{1}.spm.util.imcalc.options.interp = 1;
                matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
%                 matlabbatch{2}.spm.util.imcalc.input = {[ws.ants_ppifol 'r_con_0002.nii,1']};
%                 matlabbatch{2}.spm.util.imcalc.output = 'r_con_0002.img';
%                 matlabbatch{2}.spm.util.imcalc.outdir = {ws.ants_ppifol};
%                 matlabbatch{2}.spm.util.imcalc.expression = 'i1';
%                 matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
%                 matlabbatch{2}.spm.util.imcalc.options.mask = 0;
%                 matlabbatch{2}.spm.util.imcalc.options.interp = 1;
%                 matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
                spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                matlabbatch=[];
            end
            
            % Copy to FL folder on PC
            copyfile([ws.ants_ppifol 'r_con_0001.hdr'], [ws.pc_ppifol 'con_0001.hdr']);
            copyfile([ws.ants_ppifol 'r_con_0001.img'], [ws.pc_ppifol 'con_0001.img']);
%             copyfile([ws.ants_ppifol 'r_con_0002.hdr'], [ws.pc_ppifol 'con_0002.hdr']);
%             copyfile([ws.ants_ppifol 'r_con_0002.img'], [ws.pc_ppifol 'con_0002.img']);
            
            ws=[];
    end
end





