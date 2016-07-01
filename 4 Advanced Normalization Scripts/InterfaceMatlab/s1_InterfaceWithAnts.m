% Interface with ANTs, transferring files between SPM data folder and Ants (Neurodebian shared) folder
clear all;close all hidden; clc
where.where='D:\Dropbox\SCRIPPS\1 ContextMem fMRI'; where.exp_folder='D:\1 [Context-Memory] fMRI Data'; where.data_beh=[where.where '\2 Behavioural data'];where.data_brain=[where.exp_folder filesep '1 MRI data']; where.antsfolder='D:\host\Context_study';

% Steps?
request.TransferT1sForTemplateConstruc=0;
request.TransferStrxConsForConAdjustment=0;
request.TransferAdjustedConToPC=1;

% Which subjects
log.specificsubjects={};
% log.specificsubjects={'p01_CW'}; 
log.specificsubjects={
%     'p01_CW';'p03_EA';'p04_JL'; 'p07_LH';'p08_AM';'p09_CN';'p10_AB';
%     'p11_SS';'p12_IL';  'p14_SJ';'p16_TH';'p17_RA';
% 'p18_KB';'p19_CN';'p20_JB';'p21_SH';'p22_EK';'p23_IS';'p24_LL';'p25_BS';'p26_MC';'p27_EW';'p28_CN';'p29_ET'
    };
for o1=1:1 % General settings and specifications
    log.AntsTypes={'Basic';'Landmarks';'Template';'Landmarks2';'Landmarks4';'Landmarks5';'Landmarks6'};
    % Load subjects
    addpath(where.where)
    log.w=load([where.data_brain filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    log.allsubjects=log.subjects; 
    
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

%% Step 1: Transfer structural (T1w) scans for template constcruction

if request.TransferT1sForTemplateConstruc
    where.to_T1w=[where.antsfolder '\1_Data_CreateHPCTemplate'];  % In Host Neurodebian folder, where to?
%     where.to_T1w=[where.antsfolder '\3a_Data_CreateTemplatePartial']; disp('Partial volume template!!')
    
    where.T1wFol='1 Preprocessed\MPM\Antst1w';
    disp('Transferring T1ws for Template construction -----------------------------------')
    for s=1: log.n_subjs
        disp([log.subjects{s} ' ------------------']);
        ws.subfol_pc=[where.data_brain filesep log.subjects{s} filesep ];
        
        
        % Mean functiona (ubf), is coregistered to all functionals (non-normalized,
        % prefix ubf). Convert to nii format. - used to create a partial
        % volume template!        
%         f=spm_select('List', [ws.subfol_pc '1 Preprocessed\Func_b1\Preproc_b1'], '^meanubf.*.img'); if size(f,1)~=1; error('More than 1 mean functional found!'); end
%         v=spm_vol([ws.subfol_pc '1 Preprocessed\Func_b1\Preproc_b1\' f(1,:)]);
%         ima=spm_read_vols(v);
%         v.fname=[ws.subfol_pc where.T1wFol filesep f(1,1:length(f(1,:))-4) '.nii'];
%         spm_write_vol(v,ima);
        
        try
            f=spm_select('List', [ws.subfol_pc filesep where.T1wFol], '^bsMQ.*_T1w.nii'); % Bias corrected T1w structurals
%             f=spm_select('List', [ws.subfol_pc filesep where.T1wFol], '^meanubf.*.nii'); disp('partial volume template!')
            if size(f,1)~=1; error('No. found files ~=1'); end
            copyfile([ws.subfol_pc where.T1wFol filesep f], [where.to_T1w filesep log.subjects{s} '_' f]);
        catch
            disp('Did not')
        end
        
        ws=[];
        
    end
end

%% Step 2: For FL Contrast adjustment, bring structural scans & FL contrasts over to Ants (specific structurals that were used in functional preprocessing!)
% These structurals must be co-registered to the functionals that were went into the first-level GLMs. To verify, check-reg the coregistered
%     (but not normalized) functionals with the structurals :). Note also: Structurals used in preprocessing were NOT bias-corrected, 
%     unlike the structurals that went into constructing the template in Ants. Assume ok

if request.TransferStrxConsForConAdjustment
    request.TransferStrx=0;
    request.TransferCons=1;
    request.AdjustConMethod='Landmarks5'; 
    request.WhichAnalysisThread='s4FullCardiacWithDerivAnts'; request.WhichAnalysisThread_short='s4FullCardiacWithDeriv'; 
    requests.WhichContrasts=[]; % 1:16; % Blank to process all    
    
    % Which First-level Analysis thread + Model? 
%     request.WhichFLmodel='m_ci3_ContextItemNomotor_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_ci3Hit'];
%     request.WhichFLmodel='m_ci2_ContextonlyItem_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_ci2Hit'];
%     request.WhichFLmodel='m_i4_ItemContextevent_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_i4Hit'];
%     request.WhichFLmodel='m_i5_ItemContextonly_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_i5Hit'];
    request.WhichFLmodel='m_c4_ContextallItempresonly_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_c4Hit'];
%     request.WhichFLmodel='m_i3_ItemContextpresonly_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_i3Hit'];
%     request.WhichFLmodel='m_i6_ItemMempmodContextpresonly_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_i6Hit'];
%     request.WhichFLmodel='m_ci13_ContextStickItemNomotor_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_ci13Hit'];
  
    request.WhichFLmodel='m_c7_ContexteventItempresonly_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_c7Hit'];


    for o1=1:1% General setup
        
        % Folders
        where.ants_AdjustCons=[where.antsfolder '\2_Data_AdjustCons'];  % In Host Neurodebian folder, where to? All subjects' T1s & cons go here (no subject folders)
        where.from_PreprocStruc='1 Preprocessed\MPM';
        where.from_ConImgs=['2 First level ' request.WhichAnalysisThread filesep request.WhichFLmodel ' Contrasted'];
        
        % Checks
        if isempty(strfind(request.WhichAnalysisThread, 'Ants'))==1; error(['Analysis thread doesn''t include Ants! (' request.WhichAnalysisThread ')']); end
%         if isempty(strfind(request.WhichFLmodel, request.FLmodShortname(strfind(request.FLmodShortname, '_c')+1:strfind(request.FLmodShortname, '_c')+3)))==1; 
        if strcmp(request.WhichFLmodel(3:4), request.FLmodShortname(strfind(request.FLmodShortname, '_i')+1  :  strfind(request.FLmodShortname, '_i')+2)) + strcmp(request.WhichFLmodel(3:4), request.FLmodShortname(strfind(request.FLmodShortname, '_c')+1  :  strfind(request.FLmodShortname, '_c')+2)) ~=1, error ('Short FL model name doesn''t match long FL model name!');  end 
        if sum(strcmp(log.AntsTypes, request.AdjustConMethod))~=1; error(['Invalid contrast-adjustment method selected: '  request.AdjustConMethod]); end
        
        % Redo subject selections: Only valid subjects for requested model!
        [w.s w.s1 log.koshertable]=xlsread([where.where filesep 'i8_subjectdataok_secondlevel.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable, request.WhichFLmodel);

        disp(['N subjects for this model: ' num2str(log.n_subjs) ' ------------------------------']);
        disp('Transferring coregistered-structurals and first-level contrast images over ---------------------------------------')
    end
    disp('Requesting contrasts from:'); disp(['       Analysis Thread:      ' request.WhichAnalysisThread]);  disp(['       First-Level model:   '  request.WhichFLmodel]);  disp(['       First-Level model:   '  request.FLmodShortname ' (short name)']);     disp(['       Specific constrasts requested: ' num2str(requests.WhichContrasts)]); disp(['       Adjust Con Method: '  request.AdjustConMethod]); input('OK?  ')  
    for s=1: log.n_subjs
        disp([log.subjects{s} ' ------------------']);
        ws.subfol_pc=[where.data_brain filesep log.subjects{s} filesep];
        ws.subfol_ants=[where.ants_AdjustCons filesep log.subjects{s} filesep];
        ws.submodfol_ants=[ws.subfol_ants request.AdjustConMethod filesep request.FLmodShortname filesep];
        if isdir(ws.submodfol_ants)==0; mkdir(ws.submodfol_ants); end
        if isdir([ws.subfol_ants request.AdjustConMethod filesep 'InverseTransform'])==0; mkdir([ws.subfol_ants request.AdjustConMethod filesep 'InverseTransform']); end
        
        
        % Structurals
        if request.TransferStrx
            f=spm_select('List', [ws.subfol_pc where.from_PreprocStruc], '^sMQ.*_T1w.nii'); % Bias corrected T1w structurals
            if size(f,1)~=1; error('No. found files ~=1'); end
            copyfile([ws.subfol_pc where.from_PreprocStruc filesep f], [ws.subfol_ants filesep log.subjects{s} '_T1w_coreg.nii' ]); % rename structurals
        end
        
        % First-level contrasts - Convert to nii, and transfer over
        if request.TransferCons
            if isempty(requests.WhichContrasts);  f=spm_select('List', [ws.subfol_pc where.from_ConImgs], '^con_.*.img');
            else f=num2str(10000+requests.WhichContrasts'); f=[repmat('con_', length(requests.WhichContrasts), 1) f(:,2:end) repmat('.img', length(requests.WhichContrasts), 1)];
            end
            if isempty(f)==1; error('No valid contrasts to move over!'); end
            for i=1:size(f,1) % Convert to nii
                v=spm_vol([ws.subfol_pc where.from_ConImgs filesep f(i,:)]);
                ima=spm_read_vols(v);
                v.fname=[ws.subfol_pc where.from_ConImgs filesep f(i,1:length(f(i,:))-4) '.nii'];
                spm_write_vol(v,ima);
            end
            for i=1:size(f,1) % Move over!
                copyfile([ws.subfol_pc where.from_ConImgs filesep f(i,1:length(f(i,:))-4) '.nii'] ,[ws.submodfol_ants filesep log.subjects{s} '_' f(i,1:length(f(i,:))-4) '.nii']) % .img
            end
        end
        
        ws=[];
        
    end
    
end

% Undo model-specific subjects selections 
log.subjects=log.allsubjects; log.n_subjs=length(log.subjects);

%% Step 4: Bring adjusted first-level contrast images back from Ants to PC

if request.TransferAdjustedConToPC
    request.ConvertFirst=1;
    
    % Which First-level Analysis thread + Model?
    request.WhichAnalysisThread='s4FullCardiacWithDerivAnts'; request.WhichAnalysisThread_short='s4FullCardiacWithDeriv'; 
    %
%     request.WhichFLmodel='m_ci3_ContextItemNomotor_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_ci3Hit'];
%     request.WhichFLmodel='m_ci2_ContextonlyItem_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_ci2Hit'];
%     request.WhichFLmodel='m_i4_ItemContextevent_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_i4Hit'];
%     request.WhichFLmodel='m_i5_ItemContextonly_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_i5Hit'];
    request.WhichFLmodel='m_c4_ContextallItempresonly_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_c4Hit'];
%         request.WhichFLmodel='m_i3_ItemContextpresonly_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_i3Hit'];
%     request.WhichFLmodel='m_i6_ItemMempmodContextpresonly_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_i6Hit'];
%     request.WhichFLmodel='m_ci13_ContextStickItemNomotor_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_ci13Hit'];
  request.WhichFLmodel='m_c7_ContexteventItempresonly_Hit';  request.FLmodShortname=[request.WhichAnalysisThread_short '_c7Hit'];

    
    
    request.AdjustConMethod='Landmarks5';
% request.AdjustConMethod='Template';
    requests.WhichContrasts=[]; % Blank to process all (not coded up yet!!)
    
    for o1=1:1% General setup
        
        % Checks
        if isempty(strfind(request.WhichAnalysisThread, 'Ants'))==1; error(['Analysis thread doesn''t include Ants! (' request.WhichAnalysisThread ')']); end
        if strcmp(request.WhichFLmodel(3:4), request.FLmodShortname(strfind(request.FLmodShortname, '_i')+1  :  strfind(request.FLmodShortname, '_i')+2)) + strcmp(request.WhichFLmodel(3:4), request.FLmodShortname(strfind(request.FLmodShortname, '_c')+1  :  strfind(request.FLmodShortname, '_c')+2)) ~=1, error ('Short FL model name doesn''t match long FL model name!');  end 
        if sum(strcmp(log.AntsTypes, request.AdjustConMethod))~=1; error(['Invalid contrast-adjustment method selected: '  request.AdjustConMethod]); end
        disp('Requesting contrasts from:'); disp(['       Analysis Thread:      ' request.WhichAnalysisThread]);  disp(['       First-Level model:   '  request.WhichFLmodel]);  disp(['       First-Level model:   '  request.FLmodShortname ' (short name)']);     disp(['       Specific constrasts requested: ' num2str(requests.WhichContrasts)]); disp(['       Adjust Con Method: '  request.AdjustConMethod]); input('OK?  ')
        
        % Folders
        where.ants_AdjustCons=[where.antsfolder '\2_Data_AdjustCons'];  % In Host Neurodebian folder, where to?
        
        % Redo subject selections: Only valid subjects for requested model!
        [w.s w.s1 log.koshertable]=xlsread([where.where filesep 'i8_subjectdataok_secondlevel.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable, request.WhichFLmodel);
        disp(['N subjects for this model: ' num2str(log.n_subjs) ' ------------------------------']);
        disp('Transferring Ants-adjusted first-level contrast images back to PC ---------------------------------------')
    end
    
    for s=1: log.n_subjs
        disp([log.subjects{s} ' ------------------']);
        ws.subfol_ants=[where.ants_AdjustCons filesep log.subjects{s} filesep];
        ws.submodfol_ants=[where.ants_AdjustCons filesep log.subjects{s} filesep request.AdjustConMethod filesep request.FLmodShortname filesep];
        ws.submodfol_pcorig=[where.data_brain filesep log.subjects{s} filesep '2 First level ' request.WhichAnalysisThread filesep request.WhichFLmodel ' Contrasted' filesep];
        ws.submodfol_pc=[where.data_brain filesep log.subjects{s} filesep '2 First level ' request.WhichAnalysisThread filesep request.WhichFLmodel '_' request.AdjustConMethod ' Contrasted' filesep];
        if isdir(ws.submodfol_pc)==0; copyfile(ws.submodfol_pcorig, ws.submodfol_pc); end
        
        
        % Transfer adjusted contrast images (converted back to .hdr/.img)
        if isempty(requests.WhichContrasts)==1;
            f=spm_select('List', ws.submodfol_ants,  ['^r_' log.subjects{s} '_con_.*']); % Ants-adjusted FL contrasts
        else
            f=[];
            for i=1:length(requests.WhichContrasts)
                ff=spm_select('List', ws.submodfol_ants,  ['^r_' log.subjects{s} '_con_.*0' num2str(requests.WhichContrasts(i)) '.nii']);
                f=[f; ff];
            end
        end
        if isempty(f)==1; error('No valid contrasts to move over!');  else disp(['Moving ' num2str(size(f,1)) ' contrast images']); end
        for i=1:size(f,1)
            disp(['     ' ws.submodfol_ants f(i,:)]);
            % % Don't use this code - faulty!
            % v=spm_vol([ws.submodfol_ants f(i,:)]);
            % ima=spm_read_vols(v);
            % v.fname=[ws.submodfol_pc f(i,10:end-4) '.img'];
            % spm_write_vol(v,ima);
            if request.ConvertFirst
                matlabbatch{1}.spm.util.imcalc.input = {[ws.submodfol_ants f(i,:) ',1']};
                matlabbatch{1}.spm.util.imcalc.output = [f(i,1:end-4) '.img'];
                matlabbatch{1}.spm.util.imcalc.outdir = {[ws.submodfol_ants ]};
                matlabbatch{1}.spm.util.imcalc.expression = 'i1';
                matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
                matlabbatch{1}.spm.util.imcalc.options.mask = 0;
                matlabbatch{1}.spm.util.imcalc.options.interp = 1;
                matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
                spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                matlabbatch=[];
            end
            
            % Copy to FL folder on PC
            copyfile([ws.submodfol_ants  f(i,1:end-4) '.img'], [ws.submodfol_pc f(i,10:end-4) '.img']);
            copyfile([ws.submodfol_ants  f(i,1:end-4) '.hdr'], [ws.submodfol_pc f(i,10:end-4) '.hdr']);
        end
        
        % Mark in spm.mat variable, details of this contrast's adjustments 
        %       Ongoing history - col 1=date, col 2=Adjustment type, col 3=Adjusted contrast images
        ws.s=load([ws.submodfol_pc 'SPM.mat']);
        if isfield(ws.s.SPM,'ConAdjustLog')==1;  k=size(ws.s.SPM.ConAdjustLog,1)+1;
        else; ws.s.SPM.ConAdjustLog=cell(1,3); k=1;
        end
        ws.s.SPM.ConAdjustLog{k,1}=date;
        ws.s.SPM.ConAdjustLog{k,2}=request.AdjustConMethod;
        ws.s.SPM.ConAdjustLog{k,3}=f;
        SPM=[]; SPM=ws.s.SPM; save([ws.submodfol_pc 'SPM.mat'], 'SPM');
        
        %
        ws=[];
    
    end
    
 
    try % Notify researcher
        f_sendemail('kurzlich', strcat('Analysis script is complete (ContextMem AntsInterface)'), ' ',1);
    end    
    
end

% Undo model-specific subjects selections 
log.subjects=log.allsubjects; log.n_subjs=length(log.subjects);

