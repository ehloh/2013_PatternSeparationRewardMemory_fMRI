% Plot movement, Organize, Checkreg
clear all; close all hidden; clc

where.where='I:\1 fMRI analysis';  where.data_brain='D:\1 [Context-Memory] fMRI Data\1 MRI data'; where.data_beh=[where.where filesep '2 Behavioural data'];
% where.where='/Volumes/PENNYDISK/1 fMRI analysis'; where.data_brain='/Volumes/SANDISK/1 CONTEXT Brain data';% where.data_beh=[where.where filesep '2 Behavioural data'];

% Requested analysis
log.specificsubjects={}; % BLANK to process all subjects
request.Setup1stLevelFolder=1;
request.OrganizePreprocessingFolder=0;
request.PlotMovement=0;
request.CheckReg=0;

% Type of model I'm setting up
log.AnalysisType=' s6FullCardiacWithDeriv';
log.preprocfuncs_prefix='s6wu';
log.CardiacType='_10cardiac';

for o1=1:1 % General settings and specifications    
    
    % Subjects
    addpath(where.where)
    log.w=load([where.data_brain filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
    % Things that don't change much
    request.OrganizePreprocessingFolder_zipunused=0;
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('====================================')
    w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' ')
    disp('Requested:')
    disp(request)
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0
        disp('   Subset of subjects only:')
        disp(log.specificsubjects)
    end
    disp(' ')
    disp(['Data location: ' where.data_brain])
    disp(['FIRST LEVEL SUFFIX: ' log.AnalysisType])
    disp(' ')
    input('Hit Enter to start      ')
    disp('====================================')
    
end

%% (1) Plot movement regressors

if request.PlotMovement==1
    disp('----------------- Collecting movement parameters ----------------- ')
    movement=cell(log.n_subjs,3);
    for s=1:log.n_subjs
        for b=1:2
            f=spm_select('List', [where.data_brain filesep log.subjects{s} filesep  '1 Preprocessed' filesep 'Func_b' num2str(b)], '.txt$');
            if isempty(f)==1
                movement{s,b}=[];
            else
                movement{s,b}=load([where.data_brain filesep log.subjects{s} filesep  '1 Preprocessed' filesep 'Func_b' num2str(b) filesep f]);
            end
        end
        movement{s,3}=log.subjects{s};
    end
    % Plot
    figure('Name', 'Movement parameters', 'Position', [400 100 1200 1000])
    for s=1:log.n_subjs
        subplot(ceil(log.n_subjs/5), round(log.n_subjs/ceil(log.n_subjs/5)), s);
        plot(movement{s,1}); hold on; plot(movement{s,2}); 
        axis tight
    end
end

%% (2) Organize 

if request.OrganizePreprocessingFolder==1
    disp('--------------- Organizing scans --------------------')
    for s=1:log.n_subjs
        wb.where=[where.data_brain filesep log.subjects{s} filesep  '1 Preprocessed' filesep];
        for b=1:2
            eval('java.io.File([wb.where ''Func_b'' num2str(b)]).renameTo(java.io.File([wb.where ''Preproc_b'' num2str(b)]));') % RENAME
            if isdir([wb.where 'Preproc_b' num2str(b)])==0; input('Error: Could not rename folder, in organizing preprocessing files'); end
            if isdir([wb.where 'Func_b' num2str(b)])==0; mkdir([wb.where 'Func_b' num2str(b)]); end
            % Identify target files (to be kept)
            f=spm_select('List',[wb.where 'Preproc_b' num2str(b)], ['^' log.preprocfuncs_prefix '*']);
            ff=spm_select('List',[wb.where 'Preproc_b' num2str(b)], 'txt$'); f=vertcat(f,ff);
            for i=1:size(f,1) % Move
                movefile([wb.where 'Preproc_b' num2str(b) filesep f(i,:)],[wb.where 'Func_b' num2str(b) filesep f(i,:)]);
            end
            movefile([wb.where 'Preproc_b' num2str(b)],[wb.where 'Func_b' num2str(b) filesep 'Preproc_b' num2str(b)]);
            % Zip if requested
            if request.OrganizePreprocessingFolder_zipunused==1
                zip([wb.where 'Func_b' num2str(b) filesep 'Preproc_b' num2str(b)],[wb.where 'Func_b' num2str(b) filesep 'Preproc_b' num2str(b)]);
            end
        end
        wb=[];
    end
end

%% (3) CheckReg 

if request.CheckReg==1
    disp('--------------- Collecting scans for CheckReg --------------------')
    w.nScansPerFigure=10;
    %
    images=cell(log.n_subjs*2,1); i=1;
    for s=1:log.n_subjs
        for b=1:2
            wb.where=[where.data_brain filesep log.subjects{s} filesep '1 Preprocessed' filesep 'Func_b' num2str(b) filesep];
            % Identify scans
            f=spm_select('List', wb.where, ['^' log.preprocfuncs_prefix '.*img$']);
            if isempty(f)==1
                disp(['Error: No available scans for ' log.subjects{s} '  block ' num2str(b) ' ---'])
            else
                images{i,1}=[wb.where f(randi(size(f,1)),:) ',1'];
                images{i,2}=log.subjects{s};  i=i+1;
            end
        end
    end
    
    % Collate into figures
    f=1; i=1; figs=cell(ceil(size(images,1)/w.nScansPerFigure),1);
    for d=1:size(images,1)
        figs{f}{i,1}=char(images{d,1});
        images{d,3}=f;
        if i==w.nScansPerFigure
            i=1; f=f+1;
        else
            i=i+1; 
        end
    end
    
    % Display / Instructions for display
    if w.nScansPerFigure>=size(images,1)
        matlabbatch{1}.spm.util.checkreg.data=figs{f};
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    else
        disp(' --------------------- INSTRUCTIONS FOR CHECKREG -------------------------------------')
        disp(['Multiple figures to display (' num2str(size(figs,1)) ' figures)'])
        disp('Note: To change no. of images per display, specify in script (at beginning of module)')
        disp('See variable ''images'' to identify subject scans (Col 3)')
        disp(' ')
        disp('To display each batch, specify value of ''f'' and execute following command:')
        disp(' ')
        disp('f = # ; eval(checkregcommand) ')
        disp(' ')        
        disp(' -----------------------------------------------------------------------------------------------')
        spm_jobman('initcfg');  checkregcommand='matlabbatch{1}.spm.util.checkreg.data=figs{f}; spm_jobman(''run'' , matlabbatch);';
    end
end

%% (4) Prep for First level

if request.Setup1stLevelFolder==1
    disp('--------------- Preparing folder for 1st level --------------------')
    for s=1:log.n_subjs
        disp(['Subject  ' num2str(s) '  (' log.subjects{s} ') '])
        wb.where=[where.data_brain filesep log.subjects{s} filesep];
        wb.wherefrom=[wb.where '1 Preprocessed' filesep];
        wb.whereto=[wb.where '2 First level' log.AnalysisType filesep];
        if isdir(wb.whereto)==0; mkdir(wb.whereto); end;
        if isdir([wb.whereto 'Preproc functionals'])==0; mkdir([wb.whereto 'Preproc functionals']); end
        try % Spike + Physio (Standardize name)
            copyfile([wb.wherefrom log.subjects{s} '_reg_physiomovement' log.CardiacType '.txt'],[wb.whereto log.subjects{s} '_reg_physiomovement.txt']);
        end
        f=spm_select('List', [wb.wherefrom 'Func_b1'],log.preprocfuncs_prefix);
        for i=1:size(f,1)
            copyfile([wb.wherefrom 'Func_b1' filesep f(i,:)],[wb.whereto 'Preproc functionals' filesep f(i,:)]);
        end
        f=spm_select('List', [wb.wherefrom 'Func_b2'],log.preprocfuncs_prefix);
        for i=1:size(f,1)
            copyfile([wb.wherefrom  'Func_b2' filesep f(i,:)],[wb.whereto 'Preproc functionals' filesep f(i,:)]);
        end
    end
end

%%  

disp(' ###################################################')
disp(' DONE ' )
if request.CheckReg==1
    disp('CHECKREG: See command window (earlier command) for display instructions')
end
disp(' ###################################################')



