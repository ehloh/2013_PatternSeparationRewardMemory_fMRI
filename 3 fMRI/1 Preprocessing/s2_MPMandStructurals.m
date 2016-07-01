% Process all MPMs
clear all;close all hidden; clc

where.where='I:\1 fMRI analysis'; where.data_brain='D:\1 [Context-Memory] fMRI Data\1 MRI data';  % where.data_beh=[where.where filesep '2 Behavioural data'];
% where.where='/Volumes/PENNYDISK/1 fMRI analysis'; where.data_brain='/Volumes/SANDISK/1 CONTEXT Brain data';% where.data_beh=[where.where filesep '2 Behavioural data'];

% Requested analysis
process.mpm=0;
process.mpm_1folder=0;
process.averagestructurals=1;
process.applysubjectselectionformodel={'m_ci3_ContextItemNomotor_Hit'};   % {'m_ci3_ContextItemNomotor_Hit'}; % Apply selection criteria as per 2nd level?
log.specificsubjects={}; % BLANK to process all subjects

for o1=1:1 % General settings and specifications    
    
    % Subjects
    addpath(where.where)
    log.w=load([where.data_brain filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    if ~isempty(process.applysubjectselectionformodel)
        [w.s w.s1 log.koshertable]=xlsread([where.where filesep 'i8_subjectdataok_secondlevel.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable, char(process.applysubjectselectionformodel)) ;
    end
    
    
    % Log
    cd(where.data_brain); cd ..; where.data_folder=pwd; cd(where.where);
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' (' date ')'])
    errorlog=cell(1,1); e=1;
 
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp('Requested analysis:'); disp(process)
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location: ' where.data_brain])
    disp(' ');  input('Hit Enter to start      ')
    disp('=======================================================')
    
    
end

%% (1) Process MPMs

if process.mpm==1
    disp(' ########## PROCESS MPM ##########')
    for i=1:length(log.subjects)
        try
            disp(['Subject ' num2str(i) ' -------------------'])
            wb.where=[where.data_brain filesep log.subjects{i} filesep '1 Preprocessed'];
            cd(wb.where)
            wb.nb1=log.datalog{i+1,9}; % IDs
            wb.nb0=wb.nb1+2; % B0: 1st run (2 runs)
            wb.nMT=wb.nb1+3; % MTw: 1st run
            wb.nPD=wb.nb1+5; % PDw: 1st run
            wb.nT1=wb.nb1+7; % T1w: 1st run
            % Field map
            wb.where_b0=[wb.where filesep 'MPM_b0' filesep]; % B0
            f=spm_select('List', wb.where_b0,  '^sMQ.*\.img$');
            wb.Q=[repmat(wb.where_b0,size(f,1),1) f repmat(',1',size(f,1),1)];
            wb.where_b1=[wb.where filesep 'MPM_b1' filesep]; % B1
            f=spm_select('List',wb.where_b1,'^sMQ.*\.img$');
            wb.P=[repmat(wb.where_b1,size(f,1),1) f repmat(',1',size(f,1),1)];
            B1map_v2(wb.P,wb.Q,1192)
            % Anatomicals
            f=spm_select('List',[wb.where filesep 'MPM_MTw' filesep], ['00' num2str(wb.nMT) '-.*\.img$']); % MTw
            wb.P_mt=[repmat([wb.where filesep 'MPM_MTw' filesep],size(f,1),1) f repmat(',1',size(f,1),1)];
            f=spm_select('List',[wb.where filesep 'MPM_PDw' filesep], ['00' num2str(wb.nPD) '-.*\.img$']); % PDw
            wb.P_pd=[repmat([wb.where filesep 'MPM_PDw' filesep],size(f,1),1) f repmat(',1',size(f,1),1)];
            f=spm_select('List',[wb.where filesep 'MPM_T1w' filesep], ['00' num2str(wb.nT1) '-.*\.img$']); % T1w
            wb.P_t1=[repmat([wb.where filesep 'MPM_T1w' filesep],size(f,1),1) f repmat(',1',size(f,1),1)];        
            f=spm_select('List', wb.where_b1, 'smuB1map_');  % b1 + sensitivity
            wb.P_trans=[repmat(wb.where_b1,size(f,1),1) f repmat([',1'],size(f,1),1)];        
            wb.P_receiv=[];
            MT_analysis_altered(wb.P_mt,wb.P_pd,wb.P_t1,wb.P_trans,wb.P_receiv) % Altered MT_analysis script to accept inputs
            %
            wb=[];
        catch
            errorlog{e,1}=['Failed: Process MPM --- ' log.datalog{s+1,1}];
            e=e+1;
        end
    end
end

%% (2) Re-organize MPMs into single folder

if process.mpm_1folder
    disp(' ########## ORGANIZE MPMs ##########')
    for s=1:length(log.subjects)
        try
            disp(['Subject ' num2str(s) ' -------------------'])
            wb.where=[where.data_brain filesep log.subjects{s} filesep '1 Preprocessed'];
            wb.whereto=[wb.where filesep 'MPM' filesep];
            wb.wherefrom=[wb.where filesep 'MPM_MTw' filesep];
            cd(wb.where)
            wb.mpmfolders=dir('MPM*');
            mkdir(wb.whereto) 
            wb.file='T1w'; % T1w
            wb.f=spm_select('List', wb.wherefrom, ['_' wb.file '.']);
            movefile([wb.wherefrom wb.f], [wb.whereto wb.f]);
            wb.file='PDw'; % PDw
            wb.f=spm_select('List', wb.wherefrom, ['_' wb.file '.']);
            movefile([wb.wherefrom wb.f], [wb.whereto wb.f]);
            wb.file='MTw'; % MTw
            wb.f=spm_select('List', wb.wherefrom, ['_' wb.file '.']);
            movefile([wb.wherefrom wb.f], [wb.whereto wb.f]);
            for i=1:size(wb.mpmfolders,1) % Move other folders
                movefile([wb.where filesep wb.mpmfolders(i).name], [wb.whereto wb.mpmfolders(i).name]);
            end
            wb=[];
        catch
            errorlog{e,1}=['Failed:MPM sorting --- ' log.datalog{s+1,1}];
            e=e+1;
        end
    end
end

%% (3) Calculate average structurals (T1w & MTw)

if process.averagestructurals==1
    disp('############ Constructing average structurals ###############')
    average={'T1'; 'MT'}; 
    
    for o1=1:2
        disp(['[Group ' average{o1} '] NORMALIZE-----------'])
        
        % Normalize
        matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
        matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50;  78 76 85];
        matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 1;
        matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';
        for s=1:log.n_subjs
            ws.where=[where.data_brain filesep log.subjects{s} filesep '1 Preprocessed' filesep 'MPM' filesep];
            f=spm_select('List', ws.where, '^sM.*._T1w_seg_sn.mat'); % Parameter file (always from T1)
            matlabbatch{1}.spm.spatial.normalise.write.subj(s).matname ={[ws.where f]};
            f=spm_select('List', ws.where, ['^sM.*._' average{o1} 'w.nii']); % File to write
            matlabbatch{1}.spm.spatial.normalise.write.subj(s).resample = {[ws.where f ',1']};
        end
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
        matlabbatch=[];
        
        % Imcalc
        disp(['[Group ' average{o1} '] IMCALC-ING GROUP SCAN----------'])
        matlabbatch{1}.spm.util.imcalc.output =  ['Average' average{o1} '_n' num2str(log.n_subjs) '.img'];
        matlabbatch{1}.spm.util.imcalc.outdir = {[where.data_folder filesep '3 Checks']};
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        matlabbatch{1}.spm.util.imcalc.expression='(';
        for s=1:length(log.subjects) % Generate expression
            matlabbatch{1}.spm.util.imcalc.expression=[matlabbatch{1}.spm.util.imcalc.expression 'i' num2str(s)];
            if s<log.n_subjs
                matlabbatch{1}.spm.util.imcalc.expression=[matlabbatch{1}.spm.util.imcalc.expression '+'];
            else
                matlabbatch{1}.spm.util.imcalc.expression=[matlabbatch{1}.spm.util.imcalc.expression ')/' num2str(log.n_subjs)];
            end
        end
        for s=1:length(log.subjects) % Collect scans
            wb.where=[where.data_brain filesep log.subjects{s} filesep '1 Preprocessed' filesep 'MPM' filesep];
            f=spm_select('List', wb.where, ['^wsM.*._' average{o1} 'w.nii$']);
            wb.t1s{s} =[wb.where f ',1'];
        end
        matlabbatch{1}.spm.util.imcalc.input =wb.t1s;
        %
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
        matlabbatch=[];
    end
end

%% END

disp('====================================')
w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp(' ')
disp('Analysis completed:')
disp(process)
disp(['No. of subjects: ' num2str(log.n_subjs)])
disp(' ')
disp('Error log?')
disp(errorlog)
disp(' ')
disp('====================================')

diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s2_MPMandStructurals)'), ' ',1);
catch
end
