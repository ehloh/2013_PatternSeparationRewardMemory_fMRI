% Preprocessing: realign & unwarp, coregister, segment, normalize, smooth
clear all;close all hidden; clc

where.where='I:\1 fMRI analysis'; where.data_brain='C:\Users\eloh\Desktop\1 [Context-Memory] fMRI Data\1 MRI data';  % where.data_beh=[where.where filesep '2 Behavioural data'];
% where.where='/Volumes/PENNYDISK/1 fMRI analysis'; where.data_brain='/Volumes/SANDISK/1 CONTEXT Brain data';% where.data_beh=[where.where filesep '2 Behavioural data'];

% Requested analysis
process.t1=1;
process.mt=1;
%
process.smoothingsize=6; % Smoothing size?
log.specificsubjects={}; % BLANK to process all subjects

for o1=1:1 % General settings and specifications
    
    % Subjects
    addpath(where.where)
    log.w=load([where.data_brain filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    cd(where.data_brain); cd ..; where.data_folder=pwd; cd(where.where);
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
    
%     where.spm='D:\My Documents\MATLAB\spm8'; spm fmri
end

%% Process average T1

% Normalize
matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50;  78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';
for s=1:log.n_subjs
    ws.where=[where.data_brain filesep log.subjects{s} filesep '1 Preprocessed' filesep 'MPM' filesep];
    f=spm_select('List', ws.where, '^sM.*._T1w_seg_sn.mat'); % Parameter file
    matlabbatch{1}.spm.spatial.normalise.write.subj(s).matname ={[ws.where f]};
    f=spm_select('List', ws.where, '^sM.*._T1w.nii'); % File to write
    matlabbatch{1}.spm.spatial.normalise.write.subj(s).resample = {[ws.where f ',1']};
end

spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];
            


disp('############## Imcalc-ing ######################')


% Imcalc
matlabbatch{1}.spm.util.imcalc.output =  ['AverageT1_n' num2str(log.n_subjs) '.img'];
matlabbatch{1}.spm.util.imcalc.outdir = {[where.data_folder filesep '3 Checks']};
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
matlabbatch{1}.spm.util.imcalc.input =wb.t1s;
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
    f=spm_select('List', wb.where, '^wsM.*._T1w.nii$');
    wb.t1s{s} =[wb.where f ',1'];
end


spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];

%% END

disp('=======================================================')
w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp(' ')
disp('Analysis completed:')
disp(process)
disp(['No. of subjects: ' num2str(log.n_subjs)])
disp(' ')
disp(errorlog)
disp(' ')
disp('=======================================================')
