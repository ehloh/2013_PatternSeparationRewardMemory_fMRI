% Preprocessing: realign & unwarp, coregister, segment, normalize, smooth
clear all;close all hidden; clc

where.where='I:\1 fMRI analysis'; where.exp_folder='D:\1 [Context-Memory] fMRI Data'; where.data_beh=[where.where '\2 Behavioural data'];where.data_brain=[where.exp_folder filesep '1 MRI data'];
% where.where='/Volumes/PENNYDISK/1 fMRI analysis'; where.exp_folder='/Volumes/SANDISK/1 CONTEXT Brain data'; where.data_beh='/Volumes/PENNYDISK/1 fMRI analysis/2 Behavioural data';

% Requested analysis
process.realignunwarp=0;
process.coregister=0;
process.segment=0;
process.normalize=0;
process.smooth=1;
%
process.smoothingsize=4; % Smoothing size?
log.specificsubjects={'p01_CW'}; % BLANK to process all subjects

for o1=1:1 % General settings and specifications
    
    % Subjects
    addpath(where.where)
    log.w=load([where.data_brain filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
    % Log
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
    
    where.spm='D:\My Documents\MATLAB\spm8'; spm fmri
end

%% Step 1: Realign & unwarp

if process.realignunwarp==1
    disp(' ########## (1) REALIGN & UNWARP: Create VDM files ##########')
    default.fieldmap=[where.spm filesep 'toolbox' filesep 'FieldMap' filesep 'pm_defaults_Trio_al_96.m'];
    for s=1:length(log.subjects)
        try
            disp(['Subject ' num2str(s)   '   (' log.subjects{s} ')  --------------- '])
            wb.where=[where.data_brain filesep log.subjects{s} filesep '1 Preprocessed' filesep];
            wb.whereFM=[wb.where 'Func_Fieldmap' filesep];
            % Specify
            wb.subj=[log.datalog{s+1,3} '.']; 
            wb.n_phase=log.datalog{s+1, 5}+1;
            wb.n_mag=log.datalog{s+1, 5};
            wb.n_epi1=log.datalog{s+1, 6};
            wb.n_epi2=log.datalog{s+1, 7};
            % Choose files
            f=spm_select('List', wb.whereFM, ['00' num2str(wb.n_phase) '-00001.*\.img$']); % Phase
            matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.phase=cellstr([wb.whereFM f ',1']);
            f=spm_select('List', wb.whereFM, ['00' num2str(wb.n_mag) '-00001-000001-01.img']); % Mag
            matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.magnitude=cellstr([wb.whereFM f ',1']);
            f=spm_select('List',[wb.where 'Func_b1'], '^bfM.*img$'); % EPI
            matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session(1).epi=cellstr([wb.where 'Func_b1' filesep f(1,:) ',1']);
            f=spm_select('List', [wb.where 'Func_b2'], '^bfM.*img$');
            matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session(2).epi=cellstr([wb.where 'Func_b2' filesep f(1,:) ',1']);
            matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsfile = {default.fieldmap};
            matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchvdm = 1;
            matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.sessname = 'session';
            matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.writeunwarped = 0;
            matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.anat = '';
            matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchanat = 0;
            % Run !
            spm_jobman('run',matlabbatch);
            %
            matlabbatch=[];
            wb=[];
        catch
            errorlog{e,1}=['Failed: Realign & Unwarp: Create VDM file  --- ' log.datalog{s+1,1}];
            e=e+1;
        end
    end
    
    % Execute unwarping ------------------
    disp(' ########## (1) REALIGN & UNWARP: Execute realigning & unwarping ##########')
    for o1 =1:1 % Realign & unwarp settings 
        settings.realignunwarp.eoptions.quality = 0.9;
        settings.realignunwarp.eoptions.sep = 4;
        settings.realignunwarp.eoptions.fwhm = 5;
        settings.realignunwarp.eoptions.rtm = 0;
        settings.realignunwarp.eoptions.einterp = 2;
        settings.realignunwarp.eoptions.ewrap = [0 0 0];
        settings.realignunwarp.eoptions.weight = '';
        settings.realignunwarp.uweoptions.basfcn = [12 12];
        settings.realignunwarp.uweoptions.regorder = 1;
        settings.realignunwarp.uweoptions.lambda = 100000;
        settings.realignunwarp.uweoptions.jm = 0;
        settings.realignunwarp.uweoptions.fot = [4 5];
        settings.realignunwarp.uweoptions.sot = [];
        settings.realignunwarp.uweoptions.uwfwhm = 4;
        settings.realignunwarp.uweoptions.rem = 1;
        settings.realignunwarp.uweoptions.noi = 5;
        settings.realignunwarp.uweoptions.expround = 'Average';
        settings.realignunwarp.uwroptions.uwwhich = [2 1];
        settings.realignunwarp.uwroptions.rinterp = 4;
        settings.realignunwarp.uwroptions.wrap = [0 0 0];
        settings.realignunwarp.uwroptions.mask = 1;
        settings.realignunwarp.uwroptions.uwroptions.prefix = 'u';
    end
    for s=1:length(log.subjects)
        try
            disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------------- '])
            matlabbatch{1}.spm.spatial.realignunwarp=settings.realignunwarp;
            wb.where=[where.data_brain filesep log.datalog{s+1,1} filesep '1 Preprocessed' filesep];
            wb.whereFM=[wb.where 'Func_Fieldmap' filesep];
            wb.FMrun=log.datalog{s+1,5}+1;
            % VDM files (1 per run)
            f=spm_select('List', wb.whereFM, '^vdm5_.*.session1.img');
            matlabbatch{1}.spm.spatial.realignunwarp.data(1).pmscan = {[wb.whereFM f ',1']};
            f=spm_select('List', wb.whereFM, '^vdm5_.*.session2.img');
            matlabbatch{1}.spm.spatial.realignunwarp.data(2).pmscan = {[wb.whereFM f ',1']};
            % Select EPIs
            f=spm_select('List', [wb.where 'Func_b1'], '^bfM.*.img'); 
            for i=1:size(f,1)
                wb.r1{i,1}=[wb.where 'Func_b1' filesep f(i,:) ',1'];
            end
            matlabbatch{1}.spm.spatial.realignunwarp.data(1).scans =wb.r1;
            f=spm_select('List', [wb.where 'Func_b2'], '^bfM.*.img');
            for i=1:size(f,1)
                wb.r2{i,1}=[wb.where 'Func_b2' filesep f(i,:) ',1'];
            end
            matlabbatch{1}.spm.spatial.realignunwarp.data(2).scans =wb.r2;            
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            %
            matlabbatch=[];
            wb=[];
        catch
            errorlog{e,1}=['Failed: Realign & Unwarp  --- ' log.datalog{s+1,1}];
            e=e+1;
        end
    end
end

%% Step 2: Coregister (no prefix)

if process.coregister==1
    disp(' ############### (2) COREGISTRATION ############ ##########')
    for o2=1:1 % Settings for Coregistration
        settings.coreg.other = {''};
        settings.coreg.eoptions.cost_fun = 'nmi';
        settings.coreg.eoptions.sep = [4 2];
        settings.coreg.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        settings.coreg.eoptions.fwhm = [7 7];
    end
    for s=1:length(log.subjects)
        try
            disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------------- '])
            matlabbatch{1}.spm.spatial.coreg.estimate=settings.coreg;  % Specifications 
            wb.where=[where.data_brain filesep log.datalog{s+1,1} filesep '1 Preprocessed' filesep];
            wb.where_run1=[wb.where filesep 'Func_b1' filesep];
            wb.where_run2=[wb.where filesep 'Func_b2' filesep];
            % Choose files
            f=spm_select('List', [wb.where 'MPM'], '^sMQ*.*T1w.nii'); % Reference: Structural
            matlabbatch{1}.spm.spatial.coreg.estimate.ref={[wb.where 'MPM' filesep f ',1']};
            f=spm_select('List', wb.where_run1, '^meanubfMQ*.*.img'); % Source image: meanubfM file (or, 1st volume from 1st run) 
            matlabbatch{1}.spm.spatial.coreg.estimate.source={[wb.where_run1 f ',1']};
            f  = spm_select('List', [wb.where 'Func_b1'], '^ubfMQ*.*img'); wb.b1size=size(f,1); % Functional scans
            for i=1:size(f,1)
                wb.func{i,1}=[wb.where 'Func_b1' filesep f(i,:) ',1'];
            end
            f  = spm_select('List', [wb.where 'Func_b2'], '^ubfMQ*.*img');
            for i=1:size(f,1)
                wb.func{i+wb.b1size,1}=[wb.where 'Func_b2' filesep f(i,:) ',1'];
            end
            matlabbatch{1}.spm.spatial.coreg.estimate.other = wb.func;
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            %
            matlabbatch=[];
            wb=[];
        catch
            errorlog{e,1}=['Failed: Coregister & Reslice --- ' log.datalog{s+1,1}];
            e=e+1;
        end
    end
end

%% Step 3: Segment (no prefix)

if process.segment==1 
    disp(' ############### (3) Normalization: Segmentation  ############ ##########')
    for o1=1:1 % Settings for Segmentation (performed only on T1)
        settings.segment.output.GM = [0 0 1];
        settings.segment.output.WM = [0 0 1];
        settings.segment.output.CSF = [0 0 0];
        settings.segment.output.biascor = 1;
        settings.segment.output.cleanup = 0;
        settings.segment.opts.tpm = {[where.spm filesep 'tpm' filesep 'grey.nii'];[where.spm filesep 'tpm' filesep 'white.nii'];[where.spm filesep 'tpm' filesep 'csf.nii']};
        settings.segment.opts.ngaus = [2;2;2; 4];
        settings.segment.opts.regtype = 'mni';
        settings.segment.opts.warpreg = 1;
        settings.segment.opts.warpco = 25;
        settings.segment.opts.biasreg = 0.0001;
        settings.segment.opts.biasfwhm = 60;
        settings.segment.opts.samp = 3;
        settings.segment.opts.msk = {''};
    end
    for s=1:length(log.subjects)
        try 
            disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------------- '])
            matlabbatch{1}.spm.spatial.preproc=settings.segment;
            f=spm_select('List',  [where.data_brain filesep log.datalog{s+1,1} filesep '1 Preprocessed' filesep 'MPM' ], '^sMQ*.*_T1w.nii');
            matlabbatch{1}.spm.spatial.preproc.data = {[where.data_brain filesep log.datalog{s+1,1} filesep '1 Preprocessed' filesep 'MPM' filesep f ',1' ]};
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            %
            matlabbatch=[];
            wb=[];
        catch
            errorlog{e,1}=['Failed: Segmentation of T1 --- ' log.datalog{s+1,1}];
            e=e+1;
        end
    end
end

%% Step 4: Normalization (prefix 'w')

if process.normalize==1  
    disp(' ############### (3) Normalization: Execute normalization ############ ##########')
    for o1=1:1 % Settings for Normalization 
        settings.normalization.preserve = 0;
        settings.normalization.bb = [-78 -112 -50;78 76 85];
        settings.normalization.vox = [2 2 2];
        settings.normalization.interp = 1;
        settings.normalization.wrap = [0 0 0];
        settings.normalization.prefix = 'w';
    end
    for s=1:length(log.subjects)
        try 
            disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------------- '])
            matlabbatch{1}.spm.spatial.normalise.write.roptions=settings.normalization;
            wb.where=[where.data_brain filesep log.datalog{s+1,1} filesep '1 Preprocessed' filesep];
            f=spm_select('List', [wb.where 'MPM'], '^sMQ*.*_T1w_seg_sn.mat'); % Parameter file
            matlabbatch{1}.spm.spatial.normalise.write.subj.matname = {[wb.where 'MPM' filesep f]}; 
            % Select images to normalize (EPIs)
            wb.volnum=1; 
            wb.run=1; % 1st run 
            wb.wherethisrun=[wb.where 'Func_b' num2str(wb.run) filesep];
            f=spm_select('List', wb.wherethisrun, '^ubfM*.*img');
            wb.nvols=size(f,1);
            for i=1:wb.nvols 
                wb.allvols{wb.volnum,:}=[wb.wherethisrun f(i,:) ',1'];
                wb.volnum=wb.volnum+1;
            end
            wb.run=2; % 2nd run
            wb.wherethisrun=[wb.where 'Func_b' num2str(wb.run) filesep];
            f=spm_select('List', wb.wherethisrun, '^ubfM*.*img');
            wb.nvols=size(f,1);
            for i=1:wb.nvols
                wb.allvols{wb.volnum,:}=[wb.wherethisrun f(i,:) ',1'];
                wb.volnum=wb.volnum+1;
            end
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample =wb.allvols;
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            %
            matlabbatch=[];
            wb=[];
        catch
            errorlog{e,1}=['Failed: Normalization  --- ' log.datalog{s+1,1}];
            e=e+1;
        end
    end
end

%% Step 5: Smoothing
% Always saved in different folders

if process.smooth==1  
    disp([' ############### (4) SMOOTHING (size: ' num2str(process.smoothingsize) 'x' num2str(process.smoothingsize) 'x' num2str(process.smoothingsize) ')   ############ ##########'])
    for o1=1:1 % Settings for Smoothing
        settings.smooth.fwhm = [process.smoothingsize process.smoothingsize process.smoothingsize]; % Smoothing size
        settings.smooth.dtype = 0;
        settings.smooth.im = 0;
%         settings.smooth.prefix='s';
        settings.smooth.prefix = ['s' num2str(process.smoothingsize)];
    end
    for s=1:length(log.subjects)
%         try 
            disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------------- '])
            matlabbatch{1}.spm.spatial.smooth=settings.smooth;
            wb.where=[where.data_brain filesep log.datalog{s+1,1} filesep '1 Preprocessed' filesep ];
            % Select images to smooth (EPIs)
            wb.volnum=1; 
            wb.run=1; % 1st run 
            wb.wherethisrun=[wb.where 'Func_b' num2str(wb.run) filesep 'Preproc_b' num2str(wb.run) filesep];
            f=spm_select('List', wb.wherethisrun, '^ubfM*.*img');
%             f=spm_select('List', wb.wherethisrun, '^wubfM*.*img');
            wb.nvols=size(f,1);
            for i=1:wb.nvols 
                wb.allvols{wb.volnum,:}=[wb.wherethisrun f(i,:) ',1'];
                wb.volnum=wb.volnum+1;
            end
            wb.run=2; % 2nd run
            wb.wherethisrun=[wb.where 'Func_b' num2str(wb.run) filesep 'Preproc_b' num2str(wb.run) filesep];
            f=spm_select('List', wb.wherethisrun, '^ubfM*.*img');
%             f=spm_select('List', wb.wherethisrun, '^wubfM*.*img');
            wb.nvols=size(f,1);
            for i=1:wb.nvols
                wb.allvols{wb.volnum,:}=[wb.wherethisrun f(i,:) ',1'];
                wb.volnum=wb.volnum+1;
            end
            matlabbatch{1}.spm.spatial.smooth.data =wb.allvols;
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            %
            matlabbatch=[];
            wb=[];
%         catch
%             errorlog{e,1}=['Failed: Smoothing (size: ' num2str(process.smoothingsize) ' --- ' log.datalog{s+1,1}];
%             e=e+1;
%         end
    end
end

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

diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s3_Preprocessing)'), ' ',1);
catch
end