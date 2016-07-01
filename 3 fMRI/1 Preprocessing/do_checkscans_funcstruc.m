clear all; clc
where.where='H:\1 fMRI analysis';
where.data='C:\Users\eloh\Desktop\0 [Context-Memory] fMRI Data\MRI';
load([where.data filesep 'datalog_allsubjects.mat'])

% Which subject?
s=1;

%%

prefix='ubf';
sub=datalog{s+1,1}; disp(sub)
scanID=datalog{s+1,3};
b1=num2str(datalog{s+1,6});
b2=num2str(datalog{s+1,7});


% Identify T1
strucscan=spm_select('List', [where.data filesep sub filesep '1 Preprocessed\MPM'], '^sMQ.*T1w.nii');


%
matlabbatch{1}.spm.util.checkreg.data = {
                                [where.data filesep sub filesep '1 Preprocessed\Func_b1\' prefix scanID '-000' b1 '-00002-000002-01.img,1']
                                [where.data filesep sub filesep '1 Preprocessed\Func_b2\' prefix scanID '-000' b2 '-00181-000181-01.img,1']
                                [where.data filesep sub filesep '1 Preprocessed\MPM\' strucscan ',1']
                                         };
                                     
% spm fmri
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
disp(sub)