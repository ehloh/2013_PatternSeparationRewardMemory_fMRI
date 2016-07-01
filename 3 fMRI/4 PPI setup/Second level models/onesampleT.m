function [ batch] = onesampleT(where, log, conimages)
% [ batch] = onesampleT(where, log, conimages)
% One sample t-test. Specify name of contrasts (conimages{1}).

% Execute to debug: conimages{1}='PPI Pos'

if length(conimages)~=1; error('Wrong no. of contrast images specified for this type of 2nd level model'); end

%% (1) Specify 2nd level model 

% Set up batch
disp('Running one sample ttest ------------')
matlabbatch{1}.spm.stats.factorial_design.dir ={where.secondlevelres};
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});

% Identify & Load con image
p=load([where.data_brain filesep log.subjects{1} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep 'PPI ' log.PPImodelname filesep 'SPM.mat']);
cons=cell(size(p.SPM.xCon,2),1); for i=1:size(p.SPM.xCon,2); cons{i}=p.SPM.xCon(i).name; end
if sum(strcmp(cons,conimages{1}))~=1; error(['Error: Could not find target contrast image (' conimages{1} ')']);end
log.targetcon=p.SPM.xCon(find(strcmp(cons,conimages{1}))).Vcon.fname;
%    
for s=1:log.n_subjs
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{s,1}=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.AnalysisType filesep log.firstlevelmodel ' Contrasted' filesep 'PPI ' log.PPImodelname filesep log.targetcon ',1'];
end

%
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
batch{1}=matlabbatch; matlabbatch=[]; 
log.ingoingcontrasts=conimages;
save([where.secondlevelres filesep 'details_2ndlevel.mat'], 'log'); 

%% (2) Estimate 2nd level model

disp('Estimating model ------------')
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[where.secondlevelres 'SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
batch{2}=matlabbatch; matlabbatch=[]; 

%% (3) Add positive and negative contrasts

disp('Running contrasts on 2nd level model ------------')
matlabbatch{1}.spm.stats.con.spmmat =  {[where.secondlevelres 'SPM.mat']};
matlabbatch{1}.spm.stats.con.delete = 0;
%
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'PPI Pos';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'PPI Neg';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = -1;
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
%
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
batch{3}=matlabbatch; matlabbatch=[]; 

end

