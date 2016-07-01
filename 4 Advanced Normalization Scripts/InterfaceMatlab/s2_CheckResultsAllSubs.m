allsub={'p01_CW';
    'p02_MK';
    'p03_EA';
    'p04_JL';
    'p06_YL';
    'p07_LH';
    'p08_AM';
    'p09_CN';
    'p10_AB';
    'p11_SS';
    'p12_IL';
    'p13_CB';
    'p14_SJ';
    'p16_TH';
    'p17_RA';
    'p18_KB';
    'p19_CN';
    'p20_JB';
    'p21_SH';
    'p22_EK';
    'p23_IS';
    'p24_LL';
    'p25_BS';
    'p26_MC';
    'p27_EW';
    'p28_CN';
    'p29_ET';};

sub=allsub{3};

%%

% Subject
% matlabbatch{1}.spm.stats.results.spmmat = {['D:\1 [Context-Memory] fMRI Data\1 MRI data\' sub '\2 First level s4WithDerivAnts\m_ci3_ContextItemNomotor_Hit Contrasted\SPM.mat']};


% Group 2nd level results
% matlabbatch{1}.spm.stats.results.spmmat = {['D:\1 [Context-Memory] fMRI Data\2 Second level results s4WithDerivAnts\m_ci3_ContextItemNomotor_Hit\cm_m1_2x2\SPM.mat']};


matlabbatch{1}.spm.stats.results.spmmat = {'D:\1 [Context-Memory] fMRI Data\2 Second level results s4WithDerivAnts\BasicAdjust m_ci3_ContextItemNomotor_Hit\cm_m1_2x2\SPM.mat'}

matlabbatch{1}.spm.stats.results.conspec.titlestr = 'All';
matlabbatch{1}.spm.stats.results.conspec.contrasts = 5;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 1;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = true;

    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];
