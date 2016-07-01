% Re-slice anatomical ROIs to the fxnal data spex
clear all
 

where='F:\1 Context Mem study\1 All data'; cd(where)
% , subs=cellstr(ls('p*')); 
subs_manuNLM= {'p01_CW'; 'p03_EA';'p04_JL'; 'p07_LH';'p08_AM';'p09_CN';'p10_AB';'p11_SS';'p12_IL'; 'p14_SJ';'p16_TH';'p17_RA';'p18_KB';'p19_CN';'p20_JB';'p21_SH';'p22_EK';'p23_IS';'p24_LL';'p25_BS';'p26_MC';'p27_EW';'p28_CN';'p29_ET'}; 
subs=subs_manuNLM; 
whereto='D:\1 [Context-Memory] fMRI Data\1 MRI data';



do_makesubstrux=0;
if do_makesubstrux
    for s=1:length(subs)
        
        % Create average HPC scan
        f=cellstr(spm_select('List', [where filesep subs{s} '\1 Preprocessed\HPC'], 'sMQ.*.img'));
        matlabbatch{1}.spm.util.imcalc.input =cellfun(@(x)[where filesep subs{s} '\1 Preprocessed\HPC\' x ',1'], f, 'UniformOutput',0);
        matlabbatch{1}.spm.util.imcalc.output = [subs{s} '_MTL_coreg.img'];
        matlabbatch{1}.spm.util.imcalc.outdir = {[whereto filesep subs{s} filesep '1 Preprocessed']};
        matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3+i4)/4';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch); matlabbatch=[];
        
        % Coregister to T1w scan
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[whereto filesep subs{s} '\1 Preprocessed\' subs{s} '_T1w_coreg.nii,1']};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {[whereto filesep subs{s} '\1 Preprocessed\'  subs{s} '_MTL_coreg.img,1']};
        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch); matlabbatch=[];
        
        % Convert to nii
        cd([whereto filesep subs{s} '\1 Preprocessed' ])
        v=spm_vol([subs{s} '_MTL_coreg.img,1']);
        ima=spm_read_vols(v);
        v.fname=[subs{s} '_MTL_coreg.nii'];
        spm_write_vol(v,ima);
        delete([subs{s} '_MTL_coreg.img'])
        delete([subs{s} '_MTL_coreg.hdr'])
        
        %     spm_jobman('initcfg'); spm_jobman('run' , matlabbatch); matlabbatch=[];
    end
    error('done!')    
end






%% CHECK REG?
do_checkreg=0;
if do_checkreg
    subok=[subs num2cell(nan(length(subs),1))];
    for s=1:length(subs)
        matlabbatch{1}.spm.util.checkreg.data = {
            [whereto filesep subs{s} '\1 Preprocessed\' subs{s} '_MTL_coreg.nii,1']
            [whereto filesep subs{s} '\1 Preprocessed\' subs{s} '_T1w_coreg.nii,1']
            };
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch); matlabbatch=[];
        
        disp('-----------------------------------------------------')
        disp(subs(s))
        subok{s,2}=input('Input (1=OK, 0=Not ok) :   ');
    end
end

                                     
%% Move 
do_prepants=1; 
whereants='D:\host\Context_study\1d_Data_MTLpartial'; 
if do_prepants
    for s=1:length(subs)
        copyfile([whereto filesep subs{s} '\1 Preprocessed\' subs{s} '_MTL_coreg.nii'], ['D:\host\Context_study\1d_Data_MTLpartial\' subs{s} '_MTL_coreg.nii'] )
    end
end
