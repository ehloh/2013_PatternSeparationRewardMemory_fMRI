% cd('C:\Users\eloh\Desktop\2 [Explore]\2 Second level results s4Ants\m_v46g_ChoiceXvMargChoDiff_bpji08bpji11_Basic\TaskXChoicemargvcho_2x3\ROI\v46g')
cd('D:\1 [Context-Memory] fMRI Data\2 Second level results s4FullCardiacWithDerivAnts\m_i6_ItemMempmodContextpresonly_Hit_Landmarks5\im_m1_2x2\ROI')


f=spm_select('List', pwd, '.*.img$');
for i=1:size(f, 1)
    name=f(i, 1:strfind(f(i,:), '.img')-1);
    
    % Convert 
    v=spm_vol([name '.img']);
    ima=spm_read_vols(v);
    v.fname=[name '.nii'];
    spm_write_vol(v,ima);
    
    
end

ls




