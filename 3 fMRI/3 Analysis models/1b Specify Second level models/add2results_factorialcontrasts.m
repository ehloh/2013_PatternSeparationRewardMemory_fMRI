function [ cons] = add2results_factorialcontrasts( secondlevelfolder )
% [ cons] = add2results_factorialcontrasts( secondlevelfolder )
% Add other sensible contrasts, if it is a 2x2 (Sim x Val) Factorial design
%           (2x2x2 not allowed)
%
% ---------------------------------------------------------------------------------------


% Compile instructions
cons={'Negative effect of Similarity' [-1 -1 1 1];
            'Negative effect of Valence' [-1 1 -1 1];
            'Negative Interaction: Similarity x Valence' [-1 1 1 -1];
            'ME Sim_effect'          [1 1 0 0]; % Main effects (against baseline)
            'ME Dis_effect'          [0 0 1 1];
            'ME Val_effect'          [1 0 1 0];
            'ME Neu_effect'         [0 1 0 1];
            'SR_effect'                 [1 0 0 0];  % SR
            'SR_vs_others'          [3 -1 -1 -1];
            'SR_vs_SN'               [1 -1 0 0];
            'SR_vs_DR'               [1 0 -1 0];
            'SN_effect'                 [0 1 0 0];  % SN
            'SN_vs_others'          [-1 3 -1 -1];
            'SN_vs_SR'               [-1 1 0 0];
            'SN_vs_DN'               [0 1 0 -1];
            'DR_effect'                 [0 0 1 0]; % DR
            'DR_vs_others'          [-1 -1 3 -1];
            'DR_vs_DN'               [0 0 1 -1];
            'DN_effect'                 [0 0 0 1]; % DN
            'DN_vs_others'          [-1 -1 -1 3];
            'DN_vs_DR'               [0 0 -1 1];
            };

% Compile batch
matlabbatch{1}.spm.stats.con.spmmat =  {[secondlevelfolder filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.con.delete = 1;
for c=1:size(cons,1)
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = cons{c,1};
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = cons{c,2};
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
end

% Identity contrast
matlabbatch{1}.spm.stats.con.consess{c}.fcon.name = 'identity';
matlabbatch{1}.spm.stats.con.consess{c}.fcon.convec = eye(4);
matlabbatch{1}.spm.stats.con.consess{c}.fcon.sessrep = 'none';

% Run SPM
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];

end

