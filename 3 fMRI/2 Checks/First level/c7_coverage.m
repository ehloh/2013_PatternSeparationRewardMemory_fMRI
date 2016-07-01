function [ details] = c7_coverage(spm, ContextOrItem, memtype)
% [ details] = c7_checkfactorial_factorresults(spm, ContextOrItem, memtype)
% Find contrast weights corresponding to factor contrasts
%
%       spm:                    SPM.mat variable for model (after 1st level contrasts)
%       ContextorItem:     Are we assessing Context or Item regressors?
%       memtype:            Hit/Surehit/Rem/Surerem  - Roc not written yet
%
%       details:                 Cell 
%                                       Col 1: Name of Factor/comparison
%                                       Col 2: Regressor weights corresponding to this comparison
% ----------------------------------------------------------------------------------

% Execute: spm=ws.SPM; ContextOrItem=2; memtype=log.memtype;

details{1,1}='Coverage';
details{1,2}=zeros(1,size(spm.xCon(1).c,1));

end

