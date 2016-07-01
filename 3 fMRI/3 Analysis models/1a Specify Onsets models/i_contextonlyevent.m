function [ws, r] = i_contextonly(ws, r, ccol, pmods)
% Create boxcar regressors for Context
% [ws, r] = i_context(ws, r, ccol, pmods)
%
% Duration=Entire context presentation
% ------------------------------------------------------------------------

c_design={'SimR' 1 1; 'SimN' 1 0; 'DisR' 2 1; 'DisN' 2 0}; % Col 1: Name, Col 2= Sim status, Col 3=Reward status

%%

for d=1:size(c_design,1)
    
    % Sample
    wd=ws.c(ws.c(:,ccol.SimDis)==c_design{d,2} & ws.c(:,ccol.RewNeu)==c_design{d,3},:);
    
    % Main regressors
    ws.v.names{r}=c_design{d,1};
    ws.v.onsets{r}=wd(:,ccol.Context_Onset);
    ws.v.durations{r}=wd(:,ccol.ContextOnly_Duration);
        
    % NO memory parametric modulators
    
    
    r=r+1;
end

end

