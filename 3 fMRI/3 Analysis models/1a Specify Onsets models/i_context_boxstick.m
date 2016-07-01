function [ws, r] = i_context_boxstick(ws, r, ccol, pmods, boxcarconds, stickconds)
% Create boxcar regressors for Context - specifying whether each context is
% to be modelled with a boxcar or a stick function
% [ws, r] = i_context_boxstick(ws, r, ccol, pmods, boxcarconds, stickconds)
%
% boxcarconds/stickconds: cell containing names of conditions to be
%                                       modelled as box cars or stick functions
%
% Context conditions: SimR, SimN, DisR,DisN
% Boxcar duration=Entire context presentation
% ------------------------------------------------------------------------

% Col 1: Name, Col 2= Sim & Reward status, Col 3=Boxcar?
c_design=vertcat([boxcarconds cell(length(boxcarconds),1) num2cell(ones(length(boxcarconds),1))], [stickconds cell(length(stickconds),1) num2cell(zeros(length(stickconds),1))]);
c_design{strcmp(c_design(:,1),'SimR'), 2}=[1 1];
c_design{strcmp(c_design(:,1),'SimN'), 2}=[1 0];
c_design{strcmp(c_design(:,1),'DisR'), 2}=[2 1];
c_design{strcmp(c_design(:,1),'DisN'), 2}=[2 0];
c_design=sortrows(c_design,-1);

%%

for d=1:size(c_design,1)
    
    % Sample
    wd=ws.c(ws.c(:,ccol.SimDis)==c_design{d,2}(1) & ws.c(:,ccol.RewNeu)==c_design{d,2}(2),:);
    
    % Main regressors
    ws.v.names{r}=c_design{d,1};
    ws.v.onsets{r}=wd(:,ccol.Context_Onset);
    
    % Model context with what function?
    if c_design{d,3}==1 % Boxcar 
        ws.v.durations{r}=wd(:,ccol.Context_Duration);
    elseif c_design{d,3}==0 % Stick function
        ws.v.durations{r}=zeros(size(wd,1),1);
    end
    
    % Parametric modulators
    for p=1:length(pmods)
        ws.v.pmod(r).name{p}=['CMem_' pmods{p}];
        eval(['ws.v.pmod(r).param{p}=wd(:, ccol.Memall_' pmods{p} 'score);']);
        ws.v.pmod(r).poly{p}=1;
    end
    
    r=r+1;
end


end

