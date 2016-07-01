function [ ws,r ] = i_item_memtype(ws, r, icol, memtype)
% Create item regressors (SimDis + RewNeu+ MemOrNot)
% [ ws,r ] = i_item_memtype(ws, r, icol, memtype)
%
% memtype : 'Hit', 'Surehit', 'Rem', 'Surerem'
%
% Duration=0
%-------------------------------------------------------

c_design={'SR' 1 1; 'SN' 1 0; 'DR' 2 1; 'DN' 2 0}; % Col 1: Name, Col 2= Sim status, Col 3=Reward status

%% Items by Sim x Val, but ignoring memhit or not.



for d=1:size(c_design,1)
    
    % Compile Context's Item samples
    wd=ws.i(ws.i(:,icol.SimDis)==c_design{d,2} & ws.i(:,icol.RewNeu)==c_design{d,3},:);
    
    
    ws.v.names{r}=[c_design{d,1} '_Item'];
    ws.v.onsets{r}=wd(:,icol.Item_Onset);
    ws.v.durations{r}=wd(:,icol.Item_Duration);
    ws.v.durations{r}=zeros(size(ws.v.onsets{r})); % Reset duration to 0
    
    %
    r=r+1;
    
    wd=[];
end
end

