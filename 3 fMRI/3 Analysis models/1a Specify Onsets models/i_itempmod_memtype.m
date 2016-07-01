function [ ws,r ] = i_itempmod_memtype(ws, r, icol, memtype)
% Create item regressors (SimDis + RewNeu+ MemOrNot)
% [ ws,r ] = i_item_memtype(ws, r, icol, memtype)
%
% memtype : 'Hit', 'Surehit', 'Rem', 'Surerem'
%
% Duration=0
%-------------------------------------------------------

c_design={'SR' 1 1; 'SN' 1 0; 'DR' 2 1; 'DN' 2 0}; % Col 1: Name, Col 2= Sim status, Col 3=Reward status

% Memory criteria
if strcmp(memtype, 'Roc')==1
    eval(['col_memtype=icol.Mem_' memtype ';']);
elseif strcmp(memtype, 'Hit')==1 || strcmp(memtype, 'Surehit')==1 || strcmp(memtype, 'Rem')==1 || strcmp(memtype, 'Surerem')==1
    eval(['col_memtype=icol.Mem_' memtype 'Or;']);
else
    error(['Could not find requested memory type  -- ' memtype])
end

%%

for d=1:size(c_design,1)
    
    % Compile Context's Item samples
    wd=ws.i(ws.i(:,icol.SimDis)==c_design{d,2} & ws.i(:,icol.RewNeu)==c_design{d,3},:);
    
    % Compile Item onsets with correct Context & Valence status
    ws.v.names{r}=[c_design{d,1} '_Item'];
    ws.v.onsets{r}=wd(:,icol.Item_Onset);
    ws.v.durations{r}=wd(:,icol.Item_Duration);
    ws.v.durations{r}=zeros(size(ws.v.onsets{r})); % Reset duration to 0
    
    % Parametric modulators
    ws.v.pmod(r).name{1}=['Mem_' memtype];
    ws.v.pmod(r).param{1}=wd(:, col_memtype);
%     eval(['ws.v.pmod(r).param{1}=wd(:, icol.Mem_' memtype ');']);
    ws.v.pmod(r).poly{1}=1;
    
    %
    r=r+1;
    wd=[];
end

end

