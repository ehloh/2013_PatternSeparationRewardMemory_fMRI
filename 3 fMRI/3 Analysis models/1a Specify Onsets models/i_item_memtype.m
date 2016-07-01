function [ ws,r ] = i_item_memtype(ws, r, icol, memtype)
% Create item regressors (SimDis + RewNeu+ MemOrNot)
% [ ws,r ] = i_item_memtype(ws, r, icol, memtype)
%
% memtype : 'Hit', 'Surehit', 'Rem', 'Surerem'
%
% Duration=0
%-------------------------------------------------------

c_design={'SR' 1 1; 'SN' 1 0; 'DR' 2 1; 'DN' 2 0}; % Col 1: Name, Col 2= Sim status, Col 3=Reward status

% Memory criteria
if strcmp(memtype,'Roc')==1;  error('Requested memtype==Roc. Roc memory cannot be used for 2x2x2 (Sim x Val x Mem) item regressors!'); end
eval(['col_memtype=icol.Mem_' memtype 'Or;']);
switch memtype
    case 'Hit'
        memyes='Hit';
        memno='Miss';
    case 'Surehit'
        memyes='Surehit';
        memno='NotSurehit';
    case 'Rem'
        memyes='Rem';
        memno='NotRem';
    case 'Surerem'
        memyes='Surerem';
        memno='NotSurerem';
    otherwise
        error(['Could not find requested memory type  -- ' memtype])
end


%%

for d=1:size(c_design,1)
    
    % Compile Context's Item samples
    wd=ws.i(ws.i(:,icol.SimDis)==c_design{d,2} & ws.i(:,icol.RewNeu)==c_design{d,3},:);
    
    for i=2:-1:1
        m=i-1;
        
        % Compile Items sample with correct Context & Memory status
        wm=wd(wd(:,col_memtype)==m,:);
        
        % Main regressors
        switch m
            case 0
                ws.v.names{r}=[c_design{d,1} '_' memno];
            case 1
                ws.v.names{r}=[c_design{d,1} '_' memyes];
        end
        ws.v.onsets{r}=wm(:,icol.Item_Onset);
        ws.v.durations{r}=wm(:,icol.Item_Duration);
        ws.v.durations{r}=zeros(size(ws.v.onsets{r})); % Reset duration to 0
        
        % Parametric modulators
        %         ws.v.pmod(r).name=['Mem_' pmods{p}];
        %         eval(['ws.v.pmod(r).param{1}=wd(:, dcol.Memall_' pmods{p} ');']);
        %         ws.v.pmod(r).poly{1}=1;
        
        %
        r=r+1;
        wm=[];
    end
    
    wd=[];
end
end

