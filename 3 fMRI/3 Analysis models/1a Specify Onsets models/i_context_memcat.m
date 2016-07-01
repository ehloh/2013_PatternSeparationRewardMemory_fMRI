function [ws, r] = i_context_memcat(ws, r, ccol, memtype)
% [ws, r] = i_context_memcat(ws, r, ccol, memtype)
% Create boxcar regressors for Context. 
% Contexts fall into Sim x Val x MemCat cells (as conditions)
%
% Duration=Entire context presentation
% ------------------------------------------------------------------------

% Execute to debug: ccol=col.c; memscore4context=request.memtype; memtype=memscore4context;

c_design={'SimR' 1 1; 'SimN' 1 0; 'DisR' 2 1; 'DisN' 2 0}; % Col 1: Name, Col 2= Sim status, Col 3=Reward status

if strcmp(memtype{1}, 'Hit')==0;  error('Onsets function & data have not yet been set up for requested memtype!'); end
eval(['col4memscore= ccol.Memall_' memtype{1} 'cat;'])

%%

for d=1:size(c_design,1)
    for L=1:2
                
        % Sample
        wd=ws.c(ws.c(:,ccol.SimDis)==c_design{d,2} & ws.c(:,ccol.RewNeu)==c_design{d,3} & ws.c(:,col4memscore)==L,:);
        if isempty(wd)==1; disp(['Empty SimxValxMemcat cell:   ' c_design{d,1} '_' num2str(L)]); end
        disp([c_design{d,1} '_m' num2str(L) ':   ' num2str(size(wd,1)) ' events'])
            
        % Main regressors
        switch L
            case 1
                ws.v.names{r}=[c_design{d,1} '_L'];
            case 2
                ws.v.names{r}=[c_design{d,1} '_H'];
        end
        ws.v.onsets{r}=wd(:,ccol.Context_Onset);
        ws.v.durations{r}=wd(:,ccol.Context_Duration);

%         % Parametric modulators
%         for p=1:length(memtype)
%             ws.v.pmod(r).name{p}=['CMem_' memtype{p}];
%             eval(['ws.v.pmod(r).param{p}=wd(:, ccol.Memall_' memtype{p} 'score);']);
%             ws.v.pmod(r).poly{p}=1;
%         end

        r=r+1;
        
    end
end


end

