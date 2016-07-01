function [ws, r] = i_context_memPcat(ws, r, ccol, memtype)
% [ws, r] = i_context_memPcat(ws, r, ccol, memtype)
% Create boxcar regressors for Context. 
% Contexts fall into Sim x Val cells (as conditions), with pmods describing
% whether that context trial is a good-memory (>=2 hits, label 2) or bad
% (label 1)
%
% Duration=Entire context presentation
% ------------------------------------------------------------------------

% Execute to debug: ccol=col.c; memscore4context=design.Memscore4context; memtype=memscore4context;

c_design={'SimR' 1 1; 'SimN' 1 0; 'DisR' 2 1; 'DisN' 2 0}; % Col 1: Name, Col 2= Sim status, Col 3=Reward status

if strcmp(memtype{1}, 'Hit')==0;  error('Onsets function & data have not yet been set up for requested memtype!'); end

%%

for d=1:size(c_design,1)
                
        % Sample
        wd=ws.c(ws.c(:,ccol.SimDis)==c_design{d,2} & ws.c(:,ccol.RewNeu)==c_design{d,3},:);
        
        % Main regressors
        ws.v.names{r}=c_design{d,1};
        ws.v.onsets{r}=wd(:,ccol.Context_Onset);
        ws.v.durations{r}=wd(:,ccol.Context_Duration);

        % Parametric modulators
        for p=1:length(memtype)
            ws.v.pmod(r).name{p}=['CMem_' memtype{p}];
            eval(['ws.v.pmod(r).param{p}=wd(:, ccol.Memall_' memtype{p} 'cat);']);
            ws.v.pmod(r).poly{p}=1;
        end

        r=r+1;
end


end

