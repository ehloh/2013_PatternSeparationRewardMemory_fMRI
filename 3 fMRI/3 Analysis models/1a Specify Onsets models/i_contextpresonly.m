function [ws, r] = i_contextpresonly(ws, r, ccol)
% Create boxcar regressors for Context, presentation only (No Sim or Val
% status, no mem pmods). 
% [ws, r] = i_context(ws, r, ccol, pmods)
%
% Duration=Entire context presentation
% ------------------------------------------------------------------------

ws.v.names{r}='ContextPresonly';
ws.v.onsets{r}=ws.c(:,ccol.Context_Onset);
ws.v.durations{r}=ws.c(:,ccol.Context_Duration);

r=r+1;

end

