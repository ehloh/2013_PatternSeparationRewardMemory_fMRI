function [ ws,r ] = i_item_nostatus(ws, r, icol, memtype)
% Create item regressors, without status (i.e. item presentations only)
% [ ws,r ] = i_item_nostatus(ws, r, icol, memtype)
%
% Duration=0
%-------------------------------------------------------

%%

ws.v.names{r}='Items';
ws.v.onsets{r}=ws.i(:,icol.Item_Onset);
ws.v.durations{r}=ws.i(:, icol.Item_Duration);
ws.v.durations{r}=zeros(size(ws.v.onsets{r})); % Reset duration to 0

r=r+1;

end

