function [ws, r ] =i_nointerest(ws, r, col, nointerestregressors )
% [ws, r ] =i_nointerest(ws, r, col, nointerestregressors )
% Create conditionsregressors of no interest
% 
%  Conditions are specified in 'nointerestregressors' (cell array)
%  Options: Outcome, Motor, Error    e.g. {'Outcome'; 'Motor'; 'Error'};
% --------------------------------------------------------------------------------------------

% Format data
ws.n.Outcome_Onset=ws.c(ws.c(:,col.c.Outcome_presented)==1,col.c.Outcome_Onset);
ws.n.Outcome_Duration=ws.c(ws.c(:,col.c.Outcome_presented)==1,col.c.Outcome_Duration);
ws.n.Outcome_Shown=ws.c(ws.c(:,col.c.Outcome_presented)==1,col.c.Outcome_amount);

ws.OutcomeYes = ws.c(ws.c(:,col.c.Outcome_presented)==1,:);
ws.OutcomeNo = ws.c(ws.c(:,col.c.Outcome_presented)==0,:);
ws.n.OutcomeYes_Onset=ws.OutcomeYes(:, col.c.Outcome_Onset);
ws.n.OutcomeYes_Duration=ws.OutcomeYes(:, col.c.Outcome_Duration);
ws.n.OutcomeNo_Onset=ws.OutcomeNo(:, col.c.Outcome_Onset);
ws.n.OutcomeNo_Duration=ws.OutcomeNo(:, col.c.Outcome_Duration);

ws.n.Motor_Onset=ws.i(isnan(ws.i(:,col.i.Item_MotorOnset))==0, col.i.Item_MotorOnset);
ws.n.Motor_Duration=zeros(size(ws.n.Motor_Onset));
ws.n.Error_Onset=ws.i(ws.i(:, col.i.Item_Error)==1, col.i.Item_Onset);
ws.n.Error_Duration=ws.i(ws.i(:, col.i.Item_Error)==1, col.i.Item_Duration);
ws.n.Outcome_Duration=zeros(size(ws.n.Outcome_Duration));  % Reset All Durations =0
ws.n.Error_Duration=zeros(size(ws.n.Error_Duration));

% Create regressors
for d=1:length(nointerestregressors)
    if strcmp(nointerestregressors{d}, 'Error')==1 && isempty(ws.n.Error_Onset)==1
        % If subject made no errors, exclude error regressor
    else
        ws.v.names{r}=nointerestregressors{d};
        eval(['ws.v.onsets{r}=ws.n.' nointerestregressors{d} '_Onset;'])
        eval(['ws.v.durations{r}=ws.n.' nointerestregressors{d} '_Duration;'])
        
        % Parametric modulators
        if strcmp(nointerestregressors{d}, 'Outcome')==1
            ws.v.pmod(r).name={'Outcome_shown'};
            ws.v.pmod(r).param{1}=ws.n.Outcome_Shown;
            ws.v.pmod(r).poly{1}=1;
        end
        r=r+1;
    end
end

end

