function [ ws, r] = m_ci10_ContextEventItemNomotor( ws, r,  col, memscore4context, memscore4item)
%  [ ws, r] = m_ci10_ContextEventItemNomotor( ws, r,  col, memscore4context, memscore4item)
%  Model: Construct regressors for both Context (boxcars) and items
%             Context duration: Entire presentation time (Onset to Offset)
%
%  No context memory parametric modulators included
%                                   
%  memscore4item:      'Hit', 'Surehit', 'Rem', 'Surerem'
%
% To specify other memory criteria for item/context regressors, first
% create column with this memory score (at the data formatting stage)
% ------------------------------------------------------------------------------------------

% (1) Construct Context regressors
[ws, r] = i_contextevent( ws,r , col.c);

% (2) Construct Item regressors
[ ws,r ] = i_item_memtype(ws, r, col.i, memscore4item);

% (3) Construct regressors of No Interest
[ws, r ] =i_nointerest(ws, r, col, {'Outcome'; 'Error'});

end

