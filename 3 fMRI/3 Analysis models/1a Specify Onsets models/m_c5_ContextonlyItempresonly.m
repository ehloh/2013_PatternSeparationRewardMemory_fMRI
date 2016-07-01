function [ ws, r] = m_c5_ContextonlyItempresonly( ws, r,  col, memscore4context, memscore4item)
%  [ ws, r] = m_c5_ContextonlyItempresonly( ws, r,  col, memscore4context, memscore4item)
%  Model: Construct regressors for  Context (boxcars) 
%             Context duration: Context only (4s)
%
%  memscore4context: Cell array, with (ordered) memory types to apply to
%                                   Context regressors as Parametric modulators (e.g. 'Hit')
%                                   Options: see   col.c.Memall_ **
%                                               (Hit, Surehit, Rem, Surerem)
%                                   
%  memscore4item:      'Hit', 'Surehit', 'Rem', 'Surerem'
%
% To specify other memory criteria for item/context regressors, first
% create column with this memory score (at the data formatting stage)
% ------------------------------------------------------------------------------------------

% (1) Construct Context regressors
[ws, r] = i_contextonly( ws,r , col.c, memscore4context);

% (2) Item presentation only (No item status)
[ ws,r ] = i_item_nostatus(ws, r, col.i, memscore4item);

% (3) Construct regressors of No Interest
[ws, r ] =i_nointerest(ws, r, col, {'Outcome'; 'Motor'; 'Error'});
 
end

