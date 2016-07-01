function [ ws, r] = m_i3_ItemContextpresonly( ws, r,  col, memscore4context, memscore4item)
%  [ ws, r] = m4_Item( ws, r,  col, memscore4context, memscore4item)
%  Model: Construct regressors for Items only
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

% (1) Construct Context regressors (presentation only)
[ws, r] = i_contextpresonly(ws, r, col.c);

% (2) Construct Item regressors
[ ws,r ] = i_item_memtype(ws, r, col.i, memscore4item);

% (3) Construct regressors of No Interest
[ws, r ] =i_nointerest(ws, r, col, {'Outcome'; 'Error'});
 
end

