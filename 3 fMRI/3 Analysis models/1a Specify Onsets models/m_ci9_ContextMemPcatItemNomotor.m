function [ ws, r] = m_ci9_ContextMemPcatItemNomotor( ws, r,  col, memscore4context, memscore4item)
% [ ws, r] = m_ci9_ContextMemPcatItemNomotor( ws, r,  col, memscore4context, memscore4item)
%  Model: Construct regressors for both Context (boxcars) and items
%             Context duration: Entire presentation time (Onset to Offset)
%             Context is split into Similarity x Valence, with Memory as 
%                 pmods designating memory category (Good memory: >=2 hits)
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

% Debug: memscore4context=design.Memscore4context; memscore4item=design.Memscore4item;

% (1) Construct Context regressors
[ws, r] = i_context_memPcat( ws,r , col.c, memscore4context);

% (2) Construct Item regressors
[ ws,r ] = i_item_memtype(ws, r, col.i, memscore4item);

% (3) Construct regressors of No Interest
[ws, r ] =i_nointerest(ws, r, col, {'Outcome'; 'Error'});

end

