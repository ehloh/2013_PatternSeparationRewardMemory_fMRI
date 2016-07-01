function [res] = f_calcmem(whichmeasure, responses, farate)
% [res] = f_calcmem(whichmeasure, responses, farate)
% Function takes an array of responses (from an OLD cell), and calculates the corresponding
% memory measure 
% 
% INPUTS
%       - Which measure: 
%               1=     dprime
%               2=     remember
%               3=     know
%               4=     sure remember
%               5=     sure know
%               6=     hit rate
%               7=     sure hit rate
%   Percentage scores (no FA rate needed)
%               8=     % hit sure (% of hits that are sure)
%               9=     % hit remember
%               10=    % hit sure remember
%               11=    % remember sure (% of remembers that are sure)
%               12=    % know sure (% of knows that are sure)
%       - Array of responses for the target cell (labelled 1-6; ROC points)
%             Response types
%                 1=Sure New
%                 2= Unsure New
%                 3= U K 
%                 4= S K
%                 5= U R
%                 6= S R
%       - FA rate 
%       - NOTE: numbers entered for dprime are corrected in this script for
%                  +/- infinite values (with a .001 correction)
%
% ---------------------------------------------------------------------------

% Valid measure?
if whichmeasure<1 || whichmeasure>12 || round(whichmeasure)~=whichmeasure
    input('ERROR: Invalid measure chosen. See function details for options');
end

% Valid sample?
ntrials=size(responses,1); res.n_trials=ntrials; res.rate=999; % Default=error
if ntrials~=0
    
%%    

if whichmeasure== 1
    
    % DPRIME  --------------------------------------------------------------
    res.n_hits=sum((responses(:)>2.99));
    res.hitrate=res.n_hits/ntrials;
    switch res.hitrate % Correct for infinite values
        case 0; res.hitrate=res.hitrate+0.001;
        case 1;res.hitrate=res.hitrate-0.001;
    end
    switch farate
        case 0; farate=farate+0.001;
        case 1;farate=farate-0.001;
    end
    %
    res.rate=norminv(res.hitrate)-norminv(farate);
    
elseif whichmeasure>1 && whichmeasure<3.001
    
    % REMEMBER/KNOW -------------------------------------------------------
    switch whichmeasure
        case 2 % remember;
            resp1=6; resp2=5;
        case 3 % know
            resp1=4; resp2=3;
    end
    res.n_hits=sum(responses(:)==resp1) + sum(responses(:)==resp2);
    res.rate_uncorrected=res.n_hits/ntrials;
    res.rate=res.rate_uncorrected-farate;
    
elseif whichmeasure>3 && whichmeasure<5.001
    
    % SURE REMEMBER/KNOW --------------------------------------------------
    switch whichmeasure; case 4 % sure remember
        resp=6; case 5 % sure know
        resp=4;
    end
    res.n_hits=sum(responses(:)==resp);
    res.rate_uncorrected=res.n_hits/ntrials;
    res.rate=res.rate_uncorrected-farate;
    
elseif whichmeasure==6
    
    % HIT RATE (CORRECTED) ---------------------------------------------------
    res.n_hits=sum(responses(:)>2);
    res.rate_uncorrected=res.n_hits/ntrials;
    res.rate=res.rate_uncorrected-farate;
    
elseif whichmeasure==7
    
    % SURE HIT RATE (CORRECTED) ---------------------------------------------------
    res.n_hits=sum(responses(:)==4) + sum(responses(:)==6);
    res.rate_uncorrected=res.n_hits/ntrials;
    res.rate=res.rate_uncorrected-farate;
    
%% Percentages ###################################

elseif whichmeasure==8 || whichmeasure==9 || whichmeasure==10
    
    %  P(Sure/Rem/SureRem | Hit)  ---------------------------------------------------
    res.n_base=sum(responses(:)>2);
    switch whichmeasure; case 8 % proportion of hits that are sure
        res.n_subset=sum(responses(:)==4)+sum(responses(:)==6);
        case 9 % proportion of hits that are remembered
            res.n_subset=sum(responses(:)==5)+sum(responses(:)==6);
        case 10 % proportion of hits that are sure-remembered
            res.n_subset=sum(responses(:)==6);
    end
    res.rate=res.n_subset/res.n_base;
    
elseif whichmeasure==11 || whichmeasure==12
    
    %  P(Sure | Rem/Know)  ---------------------------------------------------
    switch whichmeasure
        case 11  % proportion of remembers that are sure
            resp=5; respsure=6;
        case 12 % proportion of knows that are sure
            resp=3; respsure=4;
    end
    res.n_base=sum(responses(:)==resp)+sum(responses(:)==respsure);
    res.n_subset=sum(responses(:)==respsure);
    res.rate=res.n_subset/res.n_base;
    
end
end

end

