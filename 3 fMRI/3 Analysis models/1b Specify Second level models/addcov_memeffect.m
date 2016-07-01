function [cov] = addcov_memeffect(subjectlist, ncells, memdata)
% [cov] = addcov_memeffect(subjectlist, ncells, memdata)
% Include each subject's memory effect as covariates
%
% ------------------------------------------------------------------------------

%% Which memory types? 

memorytype={'dprime'; 'hitrate'; 'surehitrate'; 'rem'; 'know'; 'surerem'; 'sureknow';
                        'p_hitrem'; 'p_hitsure'; 'p_hitsurerem'; 'p_knowsure';'p_remsure'};
memeffect={'Overall'; 'Overall_Simfx'; 'Overall_Valfx'; 
                    'Sim_valfx'; 'Dis_valfx'; 'R_simfx'; 'N_simfx'}; 

% Input covariate types
disp('Which type of memory to use?'); disp([num2cell(1:length(memorytype))' memorytype]); whichmemtype=[]; while isempty(whichmemtype); whichmemtype=input('Input no:    '); disp(' '); end
disp('Which memory effect to use as the covariate?'); disp([num2cell(1:length(memeffect))' memeffect]); whichmemeffectype=[]; while isempty(whichmemeffectype); whichmemeffectype=input('Input no:    '); disp(' '); end
disp(['COVARIATE:     Memory type = ' char(memorytype(whichmemtype)) '      Memory effect type: ' char(memeffect(whichmemeffectype))])

% Name for this model
fxname=char(memeffect(whichmemeffectype));
cov.name=['memfx_' char(memorytype(whichmemtype)) '-' fxname];

%% Fetch requested covariate parameters

eval(['d=memdata.d_' char(memorytype(whichmemtype)) ';']);

% Alter name (to match behavioural data format)
if strcmp(fxname, 'Overall')==1;  varname='Overall.Overall';
elseif strcmp(fxname(1:7), 'Overall')==1; varname=['summ.' fxname];
else  varname=['cell.' fxname];
end
eval(['covars=d.' varname ';'])

% Apply subject selection
[subjectlist n_subjs CovarTable] = f_selectsubjects(vertcat({'Subjects' 'Cov Value'},[memdata.sub_info.Subject num2cell(covars)]),subjectlist,[vertcat({'Subjects' 'Cov Value'},[memdata.sub_info.Subject num2cell(covars)]) vertcat({'Include_OK'}, num2cell(ones(size(memdata.sub_info.Subject))))],'Include_OK');

%% Construct subject covariates (each subject's score in its own covariate)

for s=1:n_subjs
    c=zeros(n_subjs*ncells,1); 
    for i=1:ncells
        c((i-1)*n_subjs+s)=CovarTable{s+1,2};
    end
    cov.covbatch.cov(s).c = c; 
    cov.covbatch.cov(s).cname = [subjectlist{s} '_Memeffect_' char(memorytype(whichmemtype)) '-' fxname];
    cov.covbatch.cov(s).iCFI = 1;
    cov.covbatch.cov(s).iCC = 1;
end

%% Put covariate into the batch (alterate way of doing it, probably wrong - makes the main cell conditions un-estimable!)
% 
% Covars=(cell2mat(CovarTable(2:end, 2))*ones(1,ncells))';
% cov.covbatch.cov(1).c = vertcat(Covars(:));
% cov.covbatch.cov(1).cname = ['Memeffect_' char(memorytype(whichmemtype)) '-' fxname];
% cov.covbatch.cov(1).iCFI = 1;
% cov.covbatch.cov(1).iCC = 1;

end

