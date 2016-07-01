% 

%% (1) Memory statistics: Run eval2_memory_statistics first

memtypes={'dprime' 'dpr';'hitrate' 'hr';'surehitrate' 'shr';'rem' 're';'surerem' 'sre'};

% Data variable
dd=cell(length(alldata.sub_info.Subject)+1,5);
dd(2:end,1)=alldata.sub_info.Subject(:); dd{1,1}='Subject'; % Subjects

% Titles
var={'Overall.Overall';
    'cell.R_simfx'; 'cell.N_simfx';'cell.Sim_valfx';'cell.Dis_valfx';
    'Cell.Sim_R';'Cell.Sim_N';'Cell.Dis_R';'Cell.Dis_N';
    };

k=2;

for m=1:size(memtypes,1)
    eval(['mem=d_' memtypes{m,1} ';'])
    
    for i=1:length(var)
        dd{1,k}=[memtypes{m,2} '_' var{i}];
        %
        eval(['w=mem.' var{i} ';'])
        
        dd(2:end,k)=num2cell(w);
        
        
        
        %
        k=k+1;
    end
end

%%
