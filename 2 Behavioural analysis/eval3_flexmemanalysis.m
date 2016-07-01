% eval3_flexmemanalysis
clear all; close all hidden; clc

for o1=1:1 % Options & Setup for analysis
    
    where.where='D:\Dropbox\SCRIPPS\2a ContextMem behaviour'; 
    where.scripts=[where.where filesep '1 Analysis scripts'];
    where.data='D:\Dropbox\SCRIPPS\1 ContextMem fMRI\2 Behavioural data';
    %
    w.minimum_ntrials=1;
    w.extravar=14;
    w.print2text=1;
    
    % Compile settings
    cd(where.data); w.dir=dir('p*'); for i=1:size(w.dir,1); details.subjects{i,1}=w.dir(i).name; end
    details.n_subjs=size(w.dir,1);
    addpath(where.scripts)
end

%% Load all data + Mark ROC (Col 31)

for o1=1:1 % Column settings
    col.et.SimDis=3; % original trialstats
    col.et.RewNeu=4;
    col.et.Trialnum=5;
    col.et.EncSession=6;
    col.et.OutcomeMag=7;
    col.et.OutcomeItem1=8;
    col.et.OutcomeItem2=9;
    col.et.OutcomeItem3=10;
    col.et.StimItem=[11 12 13];
    col.et.OutcomePres=29;
    col.et.AccItem=[26 27 28];
    col.et.ItemMemHit=[30 31 32];
    col.et.ItemMemSurehit=[33 34 35];
    col.et.ContextMemScoreHit=36;
    col.et.ContextMemScoreSurehit=37;
    col.et.OldFoil=38; % this is meaningless
    
    % modified trialstats (1 item per row)
    col.e.ItemStim=1;
    col.e.EncItemSerialPos=2;
    col.e.SimDis=3;
    col.e.RewNeu=4;
    col.e.Trialnum=5;
    col.e.EncSession=6;
    col.e.EncOutcomeMag=7;
    col.e.EncOutcomePres=9;
    col.e.AccItem=10;
    
    % Memory
    col.m.ItemStim=2;
    col.m.OldFoil=3;
    col.m.SimDis=10;
    col.m.RewNeu=11;
    col.m.OldNew=12; % response
    col.m.SureGuess=16;
    col.m.Roc=31;
    col.m.EncItemSerialPos=32;
    col.m.EncOutcomePres=33;
    col.m.EncOutcomeMag=34;
    col.m.EncCorrect=35;
    col.m.EncTrial=36;
    col.m.EncSession=37;
    
end

% 'subjdata': Col 1=Mem and enc trialstats, Col 2=Trialstats, old + correctly judged
%                   during encoding, Col 3=Trialstats for foils, Col
%                   4=Raw cell-splits, Col 5=Reorganized cells (design)
disp('Loading data + marking extra details ###################################')
subjdata=cell(details.n_subjs,3);

for s=1:details.n_subjs
    ws.d=load([where.data filesep details.subjects{s} filesep details.subjects{s} '_file_3memorytest.mat']);
    ws.dat=ws.d.memtest.trialstats;
    ws.p=ws.d.memtest.TPQ;
    ws.e1=load([where.data filesep details.subjects{s} filesep details.subjects{s} '_file_2encodingFMRI_b1.mat']);
    ws.e2=load([where.data filesep details.subjects{s} filesep details.subjects{s} '_file_2encodingFMRI_b2.mat']);
    
    % Encoding trials
    ws.e1.encoding.trialstats(:,col.et.EncSession)=1; ws.e2.encoding.trialstats(:,col.et.EncSession)=2; % 
    ws.et=[ws.e1.encoding.trialstats; ws.e2.encoding.trialstats];
    ws.enc=[]; % 1 item per row, col.e
    for e=1:3
        we.t(1:size(ws.et,1), col.e.EncItemSerialPos)=e;
        we.t(:, col.e.ItemStim)=ws.et(:,col.et.StimItem(e));
        we.t(:, col.e.SimDis)=ws.et(:,col.et.SimDis);
        we.t(:, col.e.RewNeu)=ws.et(:,col.et.RewNeu);
        we.t(:, col.e.Trialnum)=ws.et(:,col.et.Trialnum);
        we.t(:, col.e.EncSession)=ws.et(:,col.et.EncSession);
        we.t(:, col.e.EncOutcomeMag)=ws.et(:,col.et.OutcomeMag);
        we.t(:, col.e.EncOutcomePres)=ws.et(:,col.et.OutcomePres);
        we.t(:,col.e.AccItem)=ws.et(:,col.et.AccItem(e));
        ws.enc=[ws.enc;we.t];
        we=[];
    end
    
    % Mark context mem scores (in encoding trialstats)
    for t=1:size(ws.et,1) % Fetch item mem scores
        for i=1:3
            if length(find(ws.et(t,col.et.StimItem(i))==ws.dat(:,col.m.ItemStim)))==1
                ws.et(t, col.et.ItemMemHit(i))=2-ws.dat(find(ws.et(t,col.et.StimItem(i))==ws.dat(:,col.m.ItemStim)),col.m.OldNew);
                ws.et(t, col.et.ItemMemSurehit(i))=ws.et(t, col.et.ItemMemHit(i))*(2-ws.dat(find(ws.et(t,col.et.StimItem(i))==ws.dat(:,col.m.ItemStim)),col.m.SureGuess));
                if ws.et(t, col.et.AccItem(i))==0; % Exclude if Enc incorrectly
                    ws.et(t, col.et.ItemMemHit(i))=nan;
                end
            else error('wtf, theres a bug')
            end
        end
    end
    ws.et(:,col.et.ContextMemScoreHit)=nanmean(ws.et(:, col.et.ItemMemHit(:)),2);
    ws.et(:,col.et.ContextMemScoreSurehit)=nanmean(ws.et(:, col.et.ItemMemSurehit(:)),2);
    ws.et(:,col.et.OldFoil)=1;
    
    % Mark ROC type (Col 31)
    %             Response types
    %                 1=Sure New
    %                 2= Unsure New
    %                 3= U K 
    %                 4= S K
    %                 5= U R
    %                 6= S R
    ws.dat(ws.dat(:,14)==1 & ws.dat(:,16)==1,col.m.Roc)=6; % Rem
    ws.dat(ws.dat(:,14)==1 & ws.dat(:,16)==2,col.m.Roc)=5;
    ws.dat(ws.dat(:,14)==2 & ws.dat(:,16)==1,col.m.Roc)=4; % Know
    ws.dat(ws.dat(:,14)==2 & ws.dat(:,16)==2,col.m.Roc)=3;
    ws.dat(ws.dat(:,12)==2 & ws.dat(:,16)==2,col.m.Roc)=2; % New
    ws.dat(ws.dat(:,12)==2 & ws.dat(:,16)==1,col.m.Roc)=1;
    
    % Mark other encoding details
    for t=1:size(ws.dat,1)
        if ws.dat(t,3)==1
            if length(find(ws.enc(:,col.e.ItemStim) == ws.dat(t,col.m.ItemStim)))==1
                ws.dat(t,col.m.EncItemSerialPos)=ws.enc(find(ws.enc(:,col.e.ItemStim) == ws.dat(t,col.m.ItemStim)),col.e.EncItemSerialPos);
                ws.dat(t,col.m.EncOutcomePres)=ws.enc(find(ws.enc(:,col.e.ItemStim) == ws.dat(t,col.m.ItemStim)),col.e.EncOutcomePres);
                ws.dat(t,col.m.EncOutcomeMag)=ws.enc(find(ws.enc(:,col.e.ItemStim) == ws.dat(t,col.m.ItemStim)),col.e.EncOutcomeMag);
                ws.dat(t,col.m.EncCorrect)=ws.enc(find(ws.enc(:,col.e.ItemStim) == ws.dat(t,col.m.ItemStim)),col.e.AccItem);
                ws.dat(t,col.m.EncTrial)=ws.enc(find(ws.enc(:,col.e.ItemStim) == ws.dat(t,col.m.ItemStim)),col.e.Trialnum);
                ws.dat(t,col.m.EncSession)=ws.enc(find(ws.enc(:,col.e.ItemStim) == ws.dat(t,col.m.ItemStim)),col.e.EncSession);
            else error('Wrong no. of items matched up (finding items in encoding trials)')
            end
        else
            ws.dat(t,col.m.EncItemSerialPos)=nan;
            ws.dat(t,col.m.EncOutcomePres)=nan;
            ws.dat(t,col.m.EncOutcomeMag)=nan;
            ws.dat(t,col.m.EncCorrect)=nan;
        end
    end
    
    %
    subjdata{s,1}.mem=ws.dat;
    subjdata{s,1}.enc=ws.et;
    subjdata{s,1}.p.NS=[mean(ws.p.NS) ws.p.NS']; % mean, 1-2-3-4
    subjdata{s,1}.p.RD=[mean(ws.p.RD) ws.p.RD'];
    subjdata{s,2}=ws.dat(ws.dat(:,3)==1 & ws.dat(:,29)==1, :);
    subjdata{s,3}=ws.dat(ws.dat(:,3)==2,:);
    %
    ws=[];
end

%% Settings + prep 

% Calculate FA rates
for s=1:details.n_subjs
    subjdata{s,4}.fa_dprime=(sum(subjdata{s,3}(:,12)==1)+0.5)/size(subjdata{s,3},1); % Dprime
    subjdata{s,4}.fa_hitrate=sum(subjdata{s,3}(:,12)==1)/size(subjdata{s,3},1); % Hit rate
    subjdata{s,4}.fa_surehitrate=sum(subjdata{s,3}(:,12)==1 & subjdata{s,3}(:,16)==1)/size(subjdata{s,3},1);
    subjdata{s,4}.fa_rem=sum(subjdata{s,3}(:,14)==1)/size(subjdata{s,3},1); % Rem
    subjdata{s,4}.fa_surerem=sum(subjdata{s,3}(:,14)==1 & subjdata{s,3}(:,16)==1)/size(subjdata{s,3},1);
    subjdata{s,4}.fa_know=sum(subjdata{s,3}(:,14)==2)/size(subjdata{s,3},1); % Know
    subjdata{s,4}.fa_sureknow=sum(subjdata{s,3}(:,14)==2 & subjdata{s,3}(:,16)==1)/size(subjdata{s,3},1);
end

% Instructions for splitting cells 
for o1=1:1
    
    % Item memory
    instruc.cells={
        
        % For overall dprime
        'OverallMem'           {'OldFoil' 1;}

        % Sim x Dis x Serial Position
        'SR_pos1'       {'SimDis' 1;'RewNeu' 1;'EncItemSerialPos' 1};
        'SR_pos2'       {'SimDis' 1;'RewNeu' 1;'EncItemSerialPos' 2};
        'SR_pos3'       {'SimDis' 1;'RewNeu' 1;'EncItemSerialPos' 3};
        'SN_pos1'       {'SimDis' 1;'RewNeu' 0;'EncItemSerialPos' 1};
        'SN_pos2'       {'SimDis' 1;'RewNeu' 0;'EncItemSerialPos' 2};
        'SN_pos3'       {'SimDis' 1;'RewNeu' 0;'EncItemSerialPos' 3};
        'DR_pos1'       {'SimDis' 2;'RewNeu' 1;'EncItemSerialPos' 1};
        'DR_pos2'       {'SimDis' 2;'RewNeu' 1;'EncItemSerialPos' 2};
        'DR_pos3'       {'SimDis' 2;'RewNeu' 1;'EncItemSerialPos' 3};
        'DN_pos1'       {'SimDis' 2;'RewNeu' 0;'EncItemSerialPos' 1};
        'DN_pos2'       {'SimDis' 2;'RewNeu' 0;'EncItemSerialPos' 2};
        'DN_pos3'       {'SimDis' 2;'RewNeu' 0;'EncItemSerialPos' 3};

        % Sim x Dis x Outcome Presented?
        'SR_oY'        {'SimDis' 1;'RewNeu' 1;'EncOutcomePres' 1};
        'SR_oN'       {'SimDis' 1;'RewNeu' 1;'EncOutcomePres' 0};
        'SN_oY'       {'SimDis' 1;'RewNeu' 0;'EncOutcomePres' 1};
        'SN_oN'       {'SimDis' 1;'RewNeu' 0;'EncOutcomePres' 0};
        'DR_oY'       {'SimDis' 2;'RewNeu' 1;'EncOutcomePres' 1};
        'DR_oN'       {'SimDis' 2;'RewNeu' 1;'EncOutcomePres' 0};
        'DN_oY'       {'SimDis' 2;'RewNeu' 0;'EncOutcomePres' 1};
        'DN_oN'       {'SimDis' 2;'RewNeu' 0;'EncOutcomePres' 0};

        % Sim x Dis x Outcome Presented x Serial Position
        'SR_oY_pos1'        {'SimDis' 1;'RewNeu' 1;'EncOutcomePres' 1; 'EncItemSerialPos' 1};
        'SR_oY_pos2'        {'SimDis' 1;'RewNeu' 1;'EncOutcomePres' 1; 'EncItemSerialPos' 2};
        'SR_oY_pos3'        {'SimDis' 1;'RewNeu' 1;'EncOutcomePres' 1; 'EncItemSerialPos' 3};
        'SR_oN_pos1'       {'SimDis' 1;'RewNeu' 1;'EncOutcomePres' 0; 'EncItemSerialPos' 1};
        'SR_oN_pos2'       {'SimDis' 1;'RewNeu' 1;'EncOutcomePres' 0; 'EncItemSerialPos' 2};
        'SR_oN_pos3'       {'SimDis' 1;'RewNeu' 1;'EncOutcomePres' 0; 'EncItemSerialPos' 3};
        'SN_oY_pos1'       {'SimDis' 1;'RewNeu' 0;'EncOutcomePres' 1; 'EncItemSerialPos' 1};
        'SN_oY_pos2'       {'SimDis' 1;'RewNeu' 0;'EncOutcomePres' 1; 'EncItemSerialPos' 2};
        'SN_oY_pos3'       {'SimDis' 1;'RewNeu' 0;'EncOutcomePres' 1; 'EncItemSerialPos' 3};
        'SN_oN_pos1'       {'SimDis' 1;'RewNeu' 0;'EncOutcomePres' 0; 'EncItemSerialPos' 1};
        'SN_oN_pos2'       {'SimDis' 1;'RewNeu' 0;'EncOutcomePres' 0; 'EncItemSerialPos' 2};
        'SN_oN_pos3'       {'SimDis' 1;'RewNeu' 0;'EncOutcomePres' 0; 'EncItemSerialPos' 3};
        'DR_oY_pos1'       {'SimDis' 2;'RewNeu' 1;'EncOutcomePres' 1; 'EncItemSerialPos' 1};
        'DR_oY_pos2'       {'SimDis' 2;'RewNeu' 1;'EncOutcomePres' 1; 'EncItemSerialPos' 2};
        'DR_oY_pos3'       {'SimDis' 2;'RewNeu' 1;'EncOutcomePres' 1; 'EncItemSerialPos' 3};
        'DR_oN_pos1'       {'SimDis' 2;'RewNeu' 1;'EncOutcomePres' 0; 'EncItemSerialPos' 1};
        'DR_oN_pos2'       {'SimDis' 2;'RewNeu' 1;'EncOutcomePres' 0; 'EncItemSerialPos' 2};
        'DR_oN_pos3'       {'SimDis' 2;'RewNeu' 1;'EncOutcomePres' 0; 'EncItemSerialPos' 3};
        'DN_oY_pos1'       {'SimDis' 2;'RewNeu' 0;'EncOutcomePres' 1; 'EncItemSerialPos' 1};
        'DN_oY_pos2'       {'SimDis' 2;'RewNeu' 0;'EncOutcomePres' 1; 'EncItemSerialPos' 2};
        'DN_oY_pos3'       {'SimDis' 2;'RewNeu' 0;'EncOutcomePres' 1; 'EncItemSerialPos' 3};
        'DN_oN_pos1'       {'SimDis' 2;'RewNeu' 0;'EncOutcomePres' 0; 'EncItemSerialPos' 1};
        'DN_oN_pos2'       {'SimDis' 2;'RewNeu' 0;'EncOutcomePres' 0; 'EncItemSerialPos' 2};
        'DN_oN_pos3'       {'SimDis' 2;'RewNeu' 0;'EncOutcomePres' 0; 'EncItemSerialPos' 3};

        % Sim x Dis x OutcomeMag
        'SR_om0'       {'SimDis' 1;'RewNeu' 1;'EncOutcomeMag' 0};
        'SR_om1'       {'SimDis' 1;'RewNeu' 1;'EncOutcomeMag' 1};
        'SR_om2'       {'SimDis' 1;'RewNeu' 1;'EncOutcomeMag' 2};
        'SR_om3'       {'SimDis' 1;'RewNeu' 1;'EncOutcomeMag' 3};
        'SN_om0'       {'SimDis' 1;'RewNeu' 0;'EncOutcomeMag' 0};
        % 'SN_om1'       {'SimDis' 1;'RewNeu' 0;'EncOutcomeMag' 1};
        % 'SN_om2'       {'SimDis' 1;'RewNeu' 0;'EncOutcomeMag' 2};
        % 'SN_om3'       {'SimDis' 1;'RewNeu' 0;'EncOutcomeMag' 3};
        'DR_om0'       {'SimDis' 2;'RewNeu' 1;'EncOutcomeMag' 0};
        'DR_om1'       {'SimDis' 2;'RewNeu' 1;'EncOutcomeMag' 1};
        'DR_om2'       {'SimDis' 2;'RewNeu' 1;'EncOutcomeMag' 2};
        'DR_om3'       {'SimDis' 2;'RewNeu' 1;'EncOutcomeMag' 3};
        'DN_om0'       {'SimDis' 2;'RewNeu' 0;'EncOutcomeMag' 0};
        % 'DN_om1'       {'SimDis' 2;'RewNeu' 0;'EncOutcomeMag' 1};
        % 'DN_om2'       {'SimDis' 2;'RewNeu' 0;'EncOutcomeMag' 2};
        % 'DN_om3'       {'SimDis' 2;'RewNeu' 0;'EncOutcomeMag' 3};

        % Sim x Dis x Shown OutcomeMag
        'SR_som0'       {'SimDis' 1;'RewNeu' 1;'EncOutcomePres' 1; 'EncOutcomeMag' 0};  %(compare with SR/SN/DR/DN '_oN')
        'SR_som1'       {'SimDis' 1;'RewNeu' 1;'EncOutcomePres' 1; 'EncOutcomeMag' 1};
        'SR_som2'       {'SimDis' 1;'RewNeu' 1;'EncOutcomePres' 1; 'EncOutcomeMag' 2};
        'SR_som3'       {'SimDis' 1;'RewNeu' 1;'EncOutcomePres' 1; 'EncOutcomeMag' 3};
        'SN_som0'       {'SimDis' 1;'RewNeu' 0;'EncOutcomePres' 1; 'EncOutcomeMag' 0};
        % 'SN_som1'       {'SimDis' 1;'RewNeu' 0;'EncOutcomePres' 1; 'EncOutcomeMag' 1};
        % 'SN_som2'       {'SimDis' 1;'RewNeu' 0;'EncOutcomePres' 1; 'EncOutcomeMag' 2};
        % 'SN_som3'       {'SimDis' 1;'RewNeu' 0;'EncOutcomePres' 1; 'EncOutcomeMag' 3};
        'DR_som0'       {'SimDis' 2;'RewNeu' 1;'EncOutcomePres' 1; 'EncOutcomeMag' 0};
        'DR_som1'       {'SimDis' 2;'RewNeu' 1;'EncOutcomePres' 1; 'EncOutcomeMag' 1};
        'DR_som2'       {'SimDis' 2;'RewNeu' 1;'EncOutcomePres' 1; 'EncOutcomeMag' 2};
        'DR_som3'       {'SimDis' 2;'RewNeu' 1;'EncOutcomePres' 1; 'EncOutcomeMag' 3};
        'DN_som0'       {'SimDis' 2;'RewNeu' 0;'EncOutcomePres' 1; 'EncOutcomeMag' 0};
        % 'DN_som1'       {'SimDis' 2;'RewNeu' 0;'EncOutcomePres' 1; 'EncOutcomeMag' 1};
        % 'DN_som2'       {'SimDis' 2;'RewNeu' 0;'EncOutcomePres' 1; 'EncOutcomeMag' 2};
        % 'DN_som3'       {'SimDis' 2;'RewNeu' 0;'EncOutcomePres' 1; 'EncOutcomeMag' 3};

    };
    instruc.mem={'dprime';'rem';'know';'surerem';'sureknow';'hitrate';'surehitrate';'phitsure';'phitrem';'phitsurerem';'premsure';'pknowsure'}; % Order must match up with f_calcmem

% Instructions for context memory
 instruc.contextcells={        
        % For overall dprime
        'Overall'           {'OldFoil' 1;}

        % Sim x Dis x Serial Position
        'SimR'       {'SimDis' 1;'RewNeu' 1};
        'SimN'       {'SimDis' 1;'RewNeu' 0};
        'DisR'       {'SimDis' 2;'RewNeu' 1};
        'DisN'       {'SimDis' 2;'RewNeu' 0};
        };
end
    
%% Memory scores calculated from mem trialstats (i.e. each item by itself)

% Prep output data tables
w.table=cell(details.n_subjs+1, size(instruc.cells,1)+1);w.table{1,1}='Subject'; w.table(1,2:end)=instruc.cells(:,1);w.table(2:end,1)=details.subjects;
d_personality=[{'NSmean';'NS1';'NS2';'NS3';'NS4';'RDmean';'RD1';'RD2';'RD3';'RD4'}'; cell(details.n_subjs,10)];  d_overalldprime=['overalldprime'; cell(details.n_subjs,1)];
d_dprime=w.table; d_rem=w.table; d_surerem=w.table; d_know=w.table; d_sureknow=w.table; d_hit=w.table; d_surehit=w.table; d_phitsure=w.table; d_phitrem=w.table; d_phitsurerem=w.table; d_premsure=w.table; d_pknowsure=w.table;

% Calculate memory scores per cell (execute instructions)
disp('[Item memory] Executing cell instructions + calculating memory score ###################################')
for s=1:details.n_subjs
    disp(['Subject ' num2str(s) ':  ' details.subjects{s}])
    ws.d=subjdata{s,1}.mem;
    d_personality(s+1,:)=num2cell([subjdata{s,1}.p.NS subjdata{s,1}.p.RD]);
    for c=1:size(instruc.cells,1) ; % Each cell
        wc.d=ws.d(ws.d(:,col.m.EncCorrect)==1,:);        
        for x=1:size(instruc.cells{c,2},1)  % Compile cell by applying selection criteria serially
            eval(['wc.d=wc.d(wc.d(:,col.m.' instruc.cells{c,2}{x,1} ')=='  num2str(instruc.cells{c,2}{x,2}) ',:);'])
        end

        % For each memory measure
        subjdata{s,4}.fa_phitsure=0; subjdata{s,4}.fa_phitrem=0; subjdata{s,4}.fa_phitsurerem=0;subjdata{s,4}.fa_premsure=0;subjdata{s,4}.fa_pknowsure=0; % Fake FA rates to allow for execution in loop
        for m=1:length(instruc.mem)
            if isempty(wc.d)==0
                eval(['wm.mem=f_calcmem(' num2str(m) ',wc.d(:,col.m.Roc), subjdata{s,4}.fa_' instruc.mem{m} ');'])
                eval(['d_' instruc.mem{m}  '{s+1, c+1}=wm.mem.rate;'])
            else
                eval(['d_' instruc.mem{m}  '{s+1, c+1}=[];'])
            end
            wm=[];
        end
        wc=[];
    end
    
    
    % Overall dprime
    d_overalldprime{s+1,1}=d_dprime{s+1,2};
    ws=[];
end

%% Memory scores calculated from enc trialstats (i.e. context-blocked)
% only done for d' so far. could do for the other memory scores. 

% Setup data tables (prefix c_)
w.table=cell(details.n_subjs+1, size(instruc.contextcells,1)+1);w.table{1,1}='Subject'; w.table(1,2:end)=instruc.contextcells(:,1); w.table(2:end,1)=details.subjects;
c_hit=w.table;c_surehit=w.table;

% Calculate context memory scores (execute instructions)
disp('[Context memory] Executing cell instructions + calculating memory score ###################################')
instruc.contextmem={'hit';'surehit'};
for s=1:details.n_subjs
    disp(['Subject ' num2str(s) ':  ' details.subjects{s}])
    ws.d=subjdata{s,1}.enc;
    
    for c=1:size(instruc.contextcells,1) ; % Each cell
        wc.d=ws.d;
        for x=1:size(instruc.contextcells{c,2},1)  % Compile cell by applying selection criteria serially
            eval(['wc.d=wc.d(wc.d(:,col.et.' instruc.contextcells{c,2}{x,1} ')=='  num2str(instruc.contextcells{c,2}{x,2}) ',:);'])
        end
        
        % For each memory measure (manual)
        if isempty(wc.d)==0
            c_hit{s+1,c+1}=num2str(nanmean(wc.d(:,col.et.ContextMemScoreHit)));
            c_surehit{s+1,c+1}=num2str(nanmean(wc.d(:,col.et.ContextMemScoreSurehit)));
        else
            c_hit{s+1,c+1}=[];
            c_surehit{s+1,c+1}=[];
        end
        
        wc=[];
    end
    
    
    % Overall dprime
    d_overalldprime{s+1,1}=d_dprime{s+1,2};
    ws=[];
end

%% Correlation of memory scores within vs between trials

for s=1:details.n_subjs
    ws.d=sortrows(subjdata{s,2},col.m.EncTrial);
    ws.d(ws.d(:, col.m.EncSession)==2, col.m.EncTrial)=ws.d(ws.d(:, col.m.EncSession)==2, col.m.EncTrial)   +   max(ws.d(ws.d(:, col.m.EncSession)==1, col.m.EncTrial));  % Adjust trial num for block num
    ws.d=sortrows(ws.d,col.m.EncTrial);
    [b dev stats]=glmfit(ws.d(:, [col.m.EncTrial]),   abs(2-ws.d(:, col.m.OldNew)), 'binomial');

    
    error this is probably not the best way to do it!!!!!!
    Probably what i need is a chi square, or to compute a different behavioural index.
    
b
stats.p
    
    ws.d(:, [col.m.EncTrial col.m.OldNew])
    
    
    
    
    
end
    






disp('All subjects included!')

%% Output

input('Print 2 text?');

% Print item memory scores
for m=1:length(instruc.mem)
    eval(['wm.dat=d_' instruc.mem{m} ';'])
    [wm.printok]=print2txt([where.where filesep '3 Analysis inputs'], ['(' date ') Flexmem ' instruc.mem{m}],[wm.dat(:,1) d_overalldprime d_personality wm.dat(:,2:end)]);
    disp([instruc.mem{m} ' ---------------']); disp(wm.printok)
    wm=[];
end

% 

% Print context memory scores
for m=1:length(instruc.contextmem)
    eval(['wm.dat=c_' instruc.contextmem{m} ';'])
    [wm.printok]=print2txt([where.where filesep '3 Analysis inputs'], ['(' date ') Flexmem c_' instruc.contextmem{m}], [wm.dat(:,1) d_overalldprime d_personality wm.dat(:,2:end)]);
    disp(['c_' instruc.contextmem{m} ' ---------------']); disp(wm.printok)
    wm=[];
end



