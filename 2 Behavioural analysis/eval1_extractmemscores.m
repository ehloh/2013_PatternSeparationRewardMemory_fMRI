% eval1_extractmemscores
clear all; close all hidden; clc

for o1=1:1 % Documentation  
% Instructions list (for execution)
%     Col 1: Name of measure
%     Col 2: Cue status (1=Sim, 2=Dis, 999=n/a)
%     Col 3: Context status (1=Reward, 0=Neutral)
%     Col 4: Item status (1=Reward, 0=Neutral)
%     Col 5: Seed status (n/a)
%     Col 6: # of cells to corresponding to this measure
%     Col 7: Cells (combined) corresponding to this measure
%     Col 8: [diff only] Cells to deduct, for this measure
%     Col 9: Is this measure (theoretically) valid?
%     Col 10: Order for printing

% This script produces txt files for reading in SPSS. It also produces
% (un-saved) cell arrays corresponding to each measure, corresponding to:
% - cell means only (prefix 'm_')
% - difference between means only (prefix 'd_')
% - combination of both means and difference scores (prefix 'r_')

end

where.where='I:\2 fMRI behaviour analysis';
% where.where='/Volumes/PENNYDISK/2 fMRI behaviour analysis';
where.scripts=[where.where filesep '1 Analysis scripts']; where.data=[where.where filesep '2 Behavioural Data'];
where.printallmemfile='I:\1 fMRI analysis\2 Behavioural data\Group behaviour files'; % In fMRI folder - for beta correlations

for o1=1:1 % Options & Setup for analysis
    w.minimum_ntrials=1;
    w.extravar=14;
    w.print2text=1;
    
    % Compile settings
    load([where.scripts filesep 'i_eval1_mem.mat'])
    cd(where.data); w.dir=dir('p*'); for i=1:size(w.dir,1); details.subjects{i,1}=w.dir(i).name; end
    details.n_subjs=size(w.dir,1);
    addpath(where.scripts)
end

%% Load all data + Mark ROC (Col 31)

% 'subjdata': Col 1=Full trialstats, Col 2=Trialstats, old + correctly judged
%                   during encoding, Col 3=Trialstats for foils, Col
%                   4=Raw cell-splits, Col 5=Reorganized cells (design)
subjdata=cell(details.n_subjs,3);
for s=1:details.n_subjs
    ws.d=load([where.data filesep details.subjects{s} filesep details.subjects{s} '_file_3memorytest.mat']);
    ws.dat=ws.d.memtest.trialstats;
    % Mark ROC type (Col 31)
    %             Response types
    %                 1=Sure New
    %                 2= Unsure New
    %                 3= U K 
    %                 4= S K
    %                 5= U R
    %                 6= S R
    ws.dat(ws.dat(:,14)==1 & ws.dat(:,16)==1,31)=6; % Rem
    ws.dat(ws.dat(:,14)==1 & ws.dat(:,16)==2,31)=5;
    ws.dat(ws.dat(:,14)==2 & ws.dat(:,16)==1,31)=4; % Know
    ws.dat(ws.dat(:,14)==2 & ws.dat(:,16)==2,31)=3;
    ws.dat(ws.dat(:,12)==2 & ws.dat(:,16)==2,31)=2; % New
    ws.dat(ws.dat(:,12)==2 & ws.dat(:,16)==1,31)=1;
    % 
    subjdata{s,1}=ws.dat;
    subjdata{s,2}=ws.dat(ws.dat(:,3)==1 & ws.dat(:,29)==1, :);
    subjdata{s,3}=ws.dat(ws.dat(:,3)==2,:);
    %
    ws=[];
end

%% Pre-process for design

% Split data
for s=1:details.n_subjs
    subjdata{s,4}.s1=subjdata{s,2}(subjdata{s,2}(:,10)==1 & subjdata{s,2}(:,11)==1,:); % Sim R (sample 1)
    subjdata{s,4}.s3=subjdata{s,2}(subjdata{s,2}(:,10)==1 & subjdata{s,2}(:,11)==0,:); % Sim N (sample 3)
    subjdata{s,4}.s5=subjdata{s,2}(subjdata{s,2}(:,10)==2 & subjdata{s,2}(:,11)==1,:); % Dis R (sample 5)
    subjdata{s,4}.s7=subjdata{s,2}(subjdata{s,2}(:,10)==2 & subjdata{s,2}(:,11)==0,:); % Dis N (sample 7)
    for i=1:4  % Split according to actual outcome (0,1,2,3)
        eval(['subjdata{s,4}.s1_' num2str(i-1) '=subjdata{s,4}.s1(subjdata{s,4}.s1(:,5)==' num2str(i-1) ',:);   '])
        eval(['subjdata{s,4}.s3_' num2str(i-1) '=subjdata{s,4}.s3(subjdata{s,4}.s3(:,5)==' num2str(i-1) ',:);   '])
        eval(['subjdata{s,4}.s5_' num2str(i-1) '=subjdata{s,4}.s5(subjdata{s,4}.s5(:,5)==' num2str(i-1) ',:);   '])
        eval(['subjdata{s,4}.s7_' num2str(i-1) '=subjdata{s,4}.s7(subjdata{s,4}.s7(:,5)==' num2str(i-1) ',:);   '])
    end
end

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

%% Calculate memory scores (execute instructions)

% Mean scores
ws.temp=zeros(details.n_subjs,size(instruc,1)); r_dprime=ws.temp; r_hitrate=ws.temp; r_surehitrate=ws.temp; r_rem=ws.temp; r_surerem=ws.temp; r_know=ws.temp; r_sureknow=ws.temp;
r_p_hitsure=ws.temp;r_p_hitrem=ws.temp; r_p_hitsurerem=ws.temp; r_p_remsure=ws.temp; r_p_knowsure=ws.temp;
for s=1:details.n_subjs;
    for i=1:size(instruc,1)
        ws=subjdata{s,4};
        % Compile sample
        ws.s0=[];ws.s2=[];ws.s4=[];ws.s6=[];ws.s8=[]; % Dummies
        eval(['ws.sample=vertcat(ws.s' num2str(instruc{i,7}(1)) ', ws.s' num2str(instruc{i,7}(2)) ', ws.s' num2str(instruc{i,7}(3)) ', ws.s' num2str(instruc{i,7}(4)) ', ws.s' num2str(instruc{i,7}(5)) ', ws.s' num2str(instruc{i,7}(6)) ', ws.s' num2str(instruc{i,7}(7)) ', ws.s' num2str(instruc{i,7}(8))  ');'])
        if size(ws.sample,1)>0
        % Calculate + log (matlab vectors, overall matrix)
            % Dprime
            ws.res=f_calcmem(1,ws.sample(:,31), subjdata{s,4}.fa_dprime);
            eval(['d_dprime.' instruc{i,1} '(s,1)  =ws.res.rate;'])
            r_dprime(s,i)=ws.res.rate;
            % Hit rate
            ws.res=f_calcmem(6,ws.sample(:,31), subjdata{s,4}.fa_hitrate);
            eval(['d_hitrate.' instruc{i,1} '(s,1)  =ws.res.rate;'])
            r_hitrate(s,i)=ws.res.rate;
            % Sure Hitrate
            ws.res=f_calcmem(7,ws.sample(:,31), subjdata{s,4}.fa_surehitrate);
            eval(['d_surehitrate.' instruc{i,1} '(s,1)  =ws.res.rate;'])
            r_surehitrate(s,i)=ws.res.rate;
            % Remember
            ws.res=f_calcmem(2,ws.sample(:,31), subjdata{s,4}.fa_rem);
            eval(['d_rem.' instruc{i,1} '(s,1)  =ws.res.rate;'])
            r_rem(s,i)=ws.res.rate;
            % Sure Remember
            ws.res=f_calcmem(4,ws.sample(:,31), subjdata{s,4}.fa_surerem);
            eval(['d_surerem.' instruc{i,1} '(s,1)  =ws.res.rate;'])
            r_surerem(s,i)=ws.res.rate;
            % Know
            ws.res=f_calcmem(3,ws.sample(:,31), subjdata{s,4}.fa_know);
            eval(['d_know.' instruc{i,1} '(s,1)  =ws.res.rate;'])
            r_know(s,i)=ws.res.rate;
            % Sure Know
            ws.res=f_calcmem(5,ws.sample(:,31), subjdata{s,4}.fa_sureknow);
            eval(['d_sureknow.' instruc{i,1} '(s,1)  =ws.res.rate;'])
            r_sureknow(s,i)=ws.res.rate;  
            % Percentages ------------------------------------------
            % P(Sure|Hit)
            ws.res=f_calcmem(8,ws.sample(:,31), 0);
            eval(['p_hitsure.' instruc{i,1} '(s,1)=ws.res.rate;']);
            r_p_hitsure(s,i)=ws.res.rate;  
            % P(Rem|Hit)
            ws.res=f_calcmem(9,ws.sample(:,31), 0);
            eval(['p_hitrem.' instruc{i,1} '(s,1)=ws.res.rate;']);
            r_p_hitrem(s,i)=ws.res.rate;  
            % P(SureRem|Hit)
            ws.res=f_calcmem(10,ws.sample(:,31), 0);
            eval(['p_hitsurerem.' instruc{i,1} '(s,1)=ws.res.rate;']);
            r_p_hitsurerem(s,i)=ws.res.rate;
            % P(Sure|Rem)
            ws.res=f_calcmem(11,ws.sample(:,31), 0);
            eval(['p_remsure.' instruc{i,1} '(s,1)=ws.res.rate;']);
            r_p_remsure(s,i)=ws.res.rate;
            % P(Sure|Know)
            ws.res=f_calcmem(12,ws.sample(:,31), 0);
            eval(['p_knowsure.' instruc{i,1} '(s,1)=ws.res.rate;']);
            r_p_knowsure(s,i)=ws.res.rate;
        else
            eval(['d_dprime.' instruc{i,1} '(s,1)  =0;'])
            eval(['d_hitrate.' instruc{i,1} '(s,1)  =0;'])
            eval(['d_surehitrate.' instruc{i,1} '(s,1)  =0;'])
            eval(['d_rem.' instruc{i,1} '(s,1)  =0;'])
            eval(['d_surerem.' instruc{i,1} '(s,1)  =0;'])
            eval(['d_know.' instruc{i,1} '(s,1)  =0;'])
            eval(['d_sureknow.' instruc{i,1} '(s,1)  =0;'])
            eval(['p_hitsure.' instruc{i,1} '(s,1)=0;']); % Percentages
            eval(['p_hitrem.' instruc{i,1} '(s,1)=0;']); 
            eval(['p_hitsurerem.' instruc{i,1} '(s,1)=0;']);
            eval(['p_remsure.' instruc{i,1} '(s,1)=0;']);
            eval(['p_knowsure.' instruc{i,1} '(s,1)=0;']);
        end
        %
        ws=[];
    end
end

% Difference scores
w.i=size(instruc,1);
for s=1:details.n_subjs
    for i=1:size(instruc_diff,1)
        ws=subjdata{s,4};
        % Compile sample
        ws.s0=[];ws.s2=[];ws.s4=[];ws.s6=[];ws.s8=[]; % Dummies
        eval(['ws.sampleP=vertcat(ws.s' num2str(instruc_diff{i,7}(1)) ', ws.s' num2str(instruc_diff{i,7}(2)) ', ws.s' num2str(instruc_diff{i,7}(3)) ', ws.s' num2str(instruc_diff{i,7}(4)) ', ws.s' num2str(instruc_diff{i,7}(5)) ', ws.s' num2str(instruc_diff{i,7}(6)) ', ws.s' num2str(instruc_diff{i,7}(7)) ', ws.s' num2str(instruc_diff{i,7}(8))  ');'])
        eval(['ws.sampleM=vertcat(ws.s' num2str(instruc_diff{i,8}(1)) ', ws.s' num2str(instruc_diff{i,8}(2)) ', ws.s' num2str(instruc_diff{i,8}(3)) ', ws.s' num2str(instruc_diff{i,8}(4)) ', ws.s' num2str(instruc_diff{i,8}(5)) ', ws.s' num2str(instruc_diff{i,8}(6)) ', ws.s' num2str(instruc_diff{i,8}(7)) ', ws.s' num2str(instruc_diff{i,8}(8))  ');'])
        if size(ws.sampleM,1)>0 && size(ws.sampleP,1)>0
            % Calculate + log (matlab vectors, overall matrix)
            % Dprime
            ws.resP=f_calcmem(1,ws.sampleP(:,31), subjdata{s,4}.fa_dprime);
            ws.resM=f_calcmem(1,ws.sampleM(:,31), subjdata{s,4}.fa_dprime);
            ws.res.rate=ws.resP.rate-ws.resM.rate;
            eval(['d_dprime.' instruc_diff{i,1} '(s,1)  =ws.res.rate;'])
            r_dprime(s,i+w.i)=ws.res.rate;
            % Hit rate
            ws.resP=f_calcmem(6,ws.sampleP(:,31), subjdata{s,4}.fa_hitrate);
            ws.resM=f_calcmem(6,ws.sampleM(:,31), subjdata{s,4}.fa_hitrate);
            ws.res.rate=ws.resP.rate-ws.resM.rate;
            eval(['d_hitrate.' instruc_diff{i,1} '(s,1)  =ws.res.rate;'])
            r_hitrate(s,i+w.i)=ws.res.rate;
            % Sure Hitrate
            ws.resP=f_calcmem(7,ws.sampleP(:,31), subjdata{s,4}.fa_surehitrate);
            ws.resM=f_calcmem(7,ws.sampleM(:,31), subjdata{s,4}.fa_surehitrate);
            ws.res.rate=ws.resP.rate-ws.resM.rate;
            eval(['d_surehitrate.' instruc_diff{i,1} '(s,1)  =ws.res.rate;'])
            r_surehitrate(s,i+w.i)=ws.res.rate;
            % Remember
            ws.resP=f_calcmem(2,ws.sampleP(:,31), subjdata{s,4}.fa_rem);
            ws.resM=f_calcmem(2,ws.sampleM(:,31), subjdata{s,4}.fa_rem);
            ws.res.rate=ws.resP.rate-ws.resM.rate;
            eval(['d_rem.' instruc_diff{i,1} '(s,1)  =ws.res.rate;'])
            r_rem(s,i+w.i)=ws.res.rate;
            % Sure Remember
            ws.resP=f_calcmem(4,ws.sampleP(:,31), subjdata{s,4}.fa_surerem);
            ws.resM=f_calcmem(4,ws.sampleM(:,31), subjdata{s,4}.fa_surerem);
            ws.res.rate=ws.resP.rate-ws.resM.rate;
            eval(['d_surerem.' instruc_diff{i,1} '(s,1)  =ws.res.rate;'])
            r_surerem(s,i+w.i)=ws.res.rate;
            % Know
            ws.resP=f_calcmem(3,ws.sampleP(:,31), subjdata{s,4}.fa_know);
            ws.resM=f_calcmem(3,ws.sampleM(:,31), subjdata{s,4}.fa_know);
            ws.res.rate=ws.resP.rate-ws.resM.rate;
            eval(['d_know.' instruc_diff{i,1} '(s,1)  =ws.res.rate;'])
            r_know(s,i+w.i)=ws.res.rate;
            % Sure Know
            ws.resP=f_calcmem(5,ws.sampleP(:,31), subjdata{s,4}.fa_sureknow);
            ws.resM=f_calcmem(5,ws.sampleM(:,31), subjdata{s,4}.fa_sureknow);
            ws.res.rate=ws.resP.rate-ws.resM.rate;
            eval(['d_sureknow.' instruc_diff{i,1} '(s,1)  =ws.res.rate;'])
            r_sureknow(s,i+w.i)=ws.res.rate;
            % Percentages ------------------------------------------
            % P(Sure|Hit)
            ws.resP=f_calcmem(8,ws.sampleP(:,31), 0);
            ws.resM=f_calcmem(8,ws.sampleM(:,31), 0);
            ws.res.rate=ws.resP.rate-ws.resM.rate;
            eval(['p_hitsure.' instruc_diff{i,1} '(s,1)=ws.res.rate;']);
            r_p_hitsure(s,i+w.i)=ws.res.rate;  
            % P(Rem|Hit)
            ws.resP=f_calcmem(9,ws.sampleP(:,31), 0);
            ws.resM=f_calcmem(9,ws.sampleM(:,31), 0);
            ws.res.rate=ws.resP.rate-ws.resM.rate;
            eval(['p_hitrem.' instruc_diff{i,1} '(s,1)=ws.res.rate;']);
            r_p_hitrem(s,i+w.i)=ws.res.rate;  
            % P(SureRem|Hit)
            ws.resP=f_calcmem(10,ws.sampleP(:,31), 0);
            ws.resM=f_calcmem(10,ws.sampleM(:,31), 0);
            ws.res.rate=ws.resP.rate-ws.resM.rate;
            eval(['p_hitsurerem.' instruc_diff{i,1} '(s,1)=ws.res.rate;']);
            r_p_hitsurerem(s,i+w.i)=ws.res.rate;
            % P(Sure|Rem)
            ws.resP=f_calcmem(11,ws.sampleP(:,31), 0);
            ws.resM=f_calcmem(11,ws.sampleM(:,31), 0);
            ws.res.rate=ws.resP.rate-ws.resM.rate;
            eval(['p_remsure.' instruc_diff{i,1} '(s,1)=ws.res.rate;']);
            r_p_remsure(s,i+w.i)=ws.res.rate;
            % P(Sure|Know)
            ws.resP=f_calcmem(12,ws.sampleP(:,31), 0);
            ws.resM=f_calcmem(12,ws.sampleM(:,31), 0);
            ws.res.rate=ws.resP.rate-ws.resM.rate;
            eval(['p_knowsure.' instruc_diff{i,1} '(s,1)=ws.res.rate;']);
            r_p_knowsure(s,i+w.i)=ws.res.rate;
        else
            eval(['d_dprime.' instruc_diff{i,1} '(s,1)  =999;'])
            eval(['d_hitrate.' instruc_diff{i,1} '(s,1)  =999;'])
            eval(['d_surehitrate.' instruc_diff{i,1} '(s,1)  =999;'])
            eval(['d_rem.' instruc_diff{i,1} '(s,1)  =999;'])
            eval(['d_surerem.' instruc_diff{i,1} '(s,1)  =999;'])
            eval(['d_know.' instruc_diff{i,1} '(s,1)  =999;'])
            eval(['d_sureknow.' instruc_diff{i,1} '(s,1)  =999;'])
            eval(['p_hitsure.' instruc_diff{i,1} '(s,1)=999;']); % Percentages
            eval(['p_hitrem.' instruc_diff{i,1} '(s,1)=999;']); 
            eval(['p_hitsurerem.' instruc_diff{i,1} '(s,1)=999;']);
            eval(['p_remsure.' instruc_diff{i,1} '(s,1)=999;']);
            eval(['p_knowsure.' instruc_diff{i,1} '(s,1)=999;']);
        end
        %
        ws=[];
    end
end

%% Subject information

% error('Pause')

% Formatting to match statistical analysis script)
raw.subjectinfo=cell(details.n_subjs, 12); raw.subjectinfo{1,1}='Expt_version'; raw.subjectinfo{1,2}='Counterbal';raw.subjectinfo{1,3}='Subject';raw.subjectinfo{1,4}='p.NSall'; raw.subjectinfo{1,5}='p.RDall'; raw.subjectinfo{1,6}='p.NS_1'; raw.subjectinfo{1,7}='p.NS_2';raw.subjectinfo{1,8}='p.NS_3'; raw.subjectinfo{1,9}='p.NS_4';raw.subjectinfo{1,10}='p.RD_1';raw.subjectinfo{1,11}='p.RD_2';raw.subjectinfo{1,12}='p.RD_3';raw.subjectinfo{1,13}='p.RD_4';raw.subjectinfo{1,14}='dprime';
for s=1:details.n_subjs
    raw.subjectinfo{s+1,1}='2.5 fMRI';
    ws=load([where.data filesep details.subjects{s}  filesep details.subjects{s} '_file_3memorytest.mat']);
    try
        raw.subjectinfo{s+1,2}=ws.memtest.settings.counterbal_ver;
    catch
        raw.subjectinfo{s+1,2}=999;
    end
    raw.subjectinfo{s+1,3}=ws.memtest.settings.subjectlog.name;
    raw.subjectinfo{s+1,4}=mean(ws.memtest.TPQ.NS(1:4));
    raw.subjectinfo{s+1,5}=mean(ws.memtest.TPQ.RD(1:4));
    raw.subjectinfo{s+1,6}=ws.memtest.TPQ.NS(1);
    raw.subjectinfo{s+1,7}=ws.memtest.TPQ.NS(2);
    raw.subjectinfo{s+1,8}=ws.memtest.TPQ.NS(3);
    raw.subjectinfo{s+1,9}=ws.memtest.TPQ.NS(4);
    raw.subjectinfo{s+1,10}=ws.memtest.TPQ.RD(1);
    raw.subjectinfo{s+1,11}=ws.memtest.TPQ.RD(2);
    raw.subjectinfo{s+1,12}=ws.memtest.TPQ.RD(3);
    raw.subjectinfo{s+1,13}=ws.memtest.TPQ.RD(4);
    raw.subjectinfo{s+1,14}=d_dprime.Overall.Overall(s);
    %
    ws=[];
end
for i=1:w.extravar
    for s = 1:details.n_subjs
        if strcmp(raw.subjectinfo{1,i}, 'Expt_version')==1 || strcmp(raw.subjectinfo{1,i}, 'Counterbal')==1 || strcmp(raw.subjectinfo{1,i}, 'Subject')==1 
            eval(['sub_info.' raw.subjectinfo{1,i} '{s,1}=raw.subjectinfo{s+1,i};'])
        else
            eval(['sub_info.' raw.subjectinfo{1,i} '(s,1)=raw.subjectinfo{s+1,i};'])
        end
    end
end
sub_info.NS_1=sub_info.p.NS_1;sub_info.NS_2=sub_info.p.NS_2;sub_info.NS_3=sub_info.p.NS_3;sub_info.NS_4=sub_info.p.NS_4;
sub_info.RD_1=sub_info.p.RD_1;sub_info.RD_2=sub_info.p.RD_2;sub_info.RD_3=sub_info.p.RD_3;sub_info.RD_4=sub_info.p.RD_4;
sub_info.NSall=mean([sub_info.NS_1,sub_info.NS_2,sub_info.NS_3,sub_info.NS_4],2);
sub_info.RDall=mean([sub_info.RD_1,sub_info.RD_2,sub_info.RD_3,sub_info.RD_4],2);

%% Formatting

% Format for matlab statistical analysis script
raw.dprime=r_dprime; raw.hitrate=r_hitrate; raw.surehitrate=r_surehitrate; raw.rem=r_rem;  raw.surerem=r_surerem; raw.sureknow=r_sureknow;raw.know=r_know;
raw.p_hitrem=r_p_hitrem; raw.p_hitsure=r_p_hitsure;raw.p_hitsurerem=r_p_hitsurerem;raw.p_knowsure=r_p_knowsure; raw.p_remsure=r_p_remsure;
raw.headers=cell(1,size(r_dprime,2)); % Headers
for i=1:size(instruc,1)
    raw.headers{i}=instruc{i,1};
end
for i=1:size(instruc_diff,1)
    raw.headers{i+w.i}=instruc_diff{i,1};
end
save([where.where filesep '3 Analysis inputs' filesep '(' date ') Memory data.mat'], 'raw', 'sub_info', 'd_dprime', 'd_hitrate' , 'd_know' , 'd_rem' , 'd_surehitrate' , 'd_sureknow', 'd_surerem', 'p_hitrem', 'p_hitsure', 'p_hitsurerem', 'p_remsure', 'p_knowsure')

% Print to txt files
if w.print2text==1
    
    dat2print={'dprime' 'md';
                     'hitrate'  'mh';
                     'surehitrate'  'msh'; 
                     'rem'  'mr'; 
                     'know'  'mk'; 
                     'surerem'  'msr'; 
                     'sureknow'  'msk'; 
                      'p_hitrem'  'mphr'; 
                      'p_hitsure'  'mphs'; 
                      'p_hitsurerem'  'mphsr'; 
                      'p_knowsure'  'mpks'; 
                      'p_remsure'  'mprs'; 
                     }; % memtype & prefix for memtype in the overall document (below)
                     
    % Headers for individual-memtype documents
    headers=cell(1,size(instruc,1)+size(instruc_diff,1)); w.h=size(raw.subjectinfo,2); 
    for i=1:size(instruc,1); headers{i}=instruc{i,1}; end
    for i=1:size(instruc_diff,1); headers{i+w.i}=instruc_diff{i,1}; end
    
    % Automated print each memtype separately
    disp(' ---------------------------------------------------------------------'); disp('PRINTING data files to txt')
    memtypes=cell(size(dat2print,1) ,2);
    for i=1:size(dat2print,1) 
        eval(['ws.d=r_' dat2print{i,1} ';'])
        ws.t=raw.subjectinfo;
        for j=1:size(instruc,1)+size(instruc_diff,1)
            ws.t{1,w.h+j}=headers{j};
            for s=1:details.n_subjs
                ws.t{s+1,w.h+j}=ws.d(s,j);
            end
        end
        printok=print2txt([where.where filesep '3 Analysis inputs'], ['(' date ') ' num2str(i) ' ' dat2print{i,1}], ws.t);
        disp(['Printing:     ' dat2print{i,1}])
        disp(printok)
        %
        memtypes{i,1}=dat2print{i,1};
        memtypes{i,2}=ws.t;
    end
    
    % Print overall memory document (containing all memtypes - for correlation script)
    d_allmem=raw.subjectinfo(:, [3 2 4:13]); % Selected certain subject details
    for i=1:size(dat2print,1)        
    	wa.t=memtypes{i,2}(:,size(raw.subjectinfo,2)+1:end);

        %
        wa.titles=cellfun(@(x)[dat2print{i,2} '.' x], wa.t(1,:), 'UniformOutput', 0); % Alter titles
        wa.t(1,:)=wa.titles;
        
        %
        d_allmem=[d_allmem, wa.t];
        wa=[];
    end
    printok=print2txt(where.printallmemfile, ['(' date ') Memory scores all'], d_allmem);
    disp('Printing master memory file in fMRI folder:'); disp(printok)
    
    %
    disp('DONE with printing')
    disp(' ---------------------------------------------------------------------')
end
