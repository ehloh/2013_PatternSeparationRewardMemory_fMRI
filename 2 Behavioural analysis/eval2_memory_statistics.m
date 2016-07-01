% Statistical analysis script for behavioural data
clear all;close all hidden;clc
where.where='D:\Dropbox\SCRIPPS\2a ContextMem behaviour';
% where.where='/Users/EleanorL/Dropbox/sandisk/2a fMRI behaviour analysis';

% Memory measure
memtype={'dprime', 'hitrate', 'surehitrate', 'rem',  'know', 'surerem','sureknow', 'p_hitrem', 'p_hitsure', 'p_hitsurerem', 'p_knowsure', 'p_remsure'};
for i=1:length(memtype); disp([num2str(i) '  =  ' memtype{i}]);  end; w.whichmeasure=input('Which measure?    ');
alldata=load([where.where filesep '3 Analysis inputs' filesep '(29-May-2013) Memory data.mat']);
    
for o1=1:1 % Settings for analysis
    w.measure=memtype{w.whichmeasure}; disp(['Measure=' w.measure])

    % Requested results -----------------------------
    w.boxplot=0;
    w.mainstatisticalanalysis=1;
    w.correlations=1;
    mediansplit.switch=0; % Median split on or off?
    
    % Exclude for overall performance -----------------
    exclude.noexclusions=0; % If set to 1, all exclusions are turned off
    exclude.overalldprime=1; % Exclude by overall poor dprime
    exclude.overallmem=0; % Exclude by target mem measure
    % Settings (high & low cut offs)
    exclude.overalldprime_cutoff_low=0.3; % dprime
    exclude.overalldprime_cutoff_high=10.0;
    exclude.overallmem_cutoff_low=0.00; % target memory measure
    exclude.overallmem_cutoff_high=10.0;
    
    % Median splits: execute by running 2 exclusions ------------
    if mediansplit.switch==1
        disp('Options: (a) dprime (b) rem (c) know (d)surerem (e)sureknow (f)hitrate')
        disp('              (g) NSall (h) RDall (i) NS subscale (j) RD subscale')
        w.splitmesnum=' ';
        w.splitmeasure=input('Split by: ','s');
        switch w.splitmeasure
            case 'a'
                mediansplit.measure='dprime';
            case 'b'
                mediansplit.measure='rem';
            case 'c'
                mediansplit.measure='know';
            case 'd'
                mediansplit.measure='surerem';
            case 'e'
                mediansplit.measure='sureknow';
            case 'f'
                mediansplit.measure='hitrate';
            case 'g'
                mediansplit.measure='NSall';
            case 'h'
                mediansplit.measure='RDall';
            case 'i'
                w.splitmesnum=input('NS Subscale #: ');
                mediansplit.measure=['NS_' num2str(w.splitmesnum)];
            case 'j'
                w.splitmesnum=input('RD Subscale #: ');
                mediansplit.measure=['RD_' num2str(w.splitmesnum)];
        end
%         mediansplit.measure='NSall'; % Name of memory measure, or TPQ measure
        mediansplit.half=input('Median split, which half? (1=Higher, 2=Lower): ');
    disp([w.splitmeasure ' ' num2str(w.splitmesnum) '(' num2str(mediansplit.half)  ')']);
    end
    
    % Load details
    load([where.where filesep '1 Analysis scripts' filesep 'i_correlations.mat']);
    if w.whichmeasure<8
        eval(['data=alldata.d_' w.measure ';'])
    else
        eval(['data=alldata.' w.measure ';'])
    end
    
end 

%% Implement any restrictions (subject exclusions etc)
    
for o1=1:1  % Create exclusion lists 'w.exclude', and then Inclusion lists    
    w.exclude=zeros(size(alldata.sub_info.Subject,1),1);
    
    % Exclusion based on extra-data details (e.g. misinterpreting memtest)
    w.extraexcl={}; % e.g. 'p21_SH'; ..
    settings.allsubjects=alldata.sub_info.Subject;
    for i=1:size(settings.allsubjects,1)
        settings.allsubjects{i,2}=i;
    end
    if isempty(w.extraexcl)==0 && exclude.noexclusions==0
        for i=1:size(w.extraexcl,1)
            ws.found=0;
            ws.counter=1;
            while ws.found==0
                if strcmp(cellstr(w.extraexcl{i}),cellstr(settings.allsubjects{ws.counter,1}))==1
                    w.extraexclnum(i)=ws.counter;
                    ws.found=1;
                else
                    if ws.counter==size(settings.allsubjects,1)
                        disp('Error: Could not find subject marked for exclusion (for exrta-data details)')
                        ws.found=1;
                    else
                        ws.counter=ws.counter+1;
                    end
                end
            end
        end
        w.exclude(w.extraexclnum,1)=1;
    end
    
    % Exclusion based on performance
    for i=1:size(alldata.sub_info.Subject,1)
        if exclude.overalldprime==1 % Dprime exclusion
            if alldata.d_dprime.Overall.Overall(i,1)<exclude.overalldprime_cutoff_low || alldata.d_dprime.Overall.Overall(i,1)>exclude.overalldprime_cutoff_high
                w.exclude(i)=1;
            end
        end
        if exclude.overallmem==1 % Target memory score exclusion
            ws=[];
            eval(['ws.score=alldata.d_' w.measure '.Overall.Overall;'])
            if ws.score(i,1)<exclude.overallmem_cutoff_low || ws.score(i,1)>exclude.overallmem_cutoff_high
                w.exclude(i)=1;
            end
        end
            % <----- To exlude based on additional criteria, insert here
    end

    if exclude.noexclusions==1 % Cancel all exclusions?
        w.exclude=zeros(size(alldata.sub_info.Subject,1),1);
    end
end
for o1=1:1 % Median splits?
if mediansplit.switch==1
    if sum(strcmp(mediansplit.measure, {'dprime';'rem';'know';'surerem';'sureknow';'hitrate'}))>0
        ws.mediansplit.measure=['d_' mediansplit.measure '.Overall.Overall'];
    else
        ws.mediansplit.measure=['sub_info.' mediansplit.measure];
    end
    eval(['ws.measure_allsubj=alldata.'  ws.mediansplit.measure ';'])
    ws.counter=1; % Incorporate prior exclusions into median split
    for i=1:size(ws.measure_allsubj)
        if w.exclude(i)==0
            ws.measure(ws.counter,1)=ws.measure_allsubj(i);
            ws.counter=ws.counter+1;
        end
    end
    ws.measure(:,2)=1:size(ws.measure,1);
    ws.measure=sortrows(ws.measure,1);
    median.subjs_lower=ws.measure(1:ceil(size(ws.measure,1)/2),2); % if odd #, lower group has >
    median.subjs_higher=ws.measure(ceil(size(ws.measure,1)/2)+1:size(ws.measure,1),2);
    switch mediansplit.half
        case 1 % Higher scores
            ws.exclude=median.subjs_lower;
        case 2 % Lower scores
            ws.exclude=median.subjs_higher;
    end
    w.exclude(ws.exclude,1)=1;
    ws=[];
end    
end
   
% Re-read data, including only the selected subjects
w.include=1-w.exclude;
dat=data;
data=[];
w.head=alldata.raw.headers';
for j=1:size(w.head,1)
    if isnan(alldata.raw.dprime(1,j))==0
        ws.measure=w.head{j};
        eval(['ws.alldat=dat.' ws.measure ';']) % Read the correct measure
        ws.count=1;
        for i=1:size(w.include);
            if w.include(i)==1
                ws.dat(ws.count,1)=ws.alldat(i);
                ws.count=ws.count+1;
            end
        end
        eval(['data.' ws.measure '=ws.dat;']) % Return to read data
        ws=[];
    else
    end
end

%% Statistical analysis
w.n_subjs=size(data.Overall.Overall,1);

if w.mainstatisticalanalysis==1
    
% Factorial ANOVA
anova.data=horzcat(data.Cell.Sim_R,data.Cell.Sim_N,data.Cell.Dis_R,data.Cell.Dis_N);
anova.levels=[2 2]; % How many levels in each factor?
anova.labels={'Similarity';'Valence'};
[anova.res]=teg_repeated_measures_ANOVA(anova.data,anova.levels,anova.labels);
for i=1:3 % Mark ANOVA significance
    if anova.res.R(i,4)<0.05
        anova.res.p{i}=[num2str(anova.res.R(i,4)) '    *  '];
    elseif anova.res.R(i,4)<0.1
        anova.res.p{i}=[num2str(anova.res.R(i,4)) '    (t)'];
    else
        anova.res.p{i}=[num2str(anova.res.R(i,4)) '       '];
    end
end

% Simple effects (t-tests for graphs)
[graphs.diff_ttests(1,1), graphs.diff_ttests(1,2)]=ttest(data.Summ.Rew,data.Summ.Neu); % TTESTS: Context effects start
[graphs.diff_ttests(2,1), graphs.diff_ttests(2,2)]=ttest(data.Cell.Sim_R,data.Cell.Sim_N);
[graphs.diff_ttests(3,1), graphs.diff_ttests(3,2)]=ttest(data.Cell.Dis_R,data.Cell.Dis_N);
[graphs.diff_ttests(5,1), graphs.diff_ttests(5,2)]=ttest(data.Summ.Sim,data.Summ.Dis);  % Cue effects start
[graphs.diff_ttests(6,1), graphs.diff_ttests(6,2)]=ttest(data.Cell.Sim_R,data.Cell.Dis_R);
[graphs.diff_ttests(7,1), graphs.diff_ttests(7,2)]=ttest(data.Cell.Sim_N,data.Cell.Dis_N);
graphs.diff_ttests_p=cell(7,1);
for i=[1 2 3 5 6 7] % Mark trends
    if graphs.diff_ttests(i,2)<0.05
        graphs.diff_ttests_p{i}=[num2str(graphs.diff_ttests(i,2)) '    *  ' ];        
    elseif graphs.diff_ttests(i,2)<0.1
        graphs.diff_ttests_p{i}=[num2str(graphs.diff_ttests(i,2)) '    (t)'];
    else
        graphs.diff_ttests_p{i}=[num2str(graphs.diff_ttests(i,2)) '       '];
    end
end

for o1=1:1 % Generate values for graphs 

    % Group means (factorial)
    fig.labels_groupmeans=['Similar  ';'Disimilar'];
    graphs.factorial_means=[mean(data.Cell.Sim_R) mean(data.Cell.Sim_N); mean(data.Cell.Dis_R) mean(data.Cell.Dis_N)];
    graphs.factorial_errorbars=[std(data.Cell.Sim_R)/(sqrt(w.n_subjs)) std(data.Cell.Sim_N)/(sqrt(w.n_subjs)); std(data.Cell.Dis_R)/(sqrt(w.n_subjs)) std(data.Cell.Dis_N)/(sqrt(w.n_subjs))];
    graphs.factorial_errorbars95within=graphs.factorial_errorbars*1.96;

    %  Main-effect graphs
    fig.labels_maineffectValence=['Context_R'; 'Context_N'];
    graphs.maineffectValence_means(1)=mean(data.Summ.Rew);
    graphs.maineffectValence_means(2)=mean(data.Summ.Neu);
    graphs.maineffectValence_errorbars(1)=std(data.Summ.Rew)/(sqrt(w.n_subjs)); % Error bars = +/- SE
    graphs.maineffectValence_errorbars(2)=std(data.Summ.Neu)/(sqrt(w.n_subjs));
    graphs.maineffectValence_errorbars95within=graphs.maineffectValence_errorbars*1.96;
    fig.labels_maineffectSimilarity=['Sim'; 'Dis'];
    graphs.maineffectSimilarity_means(1)=mean(data.Summ.Sim);
    graphs.maineffectSimilarity_means(2)=mean(data.Summ.Dis);
    graphs.maineffectSimilarity_errorbars(1)=std(data.Summ.Sim)/(sqrt(w.n_subjs));
    graphs.maineffectSimilarity_errorbars(2)=std(data.Summ.Dis)/(sqrt(w.n_subjs));
    graphs.maineffectSimilarity_errorbars95within=graphs.maineffectSimilarity_errorbars*1.96;

    % Difference graphs
    fig.labels_differences={'Reward effect'; 'Similarity effect'};
    graphs.diff_means=[mean(data.Summ.Rew-data.Summ.Neu) mean(data.cell.Sim_valfx) mean(data.cell.Dis_valfx); mean(data.Summ.Sim-data.Summ.Dis) mean(data.cell.R_simfx) mean(data.cell.N_simfx)];
    graphs.diff_errorbars=[std(data.Summ.Rew-data.Summ.Neu)/sqrt(w.n_subjs) std(data.cell.Sim_valfx)/sqrt(w.n_subjs) std(data.cell.Dis_valfx)/sqrt(w.n_subjs); std(data.Summ.Sim-data.Summ.Dis)/sqrt(w.n_subjs) std(data.cell.R_simfx)/sqrt(w.n_subjs) std(data.cell.N_simfx)/sqrt(w.n_subjs)];
    graphs.diff_errorbars95within=graphs.diff_errorbars*1.96;
end


%% Plot graphs & Label

scrsz = get(0,'ScreenSize'); % scrsz=[x offset from left edge of screen, y offset from bottom of screen, horizontal width, vertical height];

for o2=1:1 % Boxplots
    if w.boxplot==1
    fig.boxplots=figure('Name',['[' w.measure '] Boxplots'],'NumberTitle','off','Position',[10, round(scrsz(4)/3*2)-90, round(scrsz(3)/3-10), round(scrsz(4)/3)]);
    subplot(4,4,[1 2 5 6])
    boxplot(data.Overall.Overall,'labels','Overall')
    subplot(4,4,[3 4 7 8])
    boxplot([data.Cell.Sim_R data.Cell.Sim_N data.Cell.Dis_R data.Cell.Dis_N],'labels', ['Sim_R'; 'Sim_N'; 'Dis_R'; 'Dis_N'], 'plotstyle','compact')
    subplot(4,4,[9 10 13 14])
    boxplot([data.Summ.Sim data.Summ.Dis],'labels',['Sim'; 'Dis'])
    subplot(4,4,[11 12 15 16])
    boxplot([data.Summ.Rew data.Summ.Neu],'labels',['context_R'; 'context_N'])
    end
end

% fig.maingraphs=figure('Name',['[' w.measure '] Group means (Factorial)'],'NumberTitle','off','Position',[10+(scrsz(3)/3), (scrsz(4)/3)-90, (scrsz(3)/3*2)*1000, (scrsz(4)/3*2)*1000]);
fig.maingraphs=figure('Name',['[' w.measure '] Group means (Factorial)'],'NumberTitle','off','Position',[715,220,1200,900]);
% fig.maingraphs=figure('Name',['[' w.measure '] Group means (Factorial)'],'NumberTitle','off','Position',[480,140,800,550]);
set(gcf,'Color',[1 1 1])

% Group means (factorial)
plotfac=subplot(5,5,[1 2 6 7]); % Plot
barwitherr(graphs.factorial_errorbars,graphs.factorial_means)
set(gca, 'XTickLabel', fig.labels_groupmeans); 
axis([0.5, 2.5, min(min(graphs.factorial_means))-max(max(graphs.factorial_errorbars))*2,max(max(graphs.factorial_means))+max(max(graphs.factorial_errorbars))*2])
ylabel(['Mean ' w.measure ' score'])
title('Factorial-means graph (error bars = +/-SE)')
legend(plotfac,['Reward ';'Neutral'])

% Main-effect graphs
subplot(5,5,[16 17 21 22]) % Valence
barwitherr([graphs.maineffectValence_errorbars 0 graphs.maineffectSimilarity_errorbars], [graphs.maineffectValence_means 0 graphs.maineffectSimilarity_means],'m')
set(gca, 'XTickLabel', ['cR '; 'cN '; '   '; 'Sim'; 'Dis'])
axis([0,6, min(min(graphs.factorial_means))-max(max(graphs.factorial_errorbars))*2,max(max(graphs.factorial_means))+max(max(graphs.factorial_errorbars))*2])
ylabel(['Mean ' w.measure ' score'])
xlabel('Main factor levels');

% Difference graphs
subplot(5,5, [4 5 9 10]) % Error bars = SE
barwitherr(graphs.diff_errorbars, graphs.diff_means)
set(gca, 'XTickLabel', fig.labels_differences)
ylabel(['Mean difference (' w.measure ')'])
title('Difference graphs (error bars = +/-SE)')
subplot(5,5, [19 20 24 25]) % Error bars = 95% CI
barwitherr(graphs.diff_errorbars95within, graphs.diff_means)
set(gca, 'XTickLabel', fig.labels_differences)
ylabel(['Mean difference (' w.measure ')'])
title('Difference graphs (error bars = 95% CI)')

for o1=1:1 % Print details to figure
% Footnotes
if exclude.overalldprime==1
    wp.inc1=['dprime >' num2str(exclude.overalldprime_cutoff_low) '; '];
else
    wp.inc1=' ';
end
if exclude.overallmem==1 && w.whichmeasure~=1
    wp.inc2=[w.measure ' >' num2str(exclude.overallmem_cutoff_low) ';'];
else
    wp.inc2=' ';
end
if exclude.noexclusions==1
    wp.inclusion=['All subjects included (N= ' num2str( w.n_subjs) ')'];    
else
    wp.inclusion=['Inclusion criteria: ' wp.inc1 wp.inc2 ' (N= ' num2str( w.n_subjs) ')'];    
end 
if  mediansplit.switch==1
    if mediansplit.half==1
        wp.medd='higher scores';
    else
        wp.medd='lower scores';
    end
    wp.med=['MEDIAN SPLIT by ' mediansplit.measure ' (' wp.medd ')'];
else
    wp.med=' ';
end
f(1)=text(0,0,['[' w.measure ']' wp.med]);
f(2)=text(0,0, wp.inclusion);
set(f(:),'Units','normalized')
set(f(1),'Position',[-1.8,2.95,0])
set(f(2),'Position',[-1.8,2.88,0])

% Statistical results
t(1)=text(0,0,' '); % ANOVA
t(2)=text(0,0, 'SIMILARITY X VALENCE ANOVA');
t(3)=text(0,0,  ['Similarity:     p=' anova.res.p{1}] );
t(4)=text(0,0,['Valence:       p=' anova.res.p{2}] );
t(5)=text(0,0, ['Sim x Val:     p=' anova.res.p{3}] );
line(1)=text(0,0, ' --------------------------------------------------------------------------------------------------');
line(2)=text(0,0,' --------------------------------------------------------------------------------------------------');
line(3)=text(0,0,' --------------------------------------------------------------------------------------------------');
line(4)=text(0,0,' --------------------------------------------------------------------------------------------------');
t(6)=text(0,0, 'REWARD EFFECTS'); % DIFFERENCE GRAPHS
t(7)=text(0,0,    ['Overall:     p=  ' graphs.diff_ttests_p{1}] );
t(8)=text(0,0,    ['Similar:     p=  ' graphs.diff_ttests_p{2}] );
t(9)=text(0,0,    ['Dissimilar: p=  ' graphs.diff_ttests_p{3}] );
t(10)=text(0,0, 'SIMILARITY EFFECTS');
t(11)=text(0,0,    ['Overall:  p=  ' graphs.diff_ttests_p{5}] );
t(12)=text(0,0,    ['Reward: p=  ' graphs.diff_ttests_p{6}] );
t(13)=text(0,0,    ['Neutral:  p=  ' graphs.diff_ttests_p{7}] );

% Text positions -- 
set(t(:),'Units','normalized')
set(line(:),'Units','normalized')
% Factorial
set(t(1),'Position',[-1.3,1.60,0])
set(t(2),'Position',[-1.55,1.49,0])
set(t(3),'Position',[-1.43,1.38,0])
set(t(4),'Position',[-1.43,1.29,0])
set(t(5),'Position',[-1.43,1.20,0])
set(line(1),'Position',[-1.7,1.58,0]) % Lines
set(line(2),'Position',[-1.7,1.17,0])
% Difference
set(t(6),'Position',[-0.1,1.49,0])
set(t(7),'Position',[-0.1,1.38,0])
set(t(8),'Position',[-0.1,1.29,0])
set(t(9),'Position',[-0.1,1.2,0])
set(t(10),'Position',[0.59,1.49,0]) % Similarity
set(t(11),'Position',[0.59,1.38,0])
set(t(12),'Position',[0.59,1.29,0])
set(t(13),'Position',[0.59,1.2,0])
set(line(3),'Position',[-0.15,1.58,0]) % Lines
set(line(4),'Position',[-0.15, 1.17,0])


% disp(['Memory measure: ' w.measure '-----------------------' ])
% disp(' ')
% disp('SIMILARITY x VALENCE ANOVA')      
% disp(['Similarity:     p=' num2str(anova.p(1))] )
% disp(['Valence:       p=' num2str(anova.p(2))] )
% disp( ['Sim x Val:     p=' num2str(anova.p(3))])
% disp(' ')
% disp('REWARD EFFECTS')
% disp( ['Overall:          p=' num2str(graphs.diff_ttests(1,2))] )
% disp(['Similar:          p=' num2str(graphs.diff_ttests(2,2))] )
% disp(['Dissimilar:     p=' num2str(graphs.diff_ttests(3,2))] )
% disp(' ')
% disp('SIMILARITY EFFECTS')
% disp( ['Overall:          p=' num2str(graphs.diff_ttests(6,2))] )
% disp(['Reward:         p=' num2str(graphs.diff_ttests(6,2))] )
% disp(['Neutral:          p=' num2str(graphs.diff_ttests(7,2))] )
end

end


%% Correlations
if w.correlations==1

    % Calculate subject-specific index of interaction (specifically: reward
    % effect for Sim only, but not Dis)
    data.summ.Interaction.Sim4Rew.Rew4Sim=(data.Cell.Sim_R-data.Cell.Sim_N)+(data.Cell.Sim_R-data.Cell.Dis_R);
    data.summ.Interaction.Rew.SimvsDis=(data.Cell.Sim_R-data.Cell.Sim_N)-(data.Cell.Dis_N-data.Cell.Dis_R);
%     data.summ.Overall_Valfx=data.Summ.Rew-data.Summ.Neu;
%     data.summ.Overall_Simfx=data.Summ.Sim-data.Summ.Dis;
    
    % Data prefixes
    w.corr=corr;
    for i=1:size(corr.iv_pers,1)
        ws.count=1;
        eval(['ws.alldat=alldata.sub_info.' corr.iv_pers{i} ';']); % Implement exclusion
        ws.dat=ws.alldat(find(w.include(:)));
        eval(['data.sub_info.' corr.iv_pers{i} '=ws.dat;'])
        corr.iv_pers{i}=horzcat('data.sub_info.',corr.iv_pers{i});
    end
    for i=1:size(corr.dv,1)
        if strcmp(corr.dv{i}(1),'O')==1 || strcmp(corr.dv{i}(1),'I')==1
            corr.dv{i}=horzcat('data.summ.',corr.dv{i});
        else
            corr.dv{i}=horzcat('data.cell.',corr.dv{i});
        end
    end
    for i=1:size(corr.iv_mem,1)
        eval(['ws.alldat=alldata.d_' corr.iv_mem{i} '.Overall.Overall;']);
        ws.count=1;
        for j=1:size(ws.alldat,1)
            if w.include(j)==1
                ws.dat(ws.count)=ws.alldat(j);
                ws.count=ws.count+1;
            end
        end
        eval(['d_meanmemory.' corr.iv_mem{i} '=ws.dat;']);
        ws=[];
        corr.iv_mem{i}=horzcat('d_meanmemory.',corr.iv_mem{i});
    end
    d_meanmemory.subjects=alldata.sub_info.Subject;

    % Pre-set correlations
    w.corrcount=1;
    corr.iv=vertcat(corr.iv_pers,corr.iv_mem);
    w.corr.iv=vertcat(w.corr.iv_pers,w.corr.iv_mem);
    for j=1:size(corr.dv,1) %  corr.result; Col 1: IV, Col 2: dv, Col 3=rho,Col 4=p, Col5=sig?
        wcc.dvname=corr.dv{j};
        for i=1:size(corr.iv,1)
            wc.ivname=corr.iv{i};
            corr.result{w.corrcount,1}=w.corr.iv{i};
            corr.result{w.corrcount,2}=w.corr.dv{j};
            eval(['[wc.r wc.p]=corrcoef(' wcc.dvname ',' wc.ivname ');'])
            corr.result{w.corrcount,3}=wc.r(1,2);
            corr.result{w.corrcount,4}=wc.p(1,2);
            if corr.result{w.corrcount,4}<0.051
                corr.result{w.corrcount,5}=1;
            elseif corr.result{w.corrcount,4}<0.1
                corr.result{w.corrcount,5}=0.5;
            else
                corr.result{w.corrcount,5}=0;
            end
            if corr.result{w.corrcount,5}>0.4 && corr.result{w.corrcount,3}<0
                corr.result{w.corrcount,5}=corr.result{w.corrcount,5}*-1;
            end
            w.corrcount=w.corrcount+1; 
            clear wc
        end
        corr.result{w.corrcount,1}='NEXT';
        corr.result{w.corrcount,2}=0;
        corr.result{w.corrcount,3}=1;
        corr.result{w.corrcount,4}=1;
        corr.result{w.corrcount,5}=1;
        w.corrcount=w.corrcount+1; 
        wcc=[];
    end

    % Memory effects correlated?
    for j=1:size(corr.effects,1) %  corr.result; Col 1: IV, Col 2: dv, Col 3=rho,Col 4=p, Col5=sig?
        wcc.var1=corr.effects{j};
        if  strcmp(wcc.var1(1),'O')==1 || strcmp(wcc.var1(1),'I')==1
            wcc.var1n=['data.summ.' wcc.var1];
        else
            wcc.var1n=['data.cell.' wcc.var1];
        end
        for i=1:size(corr.effects,1)
            if i>j
                    wc.var2=corr.effects{i};
                    corr.result{w.corrcount,1}=wcc.var1;
                    corr.result{w.corrcount,2}=wc.var2;
                    if  strcmp(wc.var2(1),'O')==1 || strcmp(wc.var2(1),'I')==1
                         wc.var2n=['data.summ.' wc.var2];
                    else
                        wc.var2n=['data.cell.' wc.var2];
                    end
                    eval(['[wc.r wc.p]=corrcoef(' wcc.var1n ', ' wc.var2n ');'])
                    corr.result{w.corrcount,3}=wc.r(1,2);
                    corr.result{w.corrcount,4}=wc.p(1,2);
                    if corr.result{w.corrcount,4}<0.051 % Mark sig
                        corr.result{w.corrcount,5}=1;
                    elseif corr.result{w.corrcount,4}<0.1
                        corr.result{w.corrcount,5}=0.5;
                    else
                        corr.result{w.corrcount,5}=0;
                    end
                    if corr.result{w.corrcount,5}>0.4 && corr.result{w.corrcount,3}<0
                        corr.result{w.corrcount,5}=corr.result{w.corrcount,5}*-1;
                    end
                    w.corrcount=w.corrcount+1; 
                clear wc
            end
        end
        corr.result{w.corrcount,1}='NEXT';
        corr.result{w.corrcount,2}=0;
        corr.result{w.corrcount,3}=1;
        corr.result{w.corrcount,4}=1;
        corr.result{w.corrcount,5}=1;
        w.corrcount=w.corrcount+1; 
        wcc=[];
    end
    
    % [OUTPUT] Array for printing to txt (text statements - format here to print as table)
    w.c=2; corr.sigres{1,1}=['Correlations:   '  w.measure];
    for i=1:size(corr.result,1)
        if strcmp('NEXT', corr.result{i,1})==1 % Buffer
            corr.sigres{w.c,1}='-------------------------------------------------------------------------------------------------';
            w.c=w.c+1;
        elseif abs(corr.result{i,5})==1 % Significant effects
            if corr.result{i,4}<0.001
                corr.sigres{w.c,1}=['[  ' corr.result{i,1} '  -  ' corr.result{i,2} '  ]   r= ' num2str(corr.result{i,3})  '  ,  p= ' num2str(corr.result{i,4})  '   ***'];
            elseif corr.result{i,4}<0.01
                corr.sigres{w.c,1}=['[  ' corr.result{i,1} '  -  ' corr.result{i,2} '  ]   r= ' num2str(corr.result{i,3})  '  ,  p= ' num2str(corr.result{i,4})   '   **'];
            else
                corr.sigres{w.c,1}=['[  ' corr.result{i,1} '  -  ' corr.result{i,2} '  ]   r= ' num2str(corr.result{i,3})  '  ,  p= ' num2str(corr.result{i,4})   '   *'];
            end
            w.c=w.c+1;
%         elseif abs(corr.result{i,5})==0.5 % Trends
%             corr.sigres{w.c,1}=['[  ' corr.result{i,1} '  -  ' corr.result{i,2} '  ]   r= ' num2str(corr.result{i,3})  '  ,  p= ' num2str(corr.result{i,4})  ];
%             w.c=w.c+1;
        end
    end
    w.printcorrok=print2txt([where.where filesep '3 Analysis inputs'], ['(' date ') Correlations - ' num2str(w.whichmeasure) ' ' w.measure], corr.sigres);
    disp('Printing correlations')
    disp(w.printcorrok)

end

%% Save variables (for use later with neural results)

% Only save if all subjects are included
datatypes=fieldnames(alldata); dd='''';
for i=1:size(datatypes,1)
    eval([datatypes{i} '=alldata.' datatypes{i}  ';']);
    dd=[dd datatypes{i} ];
    if i<size(datatypes,1)
        dd=[dd ''' , '''];
    else
        dd=[dd ''''];
    end
end
filename=[ '''' where.where filesep '(' date ') Memory scores.mat''' ];
eval(['save(' filename ', ' dd ');'])
