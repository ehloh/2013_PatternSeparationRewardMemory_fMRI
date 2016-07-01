% Is memory affected by outcome presentation? 
clear all; close all hidden, clc

% Subjects
log.subs_manu_v4_6_fmri={'p01_CW';'p03_EA';'p04_JL';'p07_LH';'p08_AM';'p09_CN';'p10_AB';'p11_SS';'p12_IL';'p14_SJ';'p16_TH';'p17_RA';'p18_KB';'p19_CN';'p20_JB';'p21_SH';'p22_EK';'p23_IS';'p24_LL';'p25_BS';'p26_MC';'p27_EW';'p28_CN';'p29_ET'}; 
log.subs_manu_v4_6_beh= {'p01_CW';'p03_EA';'p04_JL';'p06_YL';'p07_LH';'p08_AM';'p09_CN';'p10_AB';'p11_SS';'p12_IL';'p14_SJ';'p16_TH';'p17_RA';'p18_KB';'p19_CN';'p20_JB';'p21_SH';'p22_EK';'p23_IS';'p24_LL';'p25_BS';'p26_MC';'p27_EW';'p28_CN';'p29_ET'};
log.subs_manu_v4_6_beh_knowOK= {'p01_CW';'p03_EA';'p04_JL';'p06_YL';'p07_LH';'p08_AM';'p09_CN';'p10_AB';'p11_SS';'p12_IL';'p14_SJ';'p16_TH';'p17_RA';'p18_KB';'p19_CN';'p20_JB';'p22_EK';'p23_IS';'p24_LL';'p25_BS';'p26_MC';'p28_CN';'p29_ET'}; 
log.subs_manu_v4_6_beh_sureknowOK= {'p01_CW';'p03_EA';'p04_JL';'p06_YL';'p07_LH';'p08_AM';'p10_AB';'p11_SS';'p12_IL';'p14_SJ'; 'p17_RA';'p18_KB';'p19_CN';'p20_JB';'p21_SH'; 'p23_IS';'p24_LL';'p25_BS';'p26_MC';'p27_EW';'p28_CN';'p29_ET'}; 



% [ WHICH subjects? ] ################
% log.subjects= log.subs_manu_v4_6_fmri;   
log.subjects= log.subs_manu_v4_6_beh;   
% log.subjects= log.subs_manu_v4_6_beh_knowOK;  disp('KNOW mem score subject selection !!!!!') 
% log.subjects= log.subs_manu_v4_6_beh_sureknowOK;   disp('SUREKNOW mem score subject selection !!!!!') 

for o=1:1
log.n_subjs = length(log.subjects); 


    where.beh= 'D:\Dropbox\SCRIPPS\2a ContextMem behaviour\2 Behavioural Data';
    log.score_prefix=[]; 
    
    
    
end 



%% Load data

log.cells={'sr' 'sn' 'dr' 'dn'};  log.memtypes={'dpr';'r';'k';'sr';'sk';'h';'sh'; 'gr';'gk'; 'gh';}; % Order of al mem things 
d_mem= [log.memtypes repmat({nan(log.n_subjs,4)},length(log.memtypes),1)];  d_overallmem= [ [{' '}; log.subjects] [log.memtypes'; cell(log.n_subjs,length(log.memtypes))]];
d_encoutcome =nan(log.n_subjs,1); 
% repmat({nan(log.n_subjs,4)},length(log.memtypes),1)];
subjdata=[log.subjects cell(log.n_subjs,1)];  % 2: modified trialsts ('col'), 3: struc w fa etc
for o1=1:1 % Column settings     
    
    % Data colums 
    col.itemstim =1;
    col.simdis=2 ;
    col.rewneu=3 ;
    col.oldnew=4 ;
    col.sureguess =5;
    col.remknow=6;
    col.roc=7;
    col.enccorrect=8;
    col.outcomepres=9;

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
    col.m.RemKnow=14;     
end
for s=1:log.n_subjs
    ws.e1= load([where.beh filesep log.subjects{s} filesep log.subjects{s}  '_file_2encodingFMRI_b1.mat']); ws.e1= ws.e1.encoding.trialstats; 
    ws.e2= load([where.beh filesep log.subjects{s} filesep log.subjects{s}  '_file_2encodingFMRI_b2.mat']); ws.e2= ws.e2.encoding.trialstats; 
    ws.m = load([where.beh filesep log.subjects{s} filesep log.subjects{s}  '_file_3memorytest.mat']); ws.m = ws.m.memtest.trialstats;
    ws.e2(:, col.et.Trialnum)=ws.e2(:, col.et.Trialnum)+ size(ws.e1,1);  ws.e=[ws.e1; ws.e2]; 
    ws.ewhere= [[ws.e(:, col.et.StimItem(1)) ws.e(:, col.et.Trialnum) ws.e(:, col.et.AccItem(1))]; [ws.e(:, col.et.StimItem(2)) ws.e(:, col.et.Trialnum) ws.e(:, col.et.AccItem(2))]; [ws.e(:, col.et.StimItem(3)) ws.e(:, col.et.Trialnum) ws.e(:, col.et.AccItem(3))]];
    ws.m( ws.m(:, col.m.RemKnow)==1 & ws.m(:, col.m.SureGuess)==1, col.m.Roc) =6;  % Mark ROC 
    ws.m( ws.m(:, col.m.RemKnow)==1 & ws.m(:, col.m.SureGuess)==2, col.m.Roc) =5;
    ws.m( ws.m(:, col.m.RemKnow)==2 & ws.m(:, col.m.SureGuess)==1, col.m.Roc) =4;
    ws.m( ws.m(:, col.m.RemKnow)==2 & ws.m(:, col.m.SureGuess)==2, col.m.Roc) =3; 
    ws.m( ws.m(:, col.m.OldNew)==2 & ws.m(:, col.m.SureGuess)==2, col.m.Roc) =2;
    ws.m( ws.m(:, col.m.OldNew)==2 & ws.m(:, col.m.SureGuess)==1, col.m.Roc) =1;
    ws.mo=ws.m ( ws.m (:, col.m.OldFoil)==1, :);  ws.d(:, [col.roc col.itemstim col.simdis col.rewneu col.oldnew col.sureguess col.remknow])=  ws.mo(:, [col.m.Roc col.m.ItemStim col.m.SimDis col.m.RewNeu col.m.OldNew col.m.SureGuess  col.m.RemKnow]) ;
    ws.mf=ws.m ( ws.m (:, col.m.OldFoil)==2, :); 
    
    % Mark item-specific things from enc stage 
    for t=1:size(ws.d,1)
        ws.d(t, col.outcomepres) =  ws.e( ws.ewhere(ws.ewhere(:,1)==ws.d(t,col.itemstim),2), col.et.OutcomePres) ; 
        ws.d(t, col.enccorrect) = ws.ewhere(ws.ewhere(:,1)==ws.d(t,col.itemstim),3); 
    end  
    
    
    d_encoutcome(s,1) = sum(sum(ws.e(:, [col.et.OutcomeItem1 col.et.OutcomeItem2 col.et.OutcomeItem3]) )); 
    d_encoutcome(s,2) = d_encoutcome(s,1)/ (0.5*size(ws.e,1)*3); 
    
    % Enc correct only
    ws.d= ws.d(ws.d(:, col.enccorrect) ==1,:); 
%     ws.d= ws.d(ws.d(:, col.outcomepres) ==1,:);  if s==1; input('Outcome presented trials ONLY. Continue?'); log.score_prefix='oy_';end 
%     ws.d= ws.d(ws.d(:, col.outcomepres) ==0,:);  if s==1; input('Outcome NOT presented trials only. Continue?');  log.score_prefix='on_'; end 
    
    % FA rates 
    ws.fa.d= (sum(ws.mf(:, col.m.OldNew)==1)+0.5)/size(ws.mf,1); 
    ws.fa.rem= mean(ws.mf(:, col.m.RemKnow)==1); 
    ws.fa.know= mean(ws.mf(:, col.m.RemKnow)==2); 
    ws.fa.surerem= mean(ws.mf(:, col.m.Roc)==6); 
    ws.fa.sureknow= mean(ws.mf(:, col.m.Roc)==4); 
    ws.fa.hit = mean(ws.mf(:, col.m.OldNew)==1); 
    ws.fa.surehit = mean( ws.mf(:, col.m.OldNew)==1& ws.mf(:, col.m.SureGuess)==1 ); 
    ws.fa.guessrem= mean(ws.mf(:, col.m.Roc)==5); 
    ws.fa.guessknow= mean(ws.mf(:, col.m.Roc)==3); 
    ws.fa.guesshit = mean( ws.mf(:, col.m.OldNew)==1& ws.mf(:, col.m.SureGuess)==2 ); 
    subjdata{s,3}= [ws.fa.d  ws.fa.rem  ws.fa.know  ws.fa.surerem  ws.fa.sureknow  ws.fa.hit ws.fa.surehit]; 
%     error
    
    for m=1:length(subjdata{s,3})  % NOT for guess 
        d_mem{m,2}(s,1) = f_calcmemscore(m, ws.d(ws.d(:, col.simdis)==1 & ws.d(:, col.rewneu)==1, col.roc), subjdata{s,3}(m)); 
        d_mem{m,2}(s,2) = f_calcmemscore(m, ws.d(ws.d(:, col.simdis)==1 & ws.d(:, col.rewneu)==0, col.roc), subjdata{s,3}(m)); 
        d_mem{m,2}(s,3) = f_calcmemscore(m, ws.d(ws.d(:, col.simdis)==2 & ws.d(:, col.rewneu)==1, col.roc), subjdata{s,3}(m)); 
        d_mem{m,2}(s,4) = f_calcmemscore(m, ws.d(ws.d(:, col.simdis)==2 & ws.d(:, col.rewneu)==0, col.roc), subjdata{s,3}(m)); 
        d_mem{m,2}(s,5) = f_calcmemscore(m, ws.d(:,col.roc), subjdata{s,3}(m));  d_overallmem{s+1,m+1} = d_mem{m,2}(s,5); 
    end   
    
    
    
    
    % Unsure memory only 
    m=find(strcmp(log.memtypes, 'gh'));
    d_mem{m,2}(s,1) = sum(ws.d(:, col.simdis)==1 & ws.d(:, col.rewneu)==1    &   ws.d(:, col.oldnew)==1 & ws.d(:, col.sureguess)==2)/sum(ws.d(:, col.simdis)==1 & ws.d(:, col.rewneu)==1) - ws.fa.guesshit; 
    d_mem{m,2}(s,2) = sum(ws.d(:, col.simdis)==1 & ws.d(:, col.rewneu)==0    &   ws.d(:, col.oldnew)==1 & ws.d(:, col.sureguess)==2)/sum(ws.d(:, col.simdis)==1 & ws.d(:, col.rewneu)==0) - ws.fa.guesshit; 
    d_mem{m,2}(s,3) = sum(ws.d(:, col.simdis)==2 & ws.d(:, col.rewneu)==1    &   ws.d(:, col.oldnew)==1 & ws.d(:, col.sureguess)==2)/sum(ws.d(:, col.simdis)==2 & ws.d(:, col.rewneu)==1) - ws.fa.guesshit; 
    d_mem{m,2}(s,4) = sum(ws.d(:, col.simdis)==2 & ws.d(:, col.rewneu)==0    &   ws.d(:, col.oldnew)==1 & ws.d(:, col.sureguess)==2)/sum(ws.d(:, col.simdis)==2 & ws.d(:, col.rewneu)==0) - ws.fa.guesshit; 
    d_mem{m,2}(s,5) = sum(ws.d(:, col.oldnew)==1 & ws.d(:, col.sureguess)==2)- ws.fa.guesshit;    d_overallmem{s+1,m+1} = d_mem{m,2}(s,5); 
    ws.dg=ws.d(ws.d(:, col.sureguess)==2 & ws.d(:, col.oldnew)==1 ,:);    
    m=find(strcmp(log.memtypes, 'gr'));
    d_mem{m,2}(s,1) = f_calcmemscore(2, ws.dg(ws.dg(:, col.simdis)==1 & ws.dg(:, col.rewneu)==1, col.roc), ws.fa.guessrem);
    d_mem{m,2}(s,2) = f_calcmemscore(2, ws.dg(ws.dg(:, col.simdis)==1 & ws.dg(:, col.rewneu)==0, col.roc), ws.fa.guessrem);
    d_mem{m,2}(s,3) = f_calcmemscore(2, ws.dg(ws.dg(:, col.simdis)==2 & ws.dg(:, col.rewneu)==1, col.roc), ws.fa.guessrem);
    d_mem{m,2}(s,4) = f_calcmemscore(2, ws.dg(ws.dg(:, col.simdis)==2 & ws.dg(:, col.rewneu)==0, col.roc), ws.fa.guessrem);
    d_mem{m,2}(s,5) = f_calcmemscore(2, ws.dg(:, col.roc), ws.fa.guessrem);   d_overallmem{s+1,m+1} = d_mem{m,2}(s,5); 
    m=find(strcmp(log.memtypes, 'gk'));
    d_mem{m,2}(s,1) = f_calcmemscore(3, ws.dg(ws.dg(:, col.simdis)==1 & ws.dg(:, col.rewneu)==1, col.roc), ws.fa.guessknow);
    d_mem{m,2}(s,2) = f_calcmemscore(3, ws.dg(ws.dg(:, col.simdis)==1 & ws.dg(:, col.rewneu)==0, col.roc), ws.fa.guessknow);
    d_mem{m,2}(s,3) = f_calcmemscore(3, ws.dg(ws.dg(:, col.simdis)==2 & ws.dg(:, col.rewneu)==1, col.roc), ws.fa.guessknow);
    d_mem{m,2}(s,4) = f_calcmemscore(3, ws.dg(ws.dg(:, col.simdis)==2 & ws.dg(:, col.rewneu)==0, col.roc), ws.fa.guessknow);
    d_mem{m,2}(s,5) = f_calcmemscore(3, ws.dg(:, col.roc), ws.fa.guessknow);   d_overallmem{s+1,m+1} = d_mem{m,2}(s,5); 
    
    subjdata{s,2}=ws.d; 
    ws=[]; 
end 
% openvar t_mem, openvar d_overallmem, [r c]= find(cell2mat( d_overallmem(2:end,2:end)) <0); disp([r c]+1)

mean(d_encoutcome)
std(d_encoutcome)

% Write to table 
t_mem=[{'Subjects'};  log.subjects];  k=2; 
for m=1:length(log.memtypes)
    for c=1:length(log.cells)
        t_mem=[t_mem   [{[log.score_prefix log.memtypes{m} '_' log.cells{c}]};  num2cell(d_mem{m,2}(:,c))] ]; 
    end
end

%% Sim x Val Plot and ANOVA 

% Implement requesting of memtypes 
d_mem=[  d_mem(strcmp(d_mem(:,1),'dpr'),:); d_mem(strcmp(d_mem(:,1),'r'),:); d_mem(strcmp(d_mem(:,1),'k'),:); d_mem(strcmp(d_mem(:,1),'sh'),:); d_mem(strcmp(d_mem(:,1),'sr'),:); d_mem(strcmp(d_mem(:,1),'sk'),:); d_mem(strcmp(d_mem(:,1),'gh'),:); d_mem(strcmp(d_mem(:,1),'gr'),:); d_mem(strcmp(d_mem(:,1),'gk'),:)]; 
log.memtypes= d_mem(:,1); 

close all hidden
f.plotcols=3;  f.figwidth= 1800; f.figheight=400; f.fontsize=20; f.fontsize_title=25;f.fontname='PT Sans Caption';
f.subplot_VerHorz=[0.1 0.15]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.1 0.05];
figure('Name', 'Mem scores', 'NumberTitle', 'off', 'Position', [100 550 f.figwidth f.figheight], 'Color', 'w');  k=1;
for m=1:length(log.memtypes)
    subtightplot(ceil((length(log.memtypes)+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    barwitherr( reshape(std(d_mem{m,2}(:,1:4))./sqrt(log.n_subjs),2,2), reshape(mean(d_mem{m,2}(:,1:4)), 2,2))
    set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Similar' 'Dissimilar' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
    title(log.memtypes{m},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
    %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
    %                 ylabel('Parameter estimates (mean)', 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
%     if m==1,  legend('Similar-Reward', 'Similar-Neutral', 'Dissimilar-Reward', 'Dissimilar-Neutral'), end
    xlim([0.5 2.5])
    ylim( [mean( reshape(d_mem{m,2}(:,1:4), size(d_mem{m,2},1)*4,1) )-max(std(reshape(d_mem{m,2}(:,1:4), size(d_mem{m,2},1)*4,1)))*1.2            mean(reshape(d_mem{m,2}(:,1:4), size(d_mem{m,2},1)*4,1))+max(std(reshape(d_mem{m,2}(:,1:4), size(d_mem{m,2},1)*4,1)))*1.2] )
    
end 
for m=1:length(log.memtypes)  % STATS 
    
    disp([' #######  '  log.memtypes{m} ' ############################'])
    d_mem{m,3}=teg_repeated_measures_ANOVA(d_mem{m,2}(:,1:4) ,[2 2], {'Sim' 'Val '}); d_mem{m,3} =d_mem{m,3}.R; 
    c=1; disp( ['ME Sim: F(' num2str(d_mem{m,3}(c,2)) ','  num2str(d_mem{m,3}(c,3)) ')=' num2str(d_mem{m,3}(c,1),3) ',   p='  num2str(d_mem{m,3}(c,4),3) ])
    c=2; disp( ['ME Val: F(' num2str(d_mem{m,3}(c,2)) ','  num2str(d_mem{m,3}(c,3)) ')=' num2str(d_mem{m,3}(c,1),3) ',   p='  num2str(d_mem{m,3}(c,4),3) ])
    c=3; disp( ['Sim x Val : F(' num2str(d_mem{m,3}(c,2)) ','  num2str(d_mem{m,3}(c,3)) ')=' num2str(d_mem{m,3}(c,1),3) ',   p='  num2str(d_mem{m,3}(c,4),3) ])
    disp(' ') 
end


% Simple effects 
m=4; 
[h p ci st]=ttest( d_mem{m,2}(:,1) - d_mem{m,2}(:,4) ); 
disp([log.memtypes{m} '   -  t(' num2str(st.df ) ') = '  num2str(st.tstat) '  , p=' num2str(p)])



%% Is memory affected by outcome presentation?

    



