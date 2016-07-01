% Ad hoc behavioural analysis script. Analysis should be titled
clear all; close all hidden; clc




for o1=1:1 % General setup
    
    log.behoksubjects_n25={'p01_CW';'p03_EA';'p04_JL';'p06_YL';'p07_LH';'p08_AM';'p09_CN';'p10_AB';'p11_SS';'p12_IL';'p14_SJ';'p16_TH';'p17_RA';'p18_KB';'p19_CN';'p20_JB';'p21_SH';'p22_EK';'p23_IS';'p24_LL';'p25_BS';'p26_MC';'p27_EW';'p28_CN';'p29_ET';};  % Exclude d'<0.3
    log.specificsubjects=log.behoksubjects_n25;
    
    % Paths and subjects
    addpath('D:\Dropbox\SANDISK\1 fMRI analysis')
    where.data='D:\Dropbox\SANDISK\2a fMRI behaviour analysis\2 Behavioural Data';
    cd(where.data)
    log.allsubs=dir('p*'); log.allsubs=cellstr(char(log.allsubs(:).name));
    log.allsublog=[[{'Subjects'}; log.allsubs] [{'OK'}; num2cell(ones(length(log.allsubs),1))]];
    [log.subjects log.n_subjs newdatalog] = f_selectsubjects(log.allsublog, log.specificsubjects, log.allsublog, 'OK');
    
end


%% Load data 
% subjdata: Col 2=Encoding data, Col 3= memtest

for o1=1:1 % Encoding stage documentation
% Details of parameter and data files ('par & 'data')
%
% [Scene stim]   Col 1-2: Scene stim - Index (1-3)    
%                       Col 3:   Scene stim - Type. dissimilar (1=Similar, 2=Dissimilar)
%                       Col 4:   Scene stim - Reward available (1=Rewarded, 0=Neutral)
% 
%                       Col 5:      Trial number during encoding
%                       Col 6:     	
%                       Col 7:      Total winnings actually won on this trial (out of 3)
%
% [Item stim]     Col 8:   	Outcome Item 1 (1=win, 0=nothing, -1=mild loss)
%                       Col 9:   	Outcome Item 2
%                       Col 10:    	Outcome item 3
%
%                       Col 11-13:  Item stim - Stimulus index(Ordered in cols, by presentation)
%                       Col 14-16:  Item stim - Reward (effectively) available, for this item (1=Yes, 0=No)
%                       Col 17-19:  Item stim - Item's true semantic category, Items1-3 (1=Natural,2=Manmade)  
%                       Col 20:     Item stim - #1: Keypress
%                       Col 21:     Item stim - #1: RT
%                       Col 22:     Item stim - #2: Keypress
%                       Col 23:     Item stim - #2: RT
%                       Col 24:     Item stim - #3: Keypress
%                       Col 25:     Item stim - #3: RT
%
%                       Col 26:     Item stim - Accuracy (Item #1)
%                       Col 27:     Item stim - Accuracy (Item #2)
%                       Col 28:     Item stim - Accuracy (Item #3)
%                       Col 29:     Feedback presented? 
% [for fMRI]
%
% TRIAL EVENTS:
%       - Fixation
%       - Context only
%       - Item 1
%       - Item 2
%       - Item 3
%       
%                       Col 30:     START trial
%                       Col 31:     Fixation onset
%                       Col 32:     Context onset
%                       Col 33:     Item 1 onset
%                       Col 34:     Item 1 Keypress
%                       Col 35:     Item 2 onset
%                       Col 36:     Item 2 Keypress
%                       Col 37:     Item 3 onset
%                       Col 38:     Item 3 Keypress
%                       Col 39:     Outcome onset
%                       Col 40:     END trial
%
end

for o1=1:1 % Memtest documentation
% Details of parameters file ('par')
%       Col 1: Trial number (memory test) - if not tested, #=0
%       Col 2: Item stimulus index
%       Col 3: Old or New/foil (actual) (1=Old, 2=New)
% [Encoding-stage]
%       Col 4: [Enc] How many potentially-rewarded items in the enc trial (1-3)
%       Col 5: [Enc] How much actually won in the enc trial (1-3)
%       Col 6: Trial number (during encoding)
%       Col 7: Outcome type - for Item+Scene, conjunctive (1=Reward, 0=None)
%                  - this indicates effective outcome availability
%       Col 8-9: Scene stimuli index(1a-2b)
%       Col 10: Scene stimuli type (1=Similar, 2=Dissimilar)
%       Col 11: Outcome type - for Scene (1=Rewarded, 0=Neutral)
% [Memory test stage]
%       Col 12: Old/New Judgment - Keypress (1=Old, 2=New)
%       Col 13: Old/New Judgment - RT
%       Col 14: Remember/Know Judgment - Keypress (1=Rem, 2=Know; 0=n/a i.e.NEW)
%       Col 15: Remember/Know Judgment - RT
%       Col 16: Old/New Confidence - Keypress (1=Sure/Strong, 2=Guess/Weak)
%       Col 17: Old/New Confidence - RT
%       Col 18: [BLANK] Recognize Context - Keypress (Index of chosen scene)
%       Col 19: [BLANK] Recognize Context - RT
%       Col 20: [BLANK] Recognize Position - Keypress (Position1-4)
%       Col 21: [BLANK] Recognize Position - RT
%       ----------
%       Col 26: Accuracy - Item recognition
%       Col 27: 
%       Col 28: [Enc] Item semantic category (1=Not manmade, 2=Manmade)
%       Col 29: [Enc] Item stimulus: correctly judged? (1=Yes, 0=No)

       
end


for s=1:log.n_subjs
    
    % Load all data 
    ws.b1=load([log.subjects{s} filesep log.subjects{s} '_file_2encodingFMRI_b1.mat']);
    ws.b2=load([log.subjects{s} filesep log.subjects{s} '_file_2encodingFMRI_b2.mat']);
    ws.b1.encoding.trialstats(:, 6)=1; ws.b2.encoding.trialstats(:, 6)=2;
    ws.enc=[ws.b1.encoding.trialstats;  ws.b2.encoding.trialstats];  ws.enc(:,5)=1:size(ws.enc,1);
    ws.mem=load([log.subjects{s} filesep log.subjects{s} '_file_3memorytest.mat']);
    ws.mem=ws.mem.memtest.trialstats;
    

    % Save data 
    subjdata{s,1}=log.subjects{s};
    subjdata{s,2}=ws.enc;
    subjdata{s,3}= ws.mem;
end



%% (1)  Context-related memory effects
% Is it the case that memory for one object (in a certain context trial) is more closely related to memory for other 
% objects on the same context trial?
% subjdata Col 4: Results for this analysis

d_coherescore=zeros(log.n_subjs,4); % SR, SN, DR, SN
for s=1: log.n_subjs
    ws.m=subjdata{s,3};
    ws.e=subjdata{s,2};
    
    for t=1:size(ws.e,1)
        wt.itemmem=[nan nan nan];
        for in=1:3
            wt.itemmem(in)= ws.m(find(ws.m(:,2) == ws.e(t, 10+in)), 12)==1;
        end
        
        % Col 30: Coherent or not
        disp('Manually check which coherence measure you are checking !! ')
%         if sum(wt.itemmem) > 1.5 
        if sum(wt.itemmem)==3  ||  sum(wt.itemmem)==0 % Cohere= all hit, or all miss
            ws.m(t,  30)=1;
        else ws.m(t,  30)=0;
        end
    end
    
    % Log coherence scores per context type
    subjdata{s,4}.CohereScore_SR=sum(ws.m(ws.m(:,10)==1 &  ws.m(:,11)==1, 30));
    d_coherescore(s,1)=subjdata{s,4}.CohereScore_SR;
    subjdata{s,4}.CohereScore_SN=sum(ws.m(ws.m(:,10)==1 &  ws.m(:,11)==0, 30));
    d_coherescore(s,2)=subjdata{s,4}.CohereScore_SN;
    subjdata{s,4}.CohereScore_DR=sum(ws.m(ws.m(:,10)==2 &  ws.m(:,11)==1, 30));
    d_coherescore(s,3)=subjdata{s,4}.CohereScore_DR;
    subjdata{s,4}.CohereScore_DN=sum(ws.m(ws.m(:,10)==2 &  ws.m(:,11)==0, 30));
    d_coherescore(s,4)=subjdata{s,4}.CohereScore_DN;
end

openvar d_coherescore  % Copy numbers into SPSS



% [Scene stim]   Col 1-2: Scene stim - Index (1-3)    
%                       Col 3:   Scene stim - Type. dissimilar (1=Similar, 2=Dissimilar)
%                       Col 4:   Scene stim - Reward available (1=Rewarded, 0=Neutral)
% 
%                       Col 5:      Trial number during encoding
%                       Col 6:     	
%                       Col 7:      Total winnings actually won on this trial (out of 3)
%
% [Item stim]     Col 8:   	Outcome Item 1 (1=win, 0=nothing, -1=mild loss)
%                       Col 9:   	Outcome Item 2
%                       Col 10:    	Outcome item 3
%
%                       Col 11-13:  Item stim - Stimulus index(Ordered in cols, by presentation)
%                       Col 14-16:  Item stim - Reward (effectively) available, for this item (1=Yes, 0=No)
%                       Col 17-19:  Item stim - Item's true semantic category, Items1-3 (1=Natural,2=Manmade)  
%                       Col 20:     Item stim - #1: Keypress
%                       Col 21:     Item stim - #1: RT
%                       Col 22:     Item stim - #2: Keypress
%                       Col 23:     Item stim - #2: RT
%                       Col 24:     Item stim - #3: Keypress
%                       Col 25:     Item stim - #3: RT
%
%                       Col 26:     Item stim - Accuracy (Item #1)
%                       Col 27:     Item stim - Accuracy (Item #2)
%                       Col 28:     Item stim - Accuracy (Item #3)
%                       Col 29:     Feedback presented? 



for o1=1:1 % Documentation
% Details of parameters file ('par')
%       Col 1: Trial number (memory test) - if not tested, #=0
%       Col 2: Item stimulus index
%       Col 3: Old or New/foil (actual) (1=Old, 2=New)
% [Encoding-stage]
%       Col 4: [Enc] How many potentially-rewarded items in the enc trial (1-3)
%       Col 5: [Enc] How much actually won in the enc trial (1-3)
%       Col 6: Trial number (during encoding)
%       Col 7: Outcome type - for Item+Scene, conjunctive (1=Reward, 0=None)
%                  - this indicates effective outcome availability
%       Col 8-9: Scene stimuli index(1a-2b)
%       Col 10: Scene stimuli type (1=Similar, 2=Dissimilar)
%       Col 11: Outcome type - for Scene (1=Rewarded, 0=Neutral)
% [Memory test stage]
%       Col 12: Old/New Judgment - Keypress (1=Old, 2=New)
%       Col 13: Old/New Judgment - RT
%       Col 14: Remember/Know Judgment - Keypress (1=Rem, 2=Know; 0=n/a i.e.NEW)
%       Col 15: Remember/Know Judgment - RT
%       Col 16: Old/New Confidence - Keypress (1=Sure/Strong, 2=Guess/Weak)
%       Col 17: Old/New Confidence - RT
%       Col 18: [BLANK] Recognize Context - Keypress (Index of chosen scene)
%       Col 19: [BLANK] Recognize Context - RT
%       Col 20: [BLANK] Recognize Position - Keypress (Position1-4)
%       Col 21: [BLANK] Recognize Position - RT
%       ----------
%       Col 26: Accuracy - Item recognition
%       Col 27: 
%       Col 28: [Enc] Item semantic category (1=Not manmade, 2=Manmade)
%       Col 29: [Enc] Item stimulus: correctly judged? (1=Yes, 0=No)

       
end
