% Get onsets for set up of 1st level
clear all; close all hidden; clc;
% where.where='/Users/EleanorL/Dropbox/SCRIPPS/1 ContextMem fMRI'; where.data_brain='/Users/EleanorL/Desktop/1 CONTEXT fmri data/1 MRI data'; where.data_beh=[where.where filesep '2 Behavioural data'];
where.where='D:\Dropbox\SCRIPPS\1 ContextMem fMRI';  where.data_brain='D:\1 [Context-Memory] fMRI Data\1 MRI data';  where.data_beh=[where.where filesep '2 Behavioural data'];


% Requested analysis
log.specificsubjects={};
% log.specificsubjects={'p01_CW'};
request.simplifydata=1; % 
request.specifydata='(02-Apr-2013) Simplified data.mat'; 
request.saveonsets=1;
request.printcontextmemscores2txt=0;

% ###### Analysis details ################
% log.AnalysisType=' s4WithDeriv';  log.functionalscanprefix='s4wubf'; % NON-ANTS 
log.AnalysisType=' s4FullCardiacWithDerivAnts';  log.functionalscanprefix='s4ubf';
% request.onsetsmodel='m_c1_Contextall';
% request.onsetsmodel='m_c2_Contextonly';
% request.onsetsmodel='m_c3_ContextallNooutcome';
% request.onsetsmodel='m_c4_ContextallItempresonly';
% request.onsetsmodel='m_c5_ContextonlyItempresonly';
% request.onsetsmodel='m_c6_ContextallItemNomem'; 
request.onsetsmodel='m_c7_ContexteventItempresonly'; 
% request.onsetsmodel='m_ci1_ContextItem';
% request.onsetsmodel='m_ci2_ContextonlyItem';
% request.onsetsmodel='m_ci3_ContextItemNomotor';
% request.onsetsmodel='m_ci4_ContextItemNomotorNooutcome';
% request.onsetsmodel='m_ci5_ContextItempmod';
% request.onsetsmodel='m_ci6_ContextItempmodWithmotor';
% request.onsetsmodel='m_ci7_ContextMemcatItemNomotor';
% request.onsetsmodel='m_ci8_ContextMemcatItempresonly';
% request.onsetsmodel='m_ci9_ContextMemPcatItemNomotor';
% request.onsetsmodel='m_ci10_ContextEventItemNomotor';
% request.onsetsmodel='m_ci11_ContextEvent';
% request.onsetsmodel='m_ci12_ContextSimboxDisstickItemNomotor';
% request.onsetsmodel='m_ci13_ContextStickItemNomotor'; 
% request.onsetsmodel='m_i4_ItemContextevent';
% request.onsetsmodel='m_i5_ItemContextonly'; 
% log.AnalysisType=' s4FullCardiacAnts';  log.functionalscanprefix='s4ubf';  % FIR MODELS ##########  
% request.onsetsmodel='m_io1_ItemMem';   
% request.onsetsmodel='m_co1_Contextall';
% request.onsetsmodel='m_i1_Item';
% request.onsetsmodel='m_i2_ItemNooutcome';
% request.onsetsmodel='m_i3_ItemContextpresonly';
% request.onsetsmodel='m_i6_ItemMempmodContextpresonly';
request.memtype='Hit';   % Options: Hit\Surehit\Rem\Surerem\Roc 

for o1=1:1 % General settings and specifications
    
    % Load subjects
    addpath(where.where);
    log.w=load([where.data_brain filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
    % Scanning details
    addpath([where.where filesep '5b Scripts - Set up model' filesep '1a Specify Onsets models']);
    scan=load('i_scanningdetails.mat');
    scan.TRseconds=scan.TRms/1000;  % in seconds
    errorlog=cell(1,1); e=1;
    if strcmp(log.AnalysisType(3), log.functionalscanprefix(2))~=1
        disp('Functionals and Analysis type do not match - CHECK ?')
        disp(log.AnalysisType)
        disp(log.functionalscanprefix)
        input('Are these correct? Enter to Continue with these settings   ')
    end
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0;  disp('   Subset of subjects only:');  disp(log.subjects);  end
    disp(request);  disp(' '); disp(['Data location: ' where.data_brain])
    disp('=======================================================')
    
%     if 
    
end

%% (1) Format data (according to most basic design, Sim x Val x Mem)

for  o1=1:1 % Specifications for formatted data
    
    % CONTEXT Trial stats -----------------
    % Used for: Context, Outcome, Overall memory score
    col.c.SimDis=3; % Design variables
    col.c.RewNeu=4;
    col.c.Outcome_amount=7;
    col.c.Outcome_presented=29;
    col.c.Block=6;
    col.c.Trial=5;
    
    col.c.Item_Num(1)=11; % Event classifiers
    col.c.Item_Num(2)=12;
    col.c.Item_Num(3)=13;
    col.c.Item_EncAccuracy(1)=26;  
    col.c.Item_EncAccuracy(2)=27;
    col.c.Item_EncAccuracy(3)=28;
    col.c.Item_Outcome(1)=8;
    col.c.Item_Outcome(2)=9;
    col.c.Item_Outcome(3)=10;
    col.c.Item_Resp(1)=20;
    col.c.Item_Resp(2)=22;
    col.c.Item_Resp(3)=24;
    col.c.Item_Motor(1)=34;
    col.c.Item_Motor(2)=36;
    col.c.Item_Motor(3)=38;
    
    col.c.Item_Rocscore(1)=44; % Memory classifiers
    col.c.Item_Rocscore(2)=45;
    col.c.Item_Rocscore(3)=46;
    col.c.Memall_Hitscore=47;
    col.c.Memall_Surehitscore=48;
    col.c.Memall_Remscore=49;
    col.c.Memall_Sureremscore=50;
    col.c.Memall_Rocscore=51;
    col.c.Memall_Hitcat=52;
    
    col.c.Context_Onset=32; % Timing - Event onsets
    col.c.Item_Onset(1)=33;  
    col.c.Item_Onset(2)=35;
    col.c.Item_Onset(3)=37;
    col.c.Outcome_Onset=39;
    col.c.Endtrial_Onset=40;
    col.c.Context_Duration=41; % Event durations
    col.c.ContextOnly_Duration=42;
    col.c.Outcome_Duration=43;
    
    % ITEM trialstats (re-compiled) ---------------------
    % Trialstats used for: Item, Motor, Error
    col.i.SimDis=1;
    col.i.RewNeu=2; 
   
    col.i.Block=3;
    col.i.Trial=4;    
    col.i.ItemNo=5;
    
    col.i.Item_Resp=6;
    col.i.Item_Accuracy=7;
    col.i.Outcome_presented=8;
    col.i.Item_OutcomeOwn=9;
    col.i.Item_OutcomeTrial=10;
    
    col.i.Item_Onset=11;
    col.i.Item_Duration=12;
    col.i.Item_Error=13;
    col.i.Item_MotorOnset=14;
    col.i.Item_MotorOK=15;
   
    col.i.Mem_Roc=16;  % Memory scores
    col.i.Mem_HitOr=17;
    col.i.Mem_SurehitOr=18;
    col.i.Mem_RemOr=19;
    col.i.Mem_SureremOr=20;
    
    % Memory test specs ----------
    col.m.OldNew=3;
    col.m.HitMiss=12;
    col.m.RemKnow=14;
    col.m.SureGuess= 16;
    col.m.Roc=30;
    col.m.Item=2;
    
end


if request.simplifydata==1;
    disp('Formatting behavioural data ------------------------')
    subjdata=cell(log.n_subjs,1); w.block2offset=cell(log.n_subjs,3);
    for s=1:log.n_subjs
        disp(['Subject  ' num2str(s) '   (' log.subjects{s} ')'])
        subjdata{s,1}=log.subjects{s}; w.block2offset{s,3}=log.subjects{s};

        for o1=1:1 % How much to offset b2 by? 
            ws.wherefuncs=[where.data_brain filesep log.subjects{s} filesep '1 Preprocessed' filesep 'Func_b1'];
            
            ws.wherefuncs=['F:\1 Context Mem study\1 All data\' log.subjects{s} '\1 Preprocessed\Func_b1'];
            f=spm_select('List', ws.wherefuncs, ['^s4wubf.*img']);
            
%             f=spm_select('List', ws.wherefuncs, ['^' log.functionalscanprefix '.*img']);
            
            if isempty(f)==1
                error('Could not find any functionals to (to offset 2nd run!')
            elseif size(f,1)>270
                error('Implausibly high no. of functionals in 1st run!')
            end
            w.block2offset{s,1}=scan.TRseconds*(size(f,1)); % Dummies deleted already
            w.block2offset{s,2}=size(f,1);
        end
        
        % Load details
        ws.where=[where.data_beh filesep log.subjects{s} filesep];
        ws.w1=load([ws.where log.subjects{s} '_file_2encodingFMRI_b1.mat']); ws.e{1}=ws.w1.encoding.trialstats; ws.cogentstart{1}=ws.w1.encoding.times.start/1000;
        ws.w2=load([ws.where log.subjects{s} '_file_2encodingFMRI_b2.mat']); ws.e{2}=ws.w2.encoding.trialstats; ws.cogentstart{2}=ws.w2.encoding.times.start/1000;
        ws.m=load([where.data_beh filesep log.subjects{s} filesep log.subjects{s} '_file_3memorytest.mat']); ws.m=ws.m.memtest.trialstats;
        ws.m=ws.m(ws.m(:,col.m.OldNew)==1,:);    
        
        % Adjust timings + Mark durations
        for b=1:2 
            ws.e{b}(:,[34 36 38])=ws.e{b}(:,[34 36 38])/1000;
            ws.e{b}(:,30:40)=ws.e{b}(:,30:40)-ws.cogentstart{b};
            if b==2 % Adjust for concatenation
                ws.e{b}(:,30:40)=ws.e{b}(:,30:40)+w.block2offset{s,1};
            end
            ws.e{b}(:,col.c.Block)=b;
            
            % Durations
            ws.e{b}(:,col.c.Context_Duration)=ws.e{b}(:, col.c.Endtrial_Onset)-ws.e{b}(:, col.c.Context_Onset);
            ws.e{b}(:,col.c.ContextOnly_Duration)=ws.e{b}(:, col.c.Item_Onset(1))-ws.e{b}(:, col.c.Context_Onset);
            ws.e{b}(:,col.c.Outcome_Duration)=ws.e{b}(:, col.c.Endtrial_Onset)-ws.e{b}(:, col.c.Outcome_Onset);
        end
                
        % Mark Memory
        ws.m(ws.m(:,col.m.HitMiss)==1 & ws.m(:,col.m.RemKnow)==1 & ws.m(:,col.m.SureGuess)==1, col.m.Roc)=6;
        ws.m(ws.m(:,col.m.HitMiss)==1 & ws.m(:,col.m.RemKnow)==1 & ws.m(:,col.m.SureGuess)==2, col.m.Roc)=5;
        ws.m(ws.m(:,col.m.HitMiss)==1 & ws.m(:,col.m.RemKnow)==2 & ws.m(:,col.m.SureGuess)==1, col.m.Roc)=4;
        ws.m(ws.m(:,col.m.HitMiss)==1 & ws.m(:,col.m.RemKnow)==2 & ws.m(:,col.m.SureGuess)==2, col.m.Roc)=3;
        ws.m(ws.m(:,col.m.HitMiss)==2 & ws.m(:,col.m.SureGuess)==2, col.m.Roc)=2;
        ws.m(ws.m(:,col.m.HitMiss)==2 & ws.m(:,col.m.SureGuess)==1, col.m.Roc)=1;
        for b=1:2 % Mark memory scores in encoding data  
            for t=1:size(ws.e{b},1) % Item memory scores 
                for i=1:3
                    ws.e{b}(t, col.c.Item_Rocscore(i))=ws.m(ws.m(:,col.m.Item)==ws.e{b}(t, 10+i), col.m.Roc);
                end
            end
            for t=1:size(ws.e{b},1) % Context memory scores
                  ws.e{b}(t,col.c.Memall_Hitscore)=sum(ws.e{b}(t,col.c.Item_Rocscore(1):col.c.Item_Rocscore(3))>2.5);
                  ws.e{b}(t,col.c.Memall_Surehitscore)=sum(ws.e{b}(t,col.c.Item_Rocscore(1):col.c.Item_Rocscore(3))==6)+sum(ws.e{b}(t,col.c.Item_Rocscore(1):col.c.Item_Rocscore(3))==4);
                  ws.e{b}(t,col.c.Memall_Remscore)=sum(ws.e{b}(t,col.c.Item_Rocscore(1):col.c.Item_Rocscore(3))>4.5);
                  ws.e{b}(t,col.c.Memall_Sureremscore)=sum(ws.e{b}(t,col.c.Item_Rocscore(1):col.c.Item_Rocscore(3))==6);    
                  ws.e{b}(t,col.c.Memall_Rocscore)=sum(ws.e{b}(t,col.c.Item_Rocscore(1):col.c.Item_Rocscore(3)));    
                  ws.e{b}(t,col.c.Memall_Hitcat)=(sum(ws.e{b}(t,col.c.Item_Rocscore(1:3))>2.5)>1.5)+1; % 1 if n hits<2, otherwise 2
            end
        end
        ws.c=vertcat(ws.e{1},ws.e{2}); 
        
        % Compile Item trialstats (all Item-specific events)
        ws.i=cell(3,1);
        for i=1:3
            ws.i{i}=nan*zeros(size(ws.c,1),19);
            
            % Context-specific features
            ws.i{i}(:,col.i.SimDis)=ws.c(:, col.c.SimDis);
            ws.i{i}(:,col.i.RewNeu)=ws.c(:, col.c.RewNeu);
            ws.i{i}(:,col.i.Block)=ws.c(:, col.c.Block);
            ws.i{i}(:,col.i.Trial)=ws.c(:, col.c.Trial);
            ws.i{i}(:,col.i.Outcome_presented)=ws.c(:, col.c.Outcome_presented);
            ws.i{i}(:,col.i.Item_OutcomeTrial)=ws.c(:, col.c.Outcome_amount);
            
            % Item-specific features
            ws.i{i}(:, col.i.ItemNo)=ws.c(:,col.c.Item_Num(i));
            ws.i{i}(:, col.i.Item_Resp)=ws.c(:,col.c.Item_Resp(i));
            ws.i{i}(:, col.i.Item_Accuracy)=ws.c(:,col.c.Item_EncAccuracy(i));
            ws.i{i}(:, col.i.Item_OutcomeOwn)=ws.c(:,col.c.Item_Outcome(i));
            ws.i{i}(:, col.i.Item_Onset)=ws.c(:,col.c.Item_Onset(i));
            ws.i{i}(:, col.i.Item_Duration)=ws.c(:,col.c.Item_Onset(i)+2)-ws.c(:,col.c.Item_Onset(i));
            ws.i{i}(:, col.i.Item_MotorOnset)=ws.c(:,col.c.Item_Onset(i));
            ws.i{i}(:, col.i.Item_MotorOK)=(isnan(ws.i{i}(:, col.i.Item_MotorOnset))==0);
            ws.i{i}(:,col.i.Mem_Roc)=ws.c(:,col.c.Item_Rocscore(i));
        end
        ws.i=vertcat(ws.i{1},ws.i{2},ws.i{3}); ws.i=sortrows(ws.i,col.i.Item_Onset);
        ws.i(:,col.i.Item_Error)=0;
        ws.i(isnan(ws.i(:,col.i.Item_Resp))==1, col.i.Item_Error)=1;
        
        ws.i(:,col.i.Mem_HitOr)=(ws.i(:,col.i.Mem_Roc)>2.5); % Item memory scores
        ws.i(:,col.i.Mem_SurehitOr)=ws.i(:,col.i.Mem_Roc)==4 | ws.i(:,col.i.Mem_Roc)==6;        
        ws.i(:,col.i.Mem_RemOr)=ws.i(:,col.i.Mem_Roc)>4.5;
        ws.i(:,col.i.Mem_SureremOr)=ws.i(:,col.i.Mem_Roc)==6;
        
        %
        subjdata{s,2}.context=sortrows(ws.c, [col.c.Block col.c.Trial]);
        subjdata{s,2}.item=sortrows(ws.i, [col.i.Block col.i.Item_Onset]);
        ws=[];
        
    end
    
%     save([where.data_beh filesep '(' date ') Simplified data.mat'], 'subjdata')
    disp('Block 2 offset  (& no. scans in 1st block) -------------------------------------')
    disp(w.block2offset)
else
    try
        w.s=load([where.data_beh filesep request.specifydata]);
        subjdata=w.s.subjdata;
    catch
        disp(['Error: Could not find specified dataset   - ' request.specifydata])
        input('Request module in script to format simplified data instead     ');
    end
end

input('Construct regressors?   ')

% Print behavioural context-mem scores to text?
if request.printcontextmemscores2txt
    d_cmem=cell(log.n_subjs,6);
    d_cmem{1,1}='Subject';
    d_cmem{1,2}='Overall_cmem';
    d_cmem{1,3}='cmem_SimR';
    d_cmem{1,4}='cmem_SimN';
    d_cmem{1,5}='cmem_DisR';
    d_cmem{1,6}='cmem_DisN';
    for s=1:log.n_subjs
        d_cmem{s+1,1}=log.subjects{s};
        ws.d=subjdata{s,2}.context;
%         
%         ws.simr=ws.d(ws.d(:,col.c
%         
%         subjdata{s,2}.context(:,col.c.Memall_Hitscore)
%         
        
    end
    

end


%% (2) Model specifications (Set up Regressors - Item, Contexts & Other regressors)
% subjdata: Col 1=subject, Col 2=Formatted Item & Context dat, 
%                Col 3= Onsets variable (Combined Context & Item)
disp('Creating regressors -------------------------------')
design.Memscore4context={request.memtype};   % Cell, with ordered pmods
                                                        %    For options, see 'col.c.Memall_' ** 'score'        
design.Memscore4item=request.memtype;      % String. For options, see 'col.i.Mem_' ** 'Or'   [Hit\Surehit\Rem\Surerem]

for s=1:log.n_subjs
    ws.c=subjdata{s,2}.context;
    ws.i=subjdata{s,2}.item;
    r=1;
    
    % Construct regressors of interest according to requested model
    eval(['[ ws, r] = ' request.onsetsmodel '( ws, r,  col, design.Memscore4context, design.Memscore4item);'])
   
    %
    subjdata{s,3}=ws.v;
    ws=[];
end

%% (3) Save onsets ------------------------------

minimumcellsize=1; % Minimum no. of items in each cell?

% Save
if request.saveonsets==1
    disp('Saving onsets to First Level folders')
    subjectsok=cell(log.n_subjs,2);
    for i=1:log.n_subjs
        s=find(strcmp(subjdata, log.subjects{i})==1); % Identify correct subject
        names=subjdata{s,3}.names;
        onsets=subjdata{s,3}.onsets;
        durations=subjdata{s,3}.durations;
        try
            pmod=subjdata{s,3}.pmod;
        catch
            pmod=[];
        end
        subjectsok{i,1}=log.subjects{s};
        subjectsok{i,2}=1;
        for j=1:size(onsets,2)
            if size(onsets{j},1)<minimumcellsize
                subjectsok{i,2}=0;
            end
        end
        
        
        
        if isdir([where.data_brain filesep subjdata{s,1} filesep '2 First level' log.AnalysisType])==0; mkdir([where.data_brain filesep subjdata{s,1} filesep '2 First level' log.AnalysisType]); end
        save([where.data_brain filesep subjdata{s,1} filesep '2 First level' log.AnalysisType filesep subjdata{s,1} '_onsets_' request.onsetsmodel '_' request.memtype '.mat'], 'names', 'onsets', 'durations', 'pmod');
    end
end

%% 

disp('Log: Onsets files OK for all subjects?')
disp(subjectsok)
save([where.data_brain filesep 'Onsetslog_' request.onsetsmodel '_' request.memtype '.mat'], 'subjectsok')

% diary off
try % Notify researcher
%     f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s6_Getonsets)'), ' ',1);
end
