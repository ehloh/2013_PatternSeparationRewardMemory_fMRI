% Delete specific files/folders for all subjects 
clear all;close all hidden; clc

% where.where='I:\1 fMRI analysis'; 
where.where='/Volumes/SANDISK/1 fMRI analysis';

% where.exp_folder='D:\1 [Context-Memory] fMRI Data'; 
where.exp_folder='/Users/EleanorL/Desktop/1 CONTEXT fmri data';
where.data_beh=[where.where '\2 Behavioural data'];where.data_brain=[where.exp_folder filesep '1 MRI data'];


% Request specific
log.specificsubjects={}; 

wherefrom='/Users/EleanorL/Dropbox/SCRIPPS/1 ContextMem fMRI/2 Behavioural data/'; 
whereto='/Users/EleanorL/Dropbox/SCRIPPS/2a ContextMem behaviour/2 Behavioural Data/'; 

% fol='2 First level s6WithDeriv';
% fol='2 First level s4FullCardiacWithDeriv';
% fol='2 First level s4WithDeriv';




for o1=1:1 % General settings and specifications
    
    % Load subjects
    log.w=load([where.data_brain filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
%     scan=load([where.where filesep 'i_scanningdetails.mat']);
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' ');     ('Requested analysis: ADHOC SCRIPT'); disp('See direct code for exact actions to be executed !! '); disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location: ' where.data_brain]); disp(' ');  disp(['Data location (behaviour): ' where.data_beh])
    disp(' ');  input('Hit Enter to start      ')
    disp('=======================================================')

end

%%

for s=1: log.n_subjs
    %%% INSERT COMMANDS HERE  ###########################
    % ##############################################
    disp([log.subjects{s} ' ------------------']);
    ws.subfol=[where.data_brain filesep log.subjects{s} filesep ];
%         eval('java.io.File(ws.old).renameTo(java.io.File(ws.new));')
    
%%
% try
    
copyfile([wherefrom log.subjects{s} filesep 'Spike_6cardiac'], [whereto log.subjects{s} filesep 'Spike_6cardiac'])
copyfile([wherefrom log.subjects{s} filesep 'Spike_10cardiac'], [whereto log.subjects{s} filesep 'Spike_10cardiac'])
    


% catch
%     disp('dd')
    
% end



%% Move VOI images over to PC FL

    
%         eval('java.io.File(ws.old).renameTo(java.io.File(ws.new));')
ws=[];

end
