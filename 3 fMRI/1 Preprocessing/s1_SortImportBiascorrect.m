% SortImportDeletedummiesBiascorrect
clear all;close all hidden; clc

% where.where='I:\1 fMRI analysis'; where.data_brain='C:\Users\eloh\Desktop\1 [Context-Memory] fMRI Data\1 MRI data';  % where.data_beh=[where.where filesep '2 Behavioural data'];
where.where='/Volumes/PENNYDISK/1 fMRI analysis'; where.data_brain='/Volumes/SANDISK/1 CONTEXT Brain data'; % where.data_beh=[where.where filesep '2 Behavioural data'];

% Requested analysis
process.sortimport=1;
process.deletedummyvolumes=1;
process.biascorrect=1;
% 
log.specificsubjects={'p09_CN'}; % BLANK to process all subjects

for o1=1:1 % General settings and specifications
    
    % Subjects
    addpath(where.where)
    log.w=load([where.data filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    scan=load('i_scanningdetails.mat');
    
    % Log
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' (' date ')  ' log.model])
    errorlog=cell(1,1); e=1;


    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp('Requested analysis:'); disp(process)
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location: ' where.data_brain])
    disp(' ');  input('Hit Enter to start      ')
    disp('=======================================================')
    
    spm fmri
end

%% STEP 1: Sort and import together

if process.sortimport==1
disp('############### STEP 1: SORT & IMPORT ###############')
for s=1:log.n_subjs
    try
        disp(['Subject ' num2str(s) ' (' log.subjects{s} ')  --------- '])
        ws.from=[where.data filesep log.datalog{s+1} filesep '0 Original'];
        ws.to=[where.data filesep log.datalog{s+1} filesep '1 Preprocessed'];
        cd(ws.from)
        mkdir(ws.to)
        % Functional block 1
        ws.type='Func_b1'; 
        mkdir([ws.to filesep ws.type]) 
        ws.runname=[log.datalog{s+1,3} '.' num2str(log.datalog{s+1,6})];
        Import_Archive([ws.from filesep ws.runname '.tar'], [ws.to filesep ws.type filesep])
        cd([ws.to filesep ws.type filesep ws.runname])
        ws.files=dir('fM*');
        for i=1:size(ws.files,1)
            movefile(ws.files(i).name, [ws.to filesep ws.type])
        end
        cd ..
        rmdir([ws.to filesep ws.type filesep ws.runname],'s');
         % Functional block 2
        ws.type='Func_b2';
        mkdir([ws.to filesep ws.type]) 
        ws.runname=[log.datalog{s+1,3} '.' num2str(log.datalog{s+1,7})]; 
        Import_Archive([ws.from filesep ws.runname '.tar'], [ws.to filesep ws.type filesep])
        cd([ws.to filesep ws.type filesep ws.runname])
        ws.files=dir('fM*');
        for i=1:size(ws.files,1)
            movefile(ws.files(i).name, [ws.to filesep ws.type])
        end
        cd ..
        rmdir([ws.to filesep ws.type filesep ws.runname],'s');
         % Functional Fieldmap
        ws.type='Func_Fieldmap';
        mkdir([ws.to filesep ws.type]) 
        ws.runname=[log.datalog{s+1,3} '.' num2str(log.datalog{s+1,5})]; % 1st run
        Import_Archive([ws.from filesep ws.runname '.tar'], [ws.to filesep ws.type filesep])
        cd([ws.to filesep ws.type filesep ws.runname])
        ws.files=dir('sM*');
        for i=1:size(ws.files,1)
            movefile(ws.files(i).name, [ws.to filesep ws.type])
        end
        cd ..
        rmdir([ws.to filesep ws.type filesep ws.runname],'s');
        ws.runname=[log.datalog{s+1,3} '.' num2str(log.datalog{s+1,5}+1)]; % 2nd run
        Import_Archive([ws.from filesep ws.runname '.tar'], [ws.to filesep ws.type filesep])
        cd([ws.to filesep ws.type filesep ws.runname])
        ws.files=dir('sM*');
        for i=1:size(ws.files,1)
            movefile(ws.files(i).name, [ws.to filesep ws.type])
        end
        cd ..
        rmdir([ws.to filesep ws.type filesep ws.runname],'s');
         % MPM - b1
        ws.type='MPM_b1';
        mkdir([ws.to filesep ws.type]) 
        ws.runname=[log.datalog{s+1,8} '.' num2str(log.datalog{s+1,9})];
        Import_Archive([ws.from filesep ws.runname '.tar'], [ws.to filesep ws.type filesep])
        cd([ws.to filesep ws.type filesep ws.runname])
        ws.files=dir('sM*');
        for i=1:size(ws.files,1)
            movefile(ws.files(i).name, [ws.to filesep ws.type])
        end
        cd ..
        rmdir([ws.to filesep ws.type filesep ws.runname],'s');
         % MPM - b0
        ws.type='MPM_b0';
        mkdir([ws.to filesep ws.type]) 
        ws.runname=[log.datalog{s+1,8} '.' num2str(log.datalog{s+1,9}+1)]; % 1st run
        Import_Archive([ws.from filesep ws.runname '.tar'], [ws.to filesep ws.type filesep])
        cd([ws.to filesep ws.type filesep ws.runname])
        ws.files=dir('sM*');
        for i=1:size(ws.files,1)
            movefile(ws.files(i).name, [ws.to filesep ws.type])
        end
        cd ..
        rmdir([ws.to filesep ws.type filesep ws.runname],'s');
        ws.runname=[log.datalog{s+1,8} '.' num2str(log.datalog{s+1,9}+2)]; % 2nd run
        Import_Archive([ws.from filesep ws.runname '.tar'], [ws.to filesep ws.type filesep])
        cd([ws.to filesep ws.type filesep ws.runname])
        ws.files=dir('sM*');
        for i=1:size(ws.files,1)
            movefile(ws.files(i).name, [ws.to filesep ws.type])
        end
        cd ..
        rmdir([ws.to filesep ws.type filesep ws.runname],'s');
         % MPM - MTw
        ws.type='MPM_MTw';
        mkdir([ws.to filesep ws.type]) 
        ws.runname=[log.datalog{s+1,8} '.' num2str(log.datalog{s+1,9}+3)];
        Import_Archive([ws.from filesep ws.runname '.tar'], [ws.to filesep ws.type filesep])
        cd([ws.to filesep ws.type filesep ws.runname])
        ws.files=dir('sM*');
        for i=1:size(ws.files,1)
            movefile(ws.files(i).name, [ws.to filesep ws.type])
        end
        cd ..
        rmdir([ws.to filesep ws.type filesep ws.runname],'s');
         % MPM - PDw
        ws.type='MPM_PDw';
        mkdir([ws.to filesep ws.type]) 
        ws.runname=[log.datalog{s+1,8} '.' num2str(log.datalog{s+1,9}+5)];
        Import_Archive([ws.from filesep ws.runname '.tar'], [ws.to filesep ws.type filesep])
        cd([ws.to filesep ws.type filesep ws.runname])
        ws.files=dir('sM*');
        for i=1:size(ws.files,1)
            movefile(ws.files(i).name, [ws.to filesep ws.type])
        end
        cd ..
        rmdir([ws.to filesep ws.type filesep ws.runname],'s');
         % MPM - T1w
        ws.type='MPM_T1w';
        mkdir([ws.to filesep ws.type]) 
        ws.runname=[log.datalog{s+1,8} '.' num2str(log.datalog{s+1,9}+7)];   
        Import_Archive([ws.from filesep ws.runname '.tar'], [ws.to filesep ws.type filesep])
        cd([ws.to filesep ws.type filesep ws.runname])
        ws.files=dir('sM*');
        for i=1:size(ws.files,1)
            movefile(ws.files(i).name, [ws.to filesep ws.type])
        end
        cd ..
        rmdir([ws.to filesep ws.type filesep ws.runname],'s');
         % HPC
        ws.type='HPC';
        mkdir([ws.to filesep ws.type]) 
        ws.runname=[log.datalog{s+1,10} '.' num2str(log.datalog{s+1,11})]; % 1st run
        Import_Archive([ws.from filesep ws.runname '.tar'], [ws.to filesep ws.type filesep])
        cd([ws.to filesep ws.type filesep ws.runname])
        ws.files=dir('sM*');
        for i=1:size(ws.files,1)
            movefile(ws.files(i).name, [ws.to filesep ws.type])
        end
        cd ..
        rmdir([ws.to filesep ws.type filesep ws.runname],'s');
        ws.runname=[log.datalog{s+1,10} '.' num2str(log.datalog{s+1,11}+1)]; % 2nd run
        Import_Archive([ws.from filesep ws.runname '.tar'], [ws.to filesep ws.type filesep])
        cd([ws.to filesep ws.type filesep ws.runname])
        ws.files=dir('sM*');
        for i=1:size(ws.files,1)
            movefile(ws.files(i).name, [ws.to filesep ws.type])
        end
        cd ..
        rmdir([ws.to filesep ws.type filesep ws.runname],'s');
        ws.runname=[log.datalog{s+1,10} '.' num2str(log.datalog{s+1,11}+2)]; % 3rd run
        Import_Archive([ws.from filesep ws.runname '.tar'], [ws.to filesep ws.type filesep])
        cd([ws.to filesep ws.type filesep ws.runname])
        ws.files=dir('sM*');
        for i=1:size(ws.files,1)
            movefile(ws.files(i).name, [ws.to filesep ws.type])
        end
        cd ..
        rmdir([ws.to filesep ws.type filesep ws.runname],'s');
        ws.runname=[log.datalog{s+1,10} '.' num2str(log.datalog{s+1,11}+3)]; % 4th run
        Import_Archive([ws.from filesep ws.runname '.tar'], [ws.to filesep ws.type filesep])
        cd([ws.to filesep ws.type filesep ws.runname])
        ws.files=dir('sM*');
        for i=1:size(ws.files,1)
            movefile(ws.files(i).name, [ws.to filesep ws.type])
        end
        cd ..
        rmdir([ws.to filesep ws.type filesep ws.runname],'s');
        %
        ws=[];
    catch
        errorlog{e,1}=['Failed: Sort & Import --- ' log.datalog{s+1,1}];
        e=e+1;
    end
end
end

%% STEP 2: Delete dummy volumes + adjust for Special cases

if process.deletedummyvolumes==1
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) ' (' log.subjects{s} ')  --------- '])
        for b=1:2
            wb.where=[where.data filesep log.subjects{s} filesep '1 Preprocessed' filesep 'Func_b' num2str(b) filesep];
            f=spm_select('List', wb.where, '^fMQ.*000001-01'); if isempty(f)==1; erorrlog{e,1}=['Error. Could not find dummy to delete  - ' log.subjects{s} ' -  b' num2str(b)]; end; e=e+1;
            for i=1:size(f,1)
                delete([wb.where f(i,:)]);
            end
            wb=[];
        end
    end
end

% SPECIAL CASE: p09 block 1, cut off at 224 scans 
% (Total no. scans: 224 - 1 dummy = 223 vol =223*3s=669s)
% Block 1 data is cut off at end of Trial 45, which ends 667.364s AFTER
% start time
if sum(strcmp(log.subjects, 'p09_CN'))==1
    errorlog{e}='NOTE: p09_CN ends at volume 224'; disp(errorlog{e}); e=e+1;
    wb.where=[where.data filesep 'p09_CN' filesep '1 Preprocessed' filesep];
    delete([wb.where 'Func_b1' filesep 'fMQ01089-0006-00225-000225-01.hdr'])
    delete([wb.where 'Func_b1' filesep 'fMQ01089-0006-00225-000225-01.img'])
    wb=[];
end

%% STEP 3: Bias correction (Functionals only)

if process.biascorrect==1
disp('############### STEP 2: BIAS CORRECTION ###############')
for s=1:log.n_subjs
    try
        disp(['Subject ' num2str(s) ' --------- '])
        ws.to=[where.data filesep log.datalog{s+1} filesep '1 Preprocessed' filesep];
        spm_biascorrect([ws.to filesep 'Func_b1']);
        spm_biascorrect([ws.to filesep 'Func_b2']);
    catch
        errorlog{e,1}=['Failed: Bias correction --- ' log.datalog{s+1,1}];
        e=e+1;
    end
end
end

%% END

disp('====================================')
w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp(' ')
disp('Analysis completed:')
disp(process)
disp(['No. of subjects: ' num2str(log.n_subjs)])
disp(' ')
disp('Error log?')
disp(errorlog)
disp(' ')
disp('====================================')

diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s1_SortImportBiascorrect)'), ' ',1);
catch
end
