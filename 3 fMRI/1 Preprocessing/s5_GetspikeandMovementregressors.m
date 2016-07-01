% Get spike regressors + combine with movement parameters + Prep for 1st level
clear all;close all hidden; clc

where.where='I:\1 fMRI analysis';  where.data_brain='D:\1 [Context-Memory] fMRI Data\1 MRI data'; where.data_beh=[where.where filesep '2 Behavioural data'];
% where.where='/Volumes/PENNYDISK/1 fMRI analysis'; where.data_brain='/Volumes/SANDISK/1 CONTEXT Brain data';% where.data_beh=[where.where filesep '2 Behavioural data'];

% Requested analysis
log.AnalysisType=' s6FullCardiacWithDeriv';
log.CardiacType='_10cardiac';
log.specificsubjects={}; 

% Request procedures
request.spikelog2load=[]; % '(12-Mar-2013) spikelog.mat'; % Empty to extract spike
request.format_movementNphysio_regressors=1;

for o1=1:1 % General settings and specifications
    
    % Subjects
    addpath(where.where); addpath([where.where filesep '5a Scripts - Preprocessing'])
    log.w=load([where.data_brain filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
    % Add paths
    [spm_path,name]=fileparts(which('spm'));
    physiopath=sprintf('%s%s%s',spm_path,filesep,'toolbox',filesep,'physio');
    addpath(physiopath); sonpath=sprintf('%s%s%s%s%s',spm_path,filesep,'toolbox',filesep,'physio',filesep,'son'); addpath(sonpath);
    
    % Log
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' (' date ')'])
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('=======================================================')
    w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    disp('Requested analysis:')
    disp(request)
    if isempty(log.specificsubjects)==0
        disp('   Subset of subjects only:')
        disp(log.specificsubjects)
    end
    disp(' ')
    disp(['Data location (brain): ' where.data_brain])
    disp(['Data location (behaviour): ' where.data_beh])
    input('Hit enter to start                 ');
    disp('=======================================================')
    
end

%% STEP 1: Get Spike ------------------

if isempty(request.spikelog2load)==1
    disp(' ############## GET SPIKE ##############')
    spikelog=cell(log.n_subjs,4);
    for o1=1:1 % Scanning details ------------------------------------------
        w.scan=load('i_scanningdetails.mat');
        nslices=w.scan.nSlicesPerVol;  % Number of slices in volume
        ndummies=w.scan.nDummyVols;  % Number of scans excluded from analysis
        TRms=w.scan.TRms;     % TR in ms
        TR=TRms/1000;       % Slice TR in secs
        nsessions=1; % Number of scanning sessions in the file
        slicenum=ndummies/2;   % Slice number to time-lock regressors to
        % The above slice number can be determined from
        % data converted to nifti format. By default, slices
        % will be numbered from bottom to top but the acquisition
        % order can be ascending, descending or interleaved.
        % If slice order is descending or interleaved, the slice number
        % must be adjusted to represent the time at which the slice of
        % interest was acquired:
        %
        % For 3D acquisition sequences: Choose centre slice
        sliceorder='ascending'; % Ascending by default (verify with physics that this applies for this sequence: DONE)
        slicenum=get_slicenum(slicenum,nslices,sliceorder);
        
        % Channels
        channel.scanner=1;
        channel.cardiacTTL=2;
        channel.cardiacQRS=[];
        channel.respiration=5;
        %
    end
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------------------'])
        wb.where=[where.data_beh filesep log.subjects{s} filesep];
        spikelog{s,1}=log.subjects{s};
        for b=1:2  %
            disp(['Block ' num2str(b) ' -------' ])
            % Identify physio files
            f=spm_select('List', wb.where, [ '^' log.subjects{s} '.*.' num2str(b) '.smr']);
            physiofile=[wb.where f];
            slicenum=get_slicenum(slicenum,nslices,sliceorder); % Need slice information
            if ~isempty(f)
                disp(['Subject ' num2str(s) ' (' log.subjects{s}  ')  -   b' num2str(b) '  ---------------------------------------------------'])
                spikelog{s,b+1}=1;
                % #######################################################
                %
                % [From here, script is adapted from Physics' script 'physio_script.m'. See physics wiki for details]
                % The channel numbers must be assigned as they have been in spike. Unused channels should be set to empty using [];
                % The channel numbers can be checked using the routines  show_channels and check_channels as demonstrated below.
                % Once the channels have been set correctly, they should stay the same when using spike with the same set-up and configuration file.
                if s==1 && b==1
                    show_channels(physiofile);
                    check_channels(physiofile,channel.scanner,channel.cardiacTTL,channel.cardiacQRS,channel.respiration);
                else
                    disp('Channel allocations are assumed to be the same as 1st subject, 1st block')
                end
                
                % Call the main routine for calculating physio regressors. (Note: cardiacqrs calculation is disabled)
                [cardiac,cardiacqrs,respire,rvt]=make_physio_regressors(physiofile,nslices,ndummies,TR,...
                    slicenum,nsessions,channel.scanner,channel.cardiacTTL,channel.cardiacQRS,channel.respiration);
                
                % Save a record of parameters used for the regressors
                save([physiofile(1:end-4) '_physioparams' log.CardiacType], 'physiofile', 'nslices', 'ndummies', 'TRms','slicenum','nsessions','sliceorder');
                
                % ############## Construct physiological regressors ###############################
                % For each session, put regressors in a matrix called R.  Each individual set of regressors are saved and also all regressors are saved with the name 'physiofile_R_session%d'.
                % These files can be loaded into an SPM design matrix using the 'Multiple Regressors' option.                
                % NB motion parameters can also be concatenated with the physio regressors and saved as a set of regressors called R (see below for example)
                for sessnum=1:nsessions
                    R=[];
                    if ~isempty(cardiac{sessnum}) && ~isempty(cardiac{sessnum})
                        cardiac_sess = cardiac{sessnum};
                        filename = sprintf('%s_cardiac_session%d',spm_str_manip(physiofile,'r'),sessnum);
                        save(filename, 'cardiac_sess');
                        R=cat(2,R,cardiac{sessnum}); errorlog{e}=('NOTE Non-Standard setup: Including all cardiac parameters!'); disp(errorlog{e}); e=e+1;
                        % R=cat(2,R,cardiac{sessnum}(:,1:6)); % Original: first 6 only.
                    end
                    if ~isempty(cardiacqrs{sessnum}) && ~isempty(cardiacqrs{sessnum})
                        cardiacqrs_sess = cardiacqrs{sessnum};
                        filename = sprintf('%s_cardiacqrs_session%d',spm_str_manip(physiofile,'r'),sessnum);
                        save(filename, 'cardiacqrs_sess');
                        R=cat(2,R,cardiacqrs{sessnum}(:,1:6));
                    end
                    if ~isempty(respire) && ~isempty(respire{sessnum})
                        respire_sess = respire{sessnum};
                        filename = sprintf('%s_respire_session%d',spm_str_manip(physiofile,'r'),sessnum);
                        save(filename, 'respire_sess');
                        R=cat(2,R,respire{sessnum}(:,1:6));
                    end
                    if ~isempty(rvt) && ~isempty(rvt{sessnum})
                        rvt_sess = rvt{sessnum};
                        filename = sprintf('%s_rvt_session%d',spm_str_manip(physiofile,'r'),sessnum);
                        save(filename,'rvt_sess');
                        R=cat(2,R,rvt{sessnum}(:,1:size(rvt{sessnum},2)));
                    end
                    nfiles=size(R,1);
                    
                    % Save R for all physio only
                    if nfiles>0
                        oR=R;
%                         Rname = sprintf('%s_R_session%d',spm_str_manip(physiofile,'r'),sessnum);
                        Rname=[physiofile(1:end-4) log.CardiacType '_R_session' num2str(sessnum)];
                        R=R-repmat(mean(R),nfiles,1);
                        if isempty(log.CardiacType)==0
                            note=['Note: Non-standard cardiac type (' log.CardiacType ')'];
                            save(Rname, 'R','note');
                        else
                            save(Rname, 'R');
                        end
                    end
                    %    If required, also concatenate with motion params, e.g.
                    %  load rp_example.txt
                    %   RP=rp_example;
                    %   R=cat(2,oR,RP);
                    %   R=R-repmat(mean(R),nfiles,1);
                    %   Rname = sprintf('%s_R_session%d',spm_str_manip(physiofile,'r'),sessnum);
                    %   save(Rname, 'R');
                end
            else
                disp(['Subject ' num2str(s) ' (' log.subjects{s}  ')        -   b' num2str(b) '     [missing]   ------------------------------------------'])
                spikelog{s,b+1}=0;
            end
            % --------------------- END OF PASTED SCRIPT --------------------------------
            % ############################################
        end
        spikelog{s,4}=floor((spikelog{s,3}+spikelog{s,2})/2);
        save([where.data_beh filesep '(' date ') spikelog'], 'spikelog')
    end
else
    load([where.data_beh filesep request.spikelog2load]); % variable 'spikelog'
end

%% STEP 2: Format Spike & Movement parameters
% CORRECT HERE: Correct for subjects with irregular scan #s
% p09_CN: Block 1, scan 225 deleted to deal with re-timing (scanner shut off early, volumes deleted for last trial)
% p04_JL is missing SPIKE for block 2 - exclude entirely

if request.format_movementNphysio_regressors==1
    physiomovereg=cell(log.n_subjs,   4); % Col 1=Subject, 2=Block 1 reg, 3= Block 2 reg, 4=Concatenated Block 1 & 2;
    printregressorsok=cell(log.n_subjs,2);
    for i=1:log.n_subjs
        s=find(strcmp(spikelog, log.subjects{i})==1);
        disp(['Subject ' num2str(s) '   (' spikelog{s,1} ') -------------------------------'])
        physiomovereg{s,1}=spikelog{s,1}; 
        for b=1:2
            try
                wb.where=[where.data_brain filesep spikelog{s,1} filesep '1 Preprocessed' filesep 'Func_b' num2str(b) filesep];
                % Load Motion -----------------------------------------------------------------------------------------------
                f=spm_select('List', wb.where, '^rp.*txt$');
                wb.motion=load([wb.where f(1,:)]); if isempty(f)==1; errorlog{e}=(['ERROR: Could not find motion parameters for   ' spikelog{s,1} '  block ' num2str(b)]); disp(errorlog{e}); e=e+1; end
                % Load Respiration  ------------------------------------------------------------------------------------------
                if spikelog{s,4}==1
                    wb.wbreathe=[where.data_beh filesep spikelog{s,1} filesep];
                    f=spm_select('List', wb.wbreathe, [ num2str(b) log.CardiacType '_R_session1.mat$']);
                    wb.resp=load([wb.wbreathe f]);
                else
                    wb.resp.R=[];
                    errorlog{e}=['NOTE: Respiration regressors not included for subject  ' spikelog{s,1} '  both blocks']; disp(errorlog{e}); e=e+1;
                end
                
                % CHECK: Spike length matches no. of scans? ---------------------------------------------------------------------------------------------
                if  strcmp(spikelog{s,1},'p09_CN')==1 && b==1 % specific correction for p09_CN block 1
                    wb.resp.R=wb.resp.R(1:size(wb.motion,1),:);
                    errorlog{e}='NOTE: p09_CN block 1''s SPIKE paramters are cut short to match no. of scans'; disp(errorlog{e}); e=e+1;
                elseif spikelog{s,4}==1 && size(wb.resp.R,1) ~= size(wb.motion,1) % Unexpected mismatch?
                    wb.resp.R=wb.resp.R(1:size(wb.motion,1),:);
                    errorlog{e}=['ERROR: No. scans according to SPIKE & motion parameters does not match up     (' spikelog{s,1} '  block  ' num2str(b) ')']; disp(errorlog{e}); e=e+1;
                end
                % Save -------------------------------------------------------------------------------------------------
                wb.R=horzcat(wb.resp.R,wb.motion); % Respiration first (14 col), then Motion (6 col)
                wb.R=wb.R-repmat(mean(wb.R),size(wb.R,1),1); R=wb.R;
%                 save([wb.where spikelog{s,1} '_b' num2str(b) '_reg_physiomovement'], 'R');
                physiomovereg{s,1+b}=R; 
                wb=[];
            catch
                errorlog{e}=['ERROR: Could not format Spike+Motion regressors for  ' spikelog{s,1} '  block  ' num2str(b)]; disp(errorlog{e}); e=e+1;
            end
        end
        physiomovereg{s,4}=vertcat(horzcat(physiomovereg{s,2}, ones(size(physiomovereg{s,2},1),1)), horzcat(physiomovereg{s,3}, zeros(size(physiomovereg{s,3},1),1)));  % Append block regressor
        ws.preprocfolder=[where.data_brain filesep spikelog{s,1} filesep '1 Preprocessed'];
        [w.printok]=print2txt(ws.preprocfolder, [spikelog{s,1} '_reg_physiomovement' log.CardiacType], physiomovereg{s,4});
        printregressorsok{s,1}=w.printok.fileprint; printregressorsok{s,2}=w.printok.filemoved;
        if spikelog{s,4}==1 &&  size(physiomovereg{s,4},2)~=10+6+2+6+1 % 10 heart + 6 resp + 2 rvt + 6 motion+1 block column
            errorlog{e}=['ERROR: No. columns for regressor txt file is wrong     (' spikelog{s,1} ')']; disp(errorlog{e}); e=e+1; 
        end
    end
end

%% STEP 3: Organize spike files (behavioural folder)

request.organizespikefiles=1; % This needs only be done once
if request.organizespikefiles==1
    disp('STEP 3: Organizing SPIKE files -------------------------------------------')
    w.trashspikefiles={'respire_session1', 'rvt_session1', 'cardiac_session1', ['physioparams' log.CardiacType], [log.CardiacType(2:end) '_R_session1']};
    if isempty(request.spikelog2load)==0
        load([where.data_beh filesep request.spikelog2load]);
    end
    for s=1: log.n_subjs
        disp(['Subject ' num2str(s) '  (' spikelog{s,1} ') ----------------------'])
        wb.where=[where.data_beh filesep spikelog{s,1} filesep];
        ws.spikefol=[wb.where 'Spike' log.CardiacType];
        if isdir(ws.spikefol)==0; mkdir(ws.spikefol); end
        for b=1:2
            try
                for i=1:length(w.trashspikefiles)
                    movefile([wb.where spikelog{s,1} '_b' num2str(b) '_' w.trashspikefiles{i} '.mat'],[ws.spikefol filesep spikelog{s,1} '_b' num2str(b) '_' w.trashspikefiles{i} '.mat']);
                end
            catch
                errorlog{e}=['ERROR in organizing SPIKE files: Could not find files for ' spikelog{s,1} '  block ' num2str(b) '  --   ' w.trashspikefiles{i}]; disp(errorlog{e}); e=e+1;
            end
        end
    end
end

%%  End

% Spike log:
% Col 1=Subject, Col 2= Spike for block 1 (1=Yes,0=Np), Col 3- Spike for
% block 2, col 4=spike included in regressors?
disp('------------------------------------------------------------------------')
disp('SPIKE regressors successfully constructed [block 1, block 2; overall to include spike/not]')
if isempty(request.spikelog2load)==1; disp('SPIKELOG:'); disp(spikelog); end
disp('Error Log:');  errorlog{:}
disp('Note: spikelog is saved in same location as behavioural data, regressor files are saved with brain data')
disp('------------------------------------------------------------------------')

diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s5_GetspikeandMovementregressors)'), ' ',1);
catch
end
