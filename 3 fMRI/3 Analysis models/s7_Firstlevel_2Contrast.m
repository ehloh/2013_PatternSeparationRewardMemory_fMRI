% First level analysis
clear all;close all hidden; clc
where.where='D:\Dropbox\SCRIPPS\1 ContextMem fMRI';  where.data_brain='D:\1 [Context-Memory] fMRI Data\1 MRI data';  % where.data_beh=[where.where filesep '2 Behavioural data'];


% Requested analysis
log.specificsubjects={}; % BLANK to process all subjects
request.AnalysisType= ' s4FullCardiacWithDerivAnts';
request.ExecuteContrasts=1;
request.Ants_SaveCleanContrasts=0;
request.Ants_DeleteUnadjustedContrasts=0;
%
request.WeightDerivs=0; % Weight 2nd and 3rd derivatives for target conditions
request.DeleteExistingFLContrasts=1;

% Model details
% log.onsetsmodel='m_c1_Contextall';
% log.onsetsmodel='m_c2_Contextonly';
% log.onsetsmodel='m_c3_ContextallNooutcome';
% log.onsetsmodel='m_c4_ContextallItempresonly';
% % log.onsetsmodel='m_c5_ContextonlyItempresonly';
% log.onsetsmodel='m_c6_ContextallItemNomem'; 
log.onsetsmodel='m_c7_ContexteventItempresonly';
% log.onsetsmodel='m_ci1_ContextItem';
% log.onsetsmodel='m_ci2_ContextonlyItem';
% log.onsetsmodel='m_ci3_ContextItemNomotor';
% log.onsetsmodel='m_ci4_ContextItemNomotorNooutcome';
% log.onsetsmodel='m_ci5_ContextItempmod';
% log.onsetsmodel='m_ci6_ContextItempmodWithmotor';
% log.onsetsmodel='m_ci7_ContextMemcatItemNomotor';
% log.onsetsmodel='m_ci8_ContextMemcatItempresonly';
% log.onsetsmodel='m_ci9_ContextMemPcatItemNomotor';
% log.onsetsmodel='m_ci10_ContextEventItemNomotor';
% log.onsetsmodel='m_ci11_ContextEvent';
% log.onsetsmodel='m_ci12_ContextSimboxDisstickItemNomotor';
% log.onsetsmodel='m_ci13_ContextStickItemNomotor'; 
% log.onsetsmodel='m_i1_Item';
% log.onsetsmodel='m_i2_ItemNooutcome';
% log.onsetsmodel='m_i3_ItemContextpresonly';
% log.onsetsmodel='m_i4_ItemContextevent';
% log.onsetsmodel='m_i5_ItemContextonly'; 
% log.onsetsmodel='m_i6_ItemMempmodContextpresonly';
log.memtype='Hit';   % Options: Hit\Surehit\Rem\Surerem
log.model=[log.onsetsmodel '_' log.memtype];

for o1=1:1 % General settings and specifications
    
    % Subjects
    addpath(where.where); %  addpath([where filesep '6b Scripts
    log.w=load([where.data_brain filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    w.onsets=load([where.data_brain filesep 'Onsetslog_' log.model]); % Apply subject filter first according to validity of onsets file
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog,log.specificsubjects,vertcat({'Subject' 'onsets'}, w.onsets.subjectsok),'onsets');    
    
    % Log
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' (' date ')  ' log.model])
    errorlog=cell(1,1); e=1;
    request.contrasttable='i7_firstlevel_contrasts';
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp('Requested analysis:')
    disp(request);  disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location (brain): ' where.data_brain])
    disp(' '); disp(['Model: ' log.model])
    disp(' '); input('Hit Enter to start      ')
    disp('=======================================================')
end

%% STEP 1: Specifications for requested model?

for o1=1:1 % Specifications for requested model 
    
    % Requested onsets model?
    if sum(strcmp(log.onsetsmodel,{'m_c1_Contextall'; 'm_c3_ContextallNooutcome'; 'm_c4_ContextallItempresonly';'m_c5_ContextonlyItempresonly';'m_c7_ContexteventItempresonly'}))==1
        contrasts.contrasttable='ContextOnly';
    elseif sum(strcmp(log.onsetsmodel,{'m_i1_Item';'m_i2_ItemNooutcome';'m_i3_ItemContextpresonly'}))==1
        contrasts.contrasttable='ItemOnly'; 
    elseif sum(strcmp(log.onsetsmodel,{'m_ci4_ContextItemNomotorNooutcome'; 'm_ci2_ContextonlyItem'; 'm_ci13_ContextStickItemNomotor'; 'm_ci3_ContextItemNomotor';'m_ci12_ContextSimboxDisstickItemNomotor';'m_i4_ItemContextevent';'m_i5_ItemContextonly'}))==1
        contrasts.contrasttable='ContextItem';
        if sum(strcmp(log.onsetsmodel,{'m_i5_ItemContextonly';'m_i4_ItemContextevent'}))==1; disp('MANUALLY turn off contextmem contrasts (in ContextItem sheet)! Don''t forget to reverse this after'), input('Continue?'); , end
    elseif sum(strcmp(log.onsetsmodel,'m_ci5_ContextItempmod'))==1
        contrasts.contrasttable='ContextItempmod';
        elseif sum(strcmp(log.onsetsmodel,'m_i6_ItemMempmodContextpresonly'))==1
        contrasts.contrasttable='ItempmodOnly';
    elseif sum(strcmp(log.onsetsmodel,{'m_ci7_ContextMemcatItemNomotor';}))==1
        contrasts.contrasttable='ContextMemcatItem';
    elseif sum(strcmp(log.onsetsmodel,{'m_ci8_ContextMemcatItempresonly'}))==1
        contrasts.contrasttable='ContextMemcatOnly';
    elseif sum(strcmp(log.onsetsmodel,{'m_ci10_ContextEventItemNomotor';'m_ci11_ContextEvent'}))==1
        contrasts.contrasttable='ContextEventOnly';
    % Special cases ###### 
    elseif sum(strcmp(log.onsetsmodel,'m_ci1_ContextItem'))==1
        contrasts.contrasttable='ContextOnly'; % Item regressors cannot be contrasted, because Items correlated with Motor       
    else
        error('Error in Contrasts setup: Could not find requested model - it may not have yet been specified')
    end
    
    % Which memory score?
    switch log.memtype
        case 'Hit'
            contrasts.memYes='Hit';
            contrasts.memNo='Miss';
        case 'Surehit'
            contrasts.memYes='Surehit';
            contrasts.memNo='NotSurehit';
        case 'Rem'
            contrasts.memYes='Rem';
            contrasts.memNo='NotRem';
        case 'Surerem'
            contrasts.memYes='Surerem';
            contrasts.memNo='NotSurerem';
        case 'Roc'
%             if isempty(strfind(contrasts.contrasttable, 'pmod'))==1
%                 error('Error: Roc memory (pmods) is requested, but contrast table does not ask for parametric mods')
%             end
    end
    
    % Details for requested contrast table (excel file)
    switch contrasts.contrasttable
        case 'ContextOnly'
            col.contrastname=2;
            col.weightstart=3;
            col.num_conds=8;
            col.weightend=col.weightstart+col.num_conds-1;
            col.contrastnum=col.weightend+1;
            col.requested=col.contrastnum+1;
            col.contype=col.requested+1;
            row.conditionName=1;
            row.conditionType=2;
        case 'ItemOnly'
            col.contrastname=2;
            col.weightstart=3;
            col.num_conds=8;
            col.weightend=col.weightstart+col.num_conds-1;
            col.contrastnum=col.weightend+1;
            col.requested=col.contrastnum+1;
            col.contype=col.requested+1;
            row.conditionName=1;
            row.conditionType=2;
        case 'ItempmodOnly'
            col.contrastname=2;
            col.weightstart=3;
            col.num_conds=9;
            col.weightend=col.weightstart+col.num_conds-1;
            col.contrastnum=col.weightend+1;
            col.requested=col.contrastnum+1;
            col.contype=col.requested+1;
            row.conditionName=1;
            row.conditionType=2;
        case 'ContextItem'
            col.contrastname=2;
            col.weightstart=3;
            col.num_conds=16;
            col.weightend=col.weightstart+col.num_conds-1;
            col.contrastnum=col.weightend+1;
            col.requested=col.contrastnum+1;
            col.contype=col.requested+1;
            row.conditionName=1;
            row.conditionType=2;
        case 'ContextItempmod'
            col.contrastname=2;
            col.weightstart=3;
            col.num_conds=16;
            col.weightend=col.weightstart+col.num_conds-1;
            col.contrastnum=col.weightend+1;
            col.requested=col.contrastnum+1;
            row.conditionName=1;
            row.conditionType=2;
        case 'ContextMemcatItem'
            col.contrastname=2;
            col.weightstart=3;
            col.num_conds=16;
            col.weightend=col.weightstart+col.num_conds-1;
            col.contrastnum=col.weightend+1;
            col.requested=col.contrastnum+1;
            row.conditionName=1;
            row.conditionType=2;
        case 'ContextMemcatOnly'
            col.contrastname=2;
            col.weightstart=3;
            col.num_conds=8;
            col.weightend=col.weightstart+col.num_conds-1;
            col.contrastnum=col.weightend+1;
            col.requested=col.contrastnum+1;
            row.conditionName=1;
            row.conditionType=2;
            
        case 'ContextEventOnly'
            col.contrastname=2;
            col.weightstart=3;
            col.num_conds=4;
            col.weightend=col.weightstart+col.num_conds-1;
            col.contrastnum=col.weightend+1;
            col.requested=col.contrastnum+1;
            row.conditionName=1;
            row.conditionType=2;
        otherwise
            error('Error: Could not find details for the requested contrast table (i.e. Excel sheet)')
    end
end

% input('Continue to Set Up Contrasts?    ')

%% STEP 2: Set up + Run Contrasts

for o1=1:1 % Set up Requested Contrasts 
    [w.a,w.b,w.c]=xlsread([where.where filesep '5b Scripts - Set up model' filesep request.contrasttable '.xlsx'], contrasts.contrasttable);
    
    % Adjust names for memory type
    for i=col.weightstart:col.weightend % Regressor-condition names 
        
        % Item memory in Conditions (not pmods) --------------
        if strfind(w.c{row.conditionName,i}, '_m1')
            w.c{row.conditionName,i}=[w.c{row.conditionName,i}(1:length(w.c{row.conditionName,i})-2) contrasts.memYes];
        elseif strfind(w.c{row.conditionName,i}, '_m0')
            w.c{row.conditionName,i}=[w.c{row.conditionName,i}(1:length(w.c{row.conditionName,i})-2) contrasts.memNo];
            
        % Item memory in Pmods (not conditions) --------------
        elseif strfind(w.c{row.conditionName,i}, 'ItemxMem')
            w.c{row.conditionName,i}=[w.c{row.conditionName,i} '_' log.memtype];
            
        % Context memory score (in pmods) --------------
        elseif strfind(w.c{row.conditionName,i}, 'xCMem') 
            w.c{row.conditionName,i}=[w.c{row.conditionName,i} '_' log.memtype];
            
        % Event memory score (in pmods, combined item & context conditions) --------------
        elseif strfind(w.c{row.conditionName,i}, 'xEMem') 
            w.c{row.conditionName,i}=[w.c{row.conditionName,i} '_' log.memtype];
        end
        
    end
    for i=row.conditionType+1:row.conditionType+1+col.weightend-col.weightstart % Contrast names (Aim is to change memtypes for identity matrix only)
        
        % Item memory in Conditions (not pmods) --------------
        if strfind(w.c{i, col.contrastname}, '_m1')
            w.c{i, col.contrastname}=[w.c{i, col.contrastname}(1:length(w.c{i, col.contrastname})-2) contrasts.memYes];
        elseif strfind(w.c{i, col.contrastname}, '_m0')
            w.c{i, col.contrastname}=[w.c{i, col.contrastname}(1:length(w.c{i, col.contrastname})-2) contrasts.memNo];

        % Item memory in Pmods (not conditions) --------------
         elseif strfind(w.c{i, col.contrastname},  'ItemxMem')
            w.c{i, col.contrastname}=[w.c{i, col.contrastname} '_' log.memtype];

        % Context memory score (in pmods) --------------
        elseif strfind(w.c{i, col.contrastname}, 'xCMem') 
            w.c{i, col.contrastname}=[w.c{i, col.contrastname} '_' log.memtype];
        end
        
    end
    
    % Compile details for requested contrasts
    for i=col.weightstart:col.weightend % Condition details
        contrasts.condnames{i-col.weightstart+1,1}=cellstr(w.c(row.conditionName,i));
        contrasts.condtype(i-col.weightstart+1,1)=w.c(row.conditionType,i);
    end
    c=1;
    for i=row.conditionName+1:size(w.c,1) % Construct details for requested contrasts
        if w.c{i,col.requested}==1
            contrasts.contrastname{c,1}=w.c{i,col.contrastname};
            contrasts.contype(c,1)=w.c{i, col.contype};
            for j=col.weightstart:col.weightend
                contrasts.weights(c,j-col.weightstart+1)=cell2mat(w.c(i,j));
            end
            c=c+1;
        end
    end
    
    % If weighting derivatives as well, change contrast names!
    if request.WeightDerivs
        contrasts.contrastname=cellfun(@(x)[x '_ad'], contrasts.contrastname,'UniformOutput',0);
    end
end

disp('Names of contrasts to run:'),  disp(contrasts.contrastname), input('Execute?'); 


% Execute contrasts ---------
if request.ExecuteContrasts
    input('Execute contrasts?      ')
    for s=1:log.n_subjs
        %     try
        disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
        wb.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' request.AnalysisType filesep];
        
        % Organize folder for this model (Rename from Estimated to Contrasted)
        wb.wheremodel=[wb.where log.model  ' Contrasted'  filesep];
        wb.oldname=[wb.where filesep log.model ' Estimated'];
        if isdir([wb.where filesep log.model ' Contrasted'])==1
            disp('Found contrasted folder ok. Assumed correct.'); if s==1; input('Continue?'); end
        else
            eval('java.io.File(wb.oldname).renameTo(java.io.File(wb.wheremodel));')
        end
        
        % Grab details about this model
        f   = spm_select('List', wb.wheremodel, 'SPM.mat');
        wb.spm  =load([wb.wheremodel filesep f]);
        wb.regnamelist=wb.spm.SPM.xX.name';
        
        % Specify requested contrasts
        disp('Specifying Contrasts ---------------------------');
        matlabbatch{1}.spm.stats.con.spmmat = cellstr([wb.wheremodel f]); c=1;
        if request.DeleteExistingFLContrasts
            if s==1; input('Requested deletion of all existing FL contrasts. Continue?   '); end
            matlabbatch{1}.spm.stats.con.delete = 1;
        else; matlabbatch{1}.spm.stats.con.delete = 0;
        end
        for j=1:length(contrasts.contrastname)
            if mean(abs(contrasts.weights(j,:)))~=0
                
                % Which regressors to actively weight?
                wc.condnums=find(contrasts.weights(j,:)~=0);
                for i=1:length(wc.condnums)
                    wc.regnames{i,1}=char(contrasts.condnames{wc.condnums(i)});
                    wc.regnames{i,2}=contrasts.condtype{wc.condnums(i)};
                end
                [ wc.regdetails ] = f_GetRegressorNums(wb.regnamelist, wc.regnames );
                
                % Compile weights for this contrast
                wc.weights=zeros(1,length(wb.regnamelist));
                for i=1:size(wc.regdetails ,1)
                    for k=1:size(wc.regdetails{i,3},1)
                        wc.weights(wc.regdetails{i,2}(k))=contrasts.weights(c,wc.condnums(i));
                    end
                end
                if request.WeightDerivs
                    if s==1 && j==1; input('Weighting all 3 derivatives for each conditions. Does this FL model even have derivs? OK to continue   '); end
                    wc.whichfirstderivs=find(wc.weights);
                    for k=1:length(wc.whichfirstderivs)
                        wc.weights(wc.whichfirstderivs(k)+1)=wc.weights(wc.whichfirstderivs(k));
                        wc.weights(wc.whichfirstderivs(k)+2)=wc.weights(wc.whichfirstderivs(k));
                    end
                end
                
                switch contrasts.contype(j)
                    case 1 % t contrast
                        matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = contrasts.contrastname{j};
                        matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
                        matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec=wc.weights;
                    case 2 % f contrast 
                        matlabbatch{1}.spm.stats.con.consess{c}.fcon.name=contrasts.contrastname{j};
                        matlabbatch{1}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
                        matlabbatch{1}.spm.stats.con.consess{c}.fcon.convec = {wc.weights};
                    otherwise
                        error('Invalid contrast type requested! F or T contrast?');
                end
                
                
                % Record
                contrast_weights{j,1}=contrasts.contrastname{j};
                contrast_weights{j,2}=wc.weights;
                contrast_weights{j,3}=wc.regdetails;
                
                wc=[]; c=c+1;
            end
        end
        
        disp('Running contrasts -----------------')
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
        matlabbatch=[];wb=[];
        %     catch
        %         errorlog{e}=['ERROR: Could not execute contrasts for subject  ' log.subjects{s}]; disp(errorlog{e}); e=e+1;
        %     end
        
    end
end




%% STEP 3: Save clean version of first-level contrasts (pre-Ants)

if request.Ants_SaveCleanContrasts
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
        ws.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' request.AnalysisType filesep];
       % 
        copyfile([ws.where log.model ' Contrasted'], [ws.where log.model ' PreAnts']);
        
        % Delete unadjusted contrasts (Not from pre-ants folder, but from 'Contrasted' folder
        if request.Ants_DeleteUnadjustedContrasts
            f=spm_select('List', [ws.where log.model ' Contrasted'], '^con_.*.');
            for i=1:size(f,1)
                delete([ws.where log.model ' Contrasted' filesep f(i,:)]);
            end
        end
        
        %
        ws=[];
    end
end

%% END

disp('=======================================================')
w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp(' ')
disp('Analysis completed (s7_Firstlevel_2Contrast):')
disp(request)
disp(['No. of subjects: ' num2str(log.n_subjs)])
disp(' ')
disp(log)
disp('Errors:')
disp(errorlog')
disp(' ')
disp('=======================================================')

diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s7_Firstlevel_2Contrasts)'), ' ',1);
end

