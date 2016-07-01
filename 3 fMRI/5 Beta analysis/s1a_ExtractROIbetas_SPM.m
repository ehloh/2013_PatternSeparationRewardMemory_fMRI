% Extract ROI betas using SPM functions (subject-level). Voxel dimensions of ROIs and fxnal data must match. 
clear all;close all hidden; clc, log.betasuffix =[]; 
for o1=1:1 % Documentation (Read me)
    % Some variables in this script have been renamed. If bugs in future
    % scripts, see if changing the names here will help 
    %
    % instruc = 
    %                  log: 
    %         rois_extract: instruc.roi_files
    %         subcon_names: instruc.FLcons
    %        subcon_images: instruc.FLcon_imgs
    %         spec_connums: request.FLcon_nums
    %     req_subcon_names: request.FLcons
    % 
    % d_subroibetas --> d_betas
    %
% INSTRUCTIONS for using this script
%
%   (1) This script extracts ROIs from each subject x contrast image, and outputs 
%           it to the variable 'd_betas'. Betas are extracted for all
%           subjects included in the model (noted in output log.subjects)
%
%   (2) Outputs of script: 
%             - 'd_betas' is a cell array, with cells containing betas for all subjects (in order),
%                 corresponding to each ROI (row + 1) and each contrast (col + 1). Headers (ROI & Contrast)
%                 are included, thus displacing the beta vectors by 1. 
%             - 'd_vol.roi_masks':            Cell array, each cell holds the mask for one anatomical ROI
%             - 'd_vol.subcon':                 Cell array, each cell holding the full volume betas for each 
%                                                         subject (row) x contrast (col)
%             - 'instruc.FLcons':    Names of each (first-level, subjects')  contrast
%             - 'instruc.FLcon_imgs':   Contrast image name for each first-level contrast
%
% 	(2) Before running this script, ROIs must be defined within SPM. No need for importing via MarsBar. 
%
%   (3) Specify whether you would like ALL ROIs extracted, or only specific ones. 
%             - All ROIs: ROIs will be extracted for all images with the prefix 'roi', in the model's 
%                 2nd level results folder
%             - Specific ROIs: Only ROIs in the results folder's ROI folder will be extracted,
%                 again with prefix 'roi'
%   

% Other nifty ways of extracting betas
% spm_summarise('beta_1680.img',struct('def','sphere','xyz',[1 2 3]','spec',8),@mean)   %  Extracting from a sphere of 8mm at point 1,2,3
% spm_summarise('beta_1680.img','mask.img',@mean) % Extract mask.img from beta_1680.img (mean)
%   See also: spm_ROI
end

% Requested analysis
request.LoadSpecificROIs_Folder=1; % 1= Load only ROIs that are in the ROI folder, 0=Load all available ROIs for that results model
request.modeltype=1; % 1 = Par, 2= Flex

% Requested analysis
log.specificsubjects={}; % 'p01_CW';'p04_JL'};

% [UNIVAR - MODEL DETAILS] ############################
log.FLthread=' s4FullCardiacWithDerivAnts'; log.AntsType='_Landmarks5';log.memtype='Hit'; 
log.onsetsmodel='m_ci3_ContextItemNomotor';log.secondlevelmodel='cm_m1_2x2';
% log.onsetsmodel='m_ci2_ContextonlyItem'; log.secondlevelmodel='cm_m1_2x2';
% log.onsetsmodel='m_i4_ItemContextevent'; log.secondlevelmodel='iom_m1_2x2x2';
% log.onsetsmodel='m_i5_ItemContextonly'; log.secondlevelmodel='iom_m1_2x2x2';
% log.onsetsmodel='m_c4_ContextallItempresonly';log.secondlevelmodel='cm_m1_2x2';
% log.onsetsmodel='m_i3_ItemContextpresonly';log.secondlevelmodel='iom_m1_2x2x2';
% log.onsetsmodel='m_i6_ItemMempmodContextpresonly';log.secondlevelmodel='im_m1_2x2'; 
% log.onsetsmodel='m_ci13_ContextStickItemNomotor';  log.secondlevelmodel='cm_m1_2x2';
% log.onsetsmodel='m_c7_ContexteventItempresonly';log.secondlevelmodel='co_m1_2x2';

% [PPI comparison ] ############################
% where.branch_secondlevelwithinfirst='PPI comparison'; log.betasuffix =[]; 
% log.secondlevelmodel='coFam_sph_SNVTA_R2_cmSRvSN_psy_SR-SN';
% log.secondlevelmodel='coFam_sph_SNVTA_R1_cmSRvSN_psy_SR-SN';

% [PPI cell] ############################
    
seed='sph_SNVTA_R1_cmSRvSN'; 
log.betaname_prefix='DisN_';  log.betasuffix=['-' log.betaname_prefix(1:4)];
where.branch_secondlevelwithinfirst=['PPI cell\coFam_' seed '_psy_Cell\'];
log.secondlevelmodel=['coFam_' seed '_psy_' log.betaname_prefix(1:4)];


for o1=1:1 % General settings and specifications

    
     % Settings that don't change much
    where.where='D:\Dropbox\SCRIPPS\1 ContextMem fMRI' ; where.experiment_folder='F:\1 [Context-Memory] fMRI Data';  where.data_brain=[where.experiment_folder filesep '1 MRI data'];     where.beh='D:\Dropbox\SCRIPPS\2a ContextMem behaviour'; 
    addpath(where.where); 
    log.firstlevelmodel=[log.onsetsmodel '_' log.memtype log.AntsType];
    log.w=load([where.experiment_folder filesep '1 MRI data' filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    if strcmp(where.where(1),'/')==1; fs='/'; where.marsbar_anatomy='/Users/EleanorL/Documents/MATLAB/spm8/Anatomy masks/marsbar_imported'; else fs='\'; where.marsbar_anatomy='D:\My Documents\MATLAB\spm8\Anatomy masks\marsbar_imported'; end
%     where.apriorirois_marsbar=[where.experiment_folder filesep '3 Checks' filesep '2 A priori ROIs' filesep 'ROI analysis MarsBar'];
%     where.apriorirois=[where.experiment_folder filesep '3 Checks' filesep '3 VOI seeds'];
    
    % Subjects
    w.onsets=load([where.experiment_folder filesep '1 MRI data' filesep 'Onsetslog_' log.onsetsmodel '_' log.memtype ]); % This filter is only applied first to catch subjects with bad onsets files, where this has yet to be marked in excel sheet
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog,[],vertcat({'Subject' 'onsets'}, w.onsets.subjectsok),'onsets');
    if isempty(strfind(log.secondlevelmodel, 'Fam'))==0
        log.ppi=1;
        log.ppi_psy=[log.secondlevelmodel(strfind(log.secondlevelmodel, '_psy_')+5:end)];
        disp('Selecting subjects on the basis of PPI validity');        
        [w.s w.s1 log.koshertable]=xlsread([where.where filesep 'i8_subjectdataok_secondlevel.xlsx'],'PPIsingle');
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable, [log.firstlevelmodel ' ' log.secondlevelmodel(strfind(log.secondlevelmodel, 'Fam')+4:strfind(log.secondlevelmodel, '_psy_')-1)]);        
    else        
        log.ppi=0; where.branch_secondlevelwithinfirst=[]; 
        [w.s w.s1 log.koshertable]=xlsread([where.where filesep 'i8_subjectdataok_secondlevel.xlsx']);
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.firstlevelmodel);
    end
    
    % Model details (Full model, 1st & 2nd levels)
    if log.ppi==1, where.model=[where.experiment_folder filesep '2 Second level results' log.FLthread filesep log.firstlevelmodel filesep where.branch_secondlevelwithinfirst log.secondlevelmodel filesep];
    else, where.model=[where.experiment_folder filesep '2 Second level results' log.FLthread filesep log.firstlevelmodel filesep log.secondlevelmodel filesep];
    end
%     m=load([where.model 'SPM.mat']); m=m.SPM;
    if log.ppi==1, log.FLfolname=['2 First level' log.FLthread filesep log.firstlevelmodel ' Contrasted' filesep 'PPI ' log.secondlevelmodel ];
    else log.FLfolname=['2 First level' log.FLthread filesep log.firstlevelmodel ' Contrasted' ];  log.betasuffix=[]; 
    end
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location (brain): ' where.experiment_folder])
    disp(' '); disp('CHECK HERE: 1st and 2nd level models ---------'); disp(' ')
    disp(['             Model type:               ' log.FLthread])
    disp(['             First level model:      ' log.firstlevelmodel])
    disp('NOTE: If ROIs are not in same space as data, needs to be resliced!')
    input('Hit Enter to start      ')
    disp('=======================================================')
    
    
end
 
p=load([where.data_brain filesep log.subjects{1}  filesep log.FLfolname filesep 'SPM.mat']);p=p.SPM;


%% AD HOC SETTINGS * * * #################

instruc.roi_names=[];  % Empty to use roi actual names

% Manually define ############
where.roi_foldername=['ROI' filesep];
% where.roi_foldername=['ROI' filesep 'ELDrawn'];
% where.roi_foldername=['ROI' filesep 'i6 rois'];
where.roi_foldername=['ROI' filesep 'Anat' filesep]; 
where.roi_foldername=['ROI' filesep 'Anat Subfield' filesep]; 
% where.roi_foldername=['ROI' filesep 'ci3 rois' filesep]; 
% where.roi_foldername=['ROI' filesep 'ci3 rois' filesep 'peaks']; 
% where.roi_foldername=['ROI' filesep 'ci3 rois' filesep 'Second string']; 

%% (1) Extract anatomic ROIs and volume betas 

% Load all ROIs (n-ary or binary masks)
disp('Loading ROI images ########################')
where.model_res=[where.experiment_folder filesep '2 Second level results' log.FLthread filesep log.firstlevelmodel filesep where.branch_secondlevelwithinfirst filesep log.secondlevelmodel filesep];
if request.LoadSpecificROIs_Folder==1; where.model_resrois=[where.model_res where.roi_foldername filesep];  % Read only the ROIs in the ROI folder
else where.model_resrois=[where.model_res filesep];
end
instruc.roi_files=cellstr(spm_select('List', where.model_resrois, '.*.nii$'));   % Change format of files here!! 
if isempty(char(instruc.roi_files))==1;  instruc.roi_files=[]; end
d_vol.roi_masks=cell(size(instruc.roi_files,1),1);
for i=1:size(instruc.roi_files,1)
    d_vol.roi_masks{i}=spm_read_vols(spm_vol([where.model_resrois instruc.roi_files{i}]));
end

% 
% for i=1:size(instruc.roi_files,1)
%     disp([instruc.roi_files{i} '  ' num2str(sum(d_vol.roi_masks{i}(:)>0))])
% end

% Interface
disp('ROIs to extract:    ( !! CANNOT BE EMPTY !! )'); disp(char(instruc.roi_files)); input('OK?   '); disp('Names:'); disp(instruc.roi_names)

%% (2) Extract subjects' volume (contrast) betas & calculate mean

% Certain first-level contrast images only? ########### *** EDIT HERE *** ###########
disp('Available FL contrasts in this model:  ' ); disp([num2cell(1:size(p.xCon,2))' cellstr(char(p.xCon.name))])
request.FLcon_nums=[];
request.FLcon_nums=1;
if isempty(request.FLcon_nums)==0; disp('Requested contrasts:  '); disp(char(p.xCon(request.FLcon_nums).name));  else disp('Extracting betas from ALL FL contrasts');  end; input('Proceed?    '); 
    
for o1=1:1
    % Load all contrast images from all subjects (d_vol.subcon, row=sub, col=con image)
    disp('Loading subjects'' contrast images ########################')
    instruc.FLcons=cell(size(p.xCon,2),1); instruc.FLcon_imgs=cell(size(p.xCon,2),1);
    for i=1:length(instruc.FLcons) % Get names and available con images (sample 1st subject)
        instruc.FLcons{i}=p.xCon(i).name; instruc.FLcon_imgs{i}=p.xCon(i).Vcon.fname;
    end
  
    % Load first-level contrasts (see above for specification: extract from all or from some only?)
    if isempty(request.FLcon_nums)==0
        instruc.FLcons=instruc.FLcons(request.FLcon_nums);
        instruc.FLcon_imgs=instruc.FLcon_imgs(request.FLcon_nums);
    end
    
    % Extract
    d_vol.subcon=cell(log.n_subjs, length(instruc.FLcons));
    for s=1:log.n_subjs % Load volume betas for each contrast
        disp(['Subject ' num2str(s) '  -  ' log.subjects{s} ' --------------------------------'])
        ws.wheremod1st=[where.data_brain filesep log.subjects{s} filesep log.FLfolname filesep];
        disp(ws.wheremod1st)
        for c=1:length(instruc.FLcons)
            disp(['Contrast no. ' num2str(c) '  -  ' instruc.FLcons{c}]);
            d_vol.subcon{s,c}=spm_read_vols(spm_vol([ws.wheremod1st instruc.FLcon_imgs{c}]));
            
            
            
            
            %         v=spm_vol([ws.wheremod1st instruc.FLcon_imgs{c}]);
            %         d_vol.subcon{s,c}=spm_read_vols(v);
        end
    end
    
%     instruc
%     rm=d_vol.roi_masks{1}; 
%     sc= d_vol.subcon{s,c}~=0;
%     [rr cc vv]=find(d_vol.subcon{s,c}~=0); 
%     
%     
%     
%     
%     r
%     
%     error
    
    % Subject level (d_betas, col+1=contrast image, row+1=roi)
    disp('Calculating ROI betas for each ROI x Subject x Contrast #################')
    d_betas=cell(length(instruc.roi_files)+1, length(instruc.FLcons)+1);
    for r=1:length(d_vol.roi_masks)
        disp(['ROI no. ' num2str(r) ': ' instruc.roi_files{r} ' -------------------------'])
        d_betas{r+1,1}=instruc.roi_files{r};
        for c=1:length(instruc.FLcons)
            disp(['Contrast no. ' num2str(c) ': '  instruc.FLcons{c} '----'])
            if r==1; d_betas{1,c+1}=instruc.FLcons{c}; end % Title
            
            % Read each subject's ROI betas
            wc.s=nan*zeros(log.n_subjs,1);
            for s=1:log.n_subjs
                disp(['Subject ' num2str(s) ': ' log.subjects{s}])
%                 wc.s(s)=mean(d_vol.subcon{s,c}(d_vol.roi_masks{r}~=0));
                wc.s(s)=nanmean(d_vol.subcon{s,c}(d_vol.roi_masks{r}~=0));
            end
            d_betas{r+1, c+1}=wc.s;
            wc=[];
        end
    end
end

% Clear for workspace 
% d_vol=[];
% error


% d_vol.roi_masks{1}

%% (4) Export to txt (for SPSS?)

request.export2excel=1; request.calculateandprint_choicedifferencescores=0; if request.calculateandprint_choicedifferencescores; disp('Choice difference scores are requested!'); end

% Write to text file
if request.export2excel==1
    
    % ########### *** EDIT HERE *** ###########
    if isempty(instruc.roi_names); 
        % If crashing here, likely you have no ROIs
        instruc.roi_names=cellfun(@(x)x(1:strfind(x, '.')-1),  instruc.roi_files, 'UniformOutput',0); 
    end
    text.roi_names=instruc.roi_names; text.roi_names=cellfun(@(x)[x log.betasuffix], text.roi_names, 'UniformOutput', 0);
    request.FLcons=[]; if isempty(request.FLcons);  request.FLcons=instruc.FLcons; end
    
    

    
    % Confirm settings for printout
    disp('#### Printing data to excel file ########################################################')
    disp('(1) Short names for ROIs: '); disp([text.roi_names  repmat({'       ' }, length(instruc.roi_files), 1) instruc.roi_files]); disp(' ');
    disp('(2) FL contrast names:     ');disp(request.FLcons); disp(' '); input('Continue with these names?   '); disp(' '); disp(' -------------------------------' )
    
    % ########### *** EDIT HERE *** ###########
    
    t_betas= [[{'Subject'}; log.subjects] cell(log.n_subjs+1, length(instruc.roi_files)*length(request.FLcons)); ]; k=1;
    for r=1:size(instruc.roi_files,1)   % Print subject roi x contrast beta values 
        for n=1: length(request.FLcons)
            c=find(strcmp(instruc.FLcons, request.FLcons{n})); % speccific contrasts
            if isempty(c)==1; error(['Invalid first-level name requested!   ('  request.FLcons{n} ')']); end
            %
            t_betas{1,k+1}=[text.roi_names{r} '-' instruc.FLcons{c}]; % Title
%             disp(t_betas{1,k+1})
            for s=1:log.n_subjs
                t_betas{s+1, k+1}=d_betas{r+1, c+1}(s);
            end
            k=k+1;
        end
    end
    for o1=1:1  % Difference scores 
    
    % (if requested) Print difference scores for each ROI: cF_Reject-cF_Explore,  cF_Reject-ct_Bomb
    if request.calculateandprint_choicedifferencescores
        disp('Generating & printing difference scores ###########')
        
        % Instructions - specific difference-scores to compute, in order of request
        inclusterprefix=[]; 
        diffcontrasts={   {'cF_Accept'; 'ct_NoBomb'};       % Cross-task comparisons
                                  {'cF_Reject'; 'ct_Bomb'};
                                  {'cF_Explore'; 'ct_Explore'};
                                  {'cF_Reject'; 'cF_Explore'};      % Within task, compare choice
                                  {'cF_Reject'; 'cF_Accept'};
                                  {'cF_Accept'; 'cF_Explore'};
                                  {'ct_Bomb'; 'ct_Explore'};
                                  {'ct_Bomb'; 'ct_NoBomb'};
                                  {'ct_NoBomb'; 'ct_Explore'};
                              };
        
        inclusterprefix='in_'; 
        diffcontrasts={   {'cF_Reject'; 'ct_Bomb'};       % Cross-task comparisons
                                  {'cF_Explore'; 'ct_Explore'};
                                  {'cF_Reject'; 'cF_Explore'};      % Within task, compare choice
                                  {'ct_Bomb'; 'ct_Explore'};
                              };
        
        
        for r=1:size(instruc.roi_files,1)
            
            % Locate cell contrasts
%             c.cF_Accept=find(strcmp(instruc.FLcons, [inclusterprefix 'cF_Accept'])); if isempty(c.cF_Accept)==1; error('cannot find requested contrast, difference scores'); end
            c.cF_Reject=find(strcmp(instruc.FLcons, [inclusterprefix 'cF_Reject'])); if isempty(c.cF_Reject)==1; error('cannot find requested contrast, difference scores'); end
            c.cF_Explore=find(strcmp(instruc.FLcons, [inclusterprefix 'cF_Explore']));  if isempty(c.cF_Explore)==1; error('cannot find requested contrast, difference scores'); end
            %
%             c.ct_NoBomb=find(strcmp(instruc.FLcons, [inclusterprefix 'ct_NoBomb'])); if isempty(c.ct_NoBomb)==1; error('cannot find requested contrast, difference scores'); end
            c.ct_Bomb=find(strcmp(instruc.FLcons, [inclusterprefix 'ct_Bomb'])); if isempty(c.ct_Bomb)==1; error('cannot find requested contrast, difference scores'); end
            c.ct_Explore=find(strcmp(instruc.FLcons, [inclusterprefix 'ct_Explore']));  if isempty(c.ct_Explore)==1; error('cannot find requested contrast, difference scores'); end
            
            
            
            % Apply contrasts (instructed above in diffcontrasts)
            for jj=1:size(diffcontrasts,1)
                cons=diffcontrasts{jj};
                t_betas{1,k+1}=[text.roi_names{r} '-' inclusterprefix cons{1} '_minus_' inclusterprefix cons{2}]; disp(t_betas{1,k+1})
                eval(['a{1}=c.' cons{1} ';']); eval(['a{2}=c.' cons{2} ';']);
                for s=1:log.n_subjs
                    t_betas{s+1, k+1}=d_betas{r+1, a{1}+1}(s)-d_betas{r+1, a{2}+1}(s);
                end
                k=k+1;
            end
            
            
            
        end
    end
    end
%     t_betas=[t_betas [{'BORDER_BEHAVIOUR'}; cell(size(t_betas,1)-1,1)]];
%     xlswrite(['('  date ') Extracted betas' log.betasuffix '.xlsx'], t_betas);  movefile(['('  date ') Extracted betas' log.betasuffix '.xlsx'],  [where.model_resrois  '('  date ') Extracted betas' log.betasuffix '.xlsx'])   % Write file without the behavioural profile

    % Request behavioural profiles to append
    request.behfile=[where.beh filesep 'All behavioural data 15-08-13.xlsx'];
    [n t r]=xlsread(request.behfile);
    n=num2cell(n); n(isnan(cell2mat(n(:))))=repmat({[]}, sum(isnan(cell2mat(n(:)))),1); n=reshape(n, size(t)-1);
    [t_betas wb.headers wb.subs wb.n_subs excl] =f_aligntables([t_betas [{'BEH_'};  cell(size(t_betas,1)-1,1)]] , [t(:,1) [t(1,2:end); n]]);
    if isempty(excl.table1) + isempty(excl.table2)~=2
        disp('Subjects mismatch while appending beh');
        disp('         (subjects that are not present in both tables are excluded entirely)'); disp(excl); input('Continue?   ');
    end
    [st{1} msg{2}]=xlswrite(['('  date ') Extracted betas' log.betasuffix '.xlsx'], t_betas);  [st{2} msg{2}]=movefile(['('  date ') Extracted betas' log.betasuffix '.xlsx'],  [where.model_resrois  '('  date ') Extracted betas' log.betasuffix '.xlsx']);   % Write file without the behavioural profile
    disp(['        Excel file written:                      ' num2str(st{1})])
    disp(['        File moved to correct location:  ' num2str(st{2})])

end

%%

error('Done :)')

b_state=cell2mat(t_betas(2:end, find(strcmp(t_betas(1,:), 'st.State'))));
b_trait=cell2mat(t_betas(2:end, find(strcmp(t_betas(1,:), 'st.Trait'))));
% b_corres=cell2mat(length(instruc.roi_names)


for r=1:length(instruc.roi_names)
    [wr.r wr.p]=corr(b_state,   cell2mat(t_betas(2:end, find(strcmp(t_betas(1,:), [instruc.roi_names{r} '-cF_Reject']))))- cell2mat(t_betas(2:end, find(strcmp(t_betas(1,:), [instruc.roi_names{r} '-cF_Explore'])))));
    [wr.r wr.p]=corr(b_trait,   cell2mat(t_betas(2:end, find(strcmp(t_betas(1,:), [instruc.roi_names{r} '-cF_Reject']))))- cell2mat(t_betas(2:end, find(strcmp(t_betas(1,:), [instruc.roi_names{r} '-cF_Explore'])))));
    
    
end



openvar t_betas



t_betas(1,:)











