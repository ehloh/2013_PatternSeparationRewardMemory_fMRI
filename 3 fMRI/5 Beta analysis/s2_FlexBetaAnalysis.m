% Flexible beta analysis script 
clear all; close all hidden; clc, request.rois={}; request.whichbat={}; request.roi_rename=[]; 
% cd('F:\1 [Context-Memory] fMRI Data\2 Second level results s4FullCardiacWithDerivAnts')
cd('/Users/EleanorL/Desktop/1 CONTEXT fmri data/2 Second level results s4FullCardiacWithDerivAnts')

% request.where_FLrois=['m_ci3_ContextItemNomotor_Hit_Landmarks5' filesep 'cm_m1_2x2' filesep 'ROI']; 
% request.where_FLrois=['m_ci2_ContextonlyItem_Hit_Landmarks5' filesep 'cm_m1_2x2' filesep 'ROI']; 
% request.where_FLrois='m_i4_ItemContextevent_Hit_Landmarks5\iom_m1_2x2x2\ROI'; 
% request.where_FLrois='m_i5_ItemContextonly_Hit_Landmarks5\iom_m1_2x2x2\ROI'; 
% request.where_FLrois='m_c4_ContextallItempresonly_Hit_Landmarks5\cm_m1_2x2\ROI';
% request.where_FLrois=['m_i3_ItemContextpresonly_Hit_Landmarks5' filesep 'iom_m1_2x2x2' filesep 'ROI']; 
% request.where_FLrois='m_i6_ItemMempmodContextpresonly_Hit_Landmarks5\im_m1_2x2\ROI'; 
% request.where_FLrois='m_c7_ContexteventItempresonly_Hit_Landmarks5\co_m1_2x2\ROI'; 
% request.where_FLrois='m_ci3 PPIs\PPI comparison\coFam_sph_SNVTA_R1_cmSRvSN_psy_SR-SN\ROI';
request.where_FLrois='m_ci3 PPIs\PPI cell\coFam_sph_SNVTA_R2_cmSRvSN_psy_Cell';

% WHICH rois? #################
request.where_rois=request.where_FLrois; 
request.where_rois=[request.where_FLrois filesep 'ci3 rois'];     %  request.rois={'HPC_aL_svtp'; 'SNVTA_R_svtp'}; 
% request.where_rois=[request.where_FLrois filesep 'ci3 rois\peaks'];   
% request.where_rois=[request.where_FLrois filesep 'ci3 rois\Second string'];   
% request.where_rois=[request.where_FLrois filesep 'Anat Subfield'];   
% request.where_rois=[request.where_FLrois filesep 'Anat'];   
% request.where_rois=[request.where_FLrois filesep 'i6 rois'];   
% request.where_rois=[request.where_FLrois ];  request.whichbat=[]; 

% HPC rois only !! 
% request.whichbat=[];
% request.rois={'HPC_aL_svtp';'SNVTA_R_svtp'}; 
% request.rois={'rCA1_aL' 'rCA1_aR' 'rCA3_aL' 'rCA3_aR'};     % Subfield binary/re
request.rois={
% 'HPC_aL'; 'HPC_aR'; 'HPC_pL'; 'HPC_pR'
% 'SNVTA_L'; 'SNVTA_R' 
}; 

% Con/Beta Settings   #################
log.FLtype=[]; % Empty to auto select
% log.FLtype='cm';
% log.FLtype='co';
% log.FLtype='iom'; log.iom_analysis='SimValMemhit'; % Full S x V x Memhit analysis
% log.FLtype='iom';log.iom_analysis= 'SimVal'; % Do all analysis on mfx, Sim x Val 
% log.FLtype='iohit'; log.iom_analysis= 'SimVal'; % Do all analysis on mfx, Sim x Val 
% log.FLtype='io';log.io_analysis= 'xMem_Hit'; % Sim x Val analysis on item-mem pmods  
% log.FLtype='io';log.io_analysis= []; % Sim x Val analysis on item events 
log.FLtype='co_ppi_comp';
% log.FLtype='co_ppi_cell';



request.kendalls_correlation=0; 
% request.kendalls_correlation=1; 
request.meancentre=0;

%% Set up data files 

for o1=1:1 % ROI BATTERY requests  
    
    % [c3 battery] ##################
    request.c3batt_rois={
        'HPC_L_satc'; 'HPC_R_satc';
        };
    request.c3batt_roi_rename={
        };
     
end
if isempty(request.whichbat )==0, eval(['request.rois=request.' request.whichbat '_rois;  request.roi_rename=request.' request.whichbat '_roi_rename;']) , end
% request.roi_rename=[];  % Empty to omit % request.rois={};  % Empty to use all ROIs (without re-ordering)

% ################################################
% 
for o1=1:1   % Manual settings
    log.specificsubjects=[];
    % Get file name
    request.filename=[];
    if isempty(request.filename)
        %         f=spm_select('List',request.where_rois, 'Extracted');
        %         if size(f,1)==0; error('No beta excel files!'); end
        %         if size(f,1)~=1; error('No. of beta excel files ~=1!'); end
        %         request.filename= f(1:strfind(f, '.xls')-1);
        
        cd( request.where_rois)
        f= dir('*Extracted*');
        if size(f,1)==0; error('No beta excel files!'); end
        if size(f,1)~=1; error('No. of beta excel files ~=1!'); end
        request.filename= f.name(1:strfind(f.name, '.xls')-1) ;
    end
    
    % Other non-changing setup/details
    log.roi_con_connector='-';
    
end
for o1=1:1 % Load file + auto-read details about betas (rename beh, beta, etc)
    w=pwd;  if strcmp(w(1), '/'); where.where='/Users/EleanorL/Dropbox/SCRIPPS/1 ContextMem fMRI'; else  where.where='D:\Dropbox\SCRIPPS\1 ContextMem fMRI' ; end; addpath(where.where)
    
    % Load file + select ROIs/subjects
%     [n, t r]=xlsread([request.where_rois filesep request.filename '.xlsx']);
%     cd(request.where_rois)
    [n, t r]=xlsread([request.filename '.xlsx']);
    d_betas=[t(:,1) [strtrim( t(1, 2:end) ); num2cell(n)]]; log.orig_d_betas = d_betas; % Original table, save!
    if isempty(request.rois)==0
        w.roicon_names=d_betas(1,2:find(cellfun(@(x)~isempty(x),  strfind(d_betas(1,1:end), 'BEH')), 1, 'first')-1); w.okconcols=[];  % ROI selection
        for r=1:length(request.rois)
            if sum(cellfun(@(x)~isempty(x),   strfind(w.roicon_names, request.rois{r}))) ==0;   disp('Column titles:' );  disp(d_betas(1,:)');  error(['Could not find requested ROI. Check name: ' request.rois{r}]); end
            
            w.okconcols=[w.okconcols   find(cellfun(@(x)~isempty(x),   strfind(w.roicon_names, [request.rois{r} log.roi_con_connector]))) ];
        end
        d_betas=[d_betas(:, [1  w.okconcols+1]) d_betas(:,  find(cellfun(@(x)~isempty(x),  strfind(d_betas(1,:), 'BEH')), 1, 'first') : end)];
    end
    [log.subjects, log.n_subjs d_betas] = f_selectsubjects(d_betas, log.specificsubjects,  [d_betas(:,1) [{'All'}; num2cell(ones(size(d_betas,1)-1,1))]], 'All');  % Subject selection
    
    % Rename ROIs (proper name)
    if isempty(request.roi_rename)==0
        disp('Renaming of ROIs requested. Check names:'); openvar ans; ans=[request.rois request.roi_rename]; disp([request.rois request.roi_rename]); input('Continue?  ');
        w.roicon_names=d_betas(1,:);
        for r=1: length(request.rois)
            disp(['Assumed no. contrasts per roi: '  num2str(sum(cellfun(@(x)~isempty(x),  strfind(w.roicon_names,  request.rois{r}))))])
            w.cons2rename=find(cellfun(@(x)~isempty(x),  strfind(w.roicon_names,  request.rois{r})));
            for i=1:length(w.cons2rename)
                d_betas{1, w.cons2rename(i)}=[request.roi_rename{r} d_betas{1, w.cons2rename(i)}(length(request.rois{r})+1:end)];
            end
        end
        ans=[w.roicon_names' d_betas(1,:)']; input('Check line up of rename ROIs (in var window). Proceed?  ');
        request.rois=request.roi_rename;
    end
    
    % Manual renaming of behavioural variables
    %     d_betas{1, find(cellfun(@(x) ~isempty(x), strfind(d_betas(1,:),'st.State')))} = 'State anxiety scores';
    
    % Read out beta details from table
    log.data_names= d_betas(1,2:end)';
    log.roicons=log.data_names(1:find(1-cellfun(@(x)isempty(x), strfind(log.data_names, 'BEH_')),1)-1);
    
    % Auto-read FL con names and ROI names
    %     log.roicons=d_betas(1,:)'; log.roicons= log.roicons(2: find(cellfun(@(x)[1- isempty(x)],  strfind(log.roicons,'BEH_')),1,'first')-1);
    wc.rc=cellfun(@(x)fliplr(x), log.roicons, 'UniformOutput',0);
    wc.firstconcharnum=strfind(wc.rc{1}, log.roi_con_connector);
    % [Identify based on Cons (if roi names are non-unique)] ------------------------
    %     wc.firstconname=  wc.rc{1}(1:   wc.firstconcharnum(end)  );    
    %     wc.firstconnums_perroi= find(cellfun(@(x)[1- isempty(x)],  strfind(wc.rc,wc.firstconname)));
    %     wc.firstcons_perroi= log.roicons(wc.firstconnums_perroi);
    %     wc.charsin1stroiconname=cell2mat(strfind(wc.firstcons_perroi, fliplr(wc.firstconname)))-1;
    %     for i=1:size(wc.firstcons_perroi,1)
    %         log.rois_available{i,1}=  wc.firstcons_perroi{i}( 1: wc.charsin1stroiconname(i));
    %     end
    %     disp('Available ROIs:'), disp(log.rois_available)
    %     log.n_cons = wc.firstconnums_perroi(2) - wc.firstconnums_perroi(1);
    %     wc.roi1cons= log.roicons(1:log.n_cons );
    %     wc.roi1_name= wc.roi1cons{1}(1:strfind(wc.roi1cons{1}, log.roi_con_connector)-1 );
    %     log.cons=cellfun(@(x)x(length(wc.roi1_name) +2:end),  wc.roi1cons, 'UniformOutput',0);
    %     log.n_rois    = sum(1-cellfun(@(x)isempty(x), strfind(log.roicons, log.cons{1})));  % If there are repeats here, the no. of ROIs will be wrong!
    %     log.n_cons=length(log.cons);
    %     log.rois= cellfun(@(x)x(  1:  length(x)-length(log.cons{1})-length(log.roi_con_connector)   ), log.roicons(find(1-cellfun(@(x)isempty(x), strfind(log.roicons, log.cons{1})))), 'UniformOutput',0) ;
    
    % [Identify based on ROIs (if con names are non-unique] ------------------------
    wc.firstroiname =  fliplr(wc.rc{1}(wc.firstconcharnum(end)+1:end  ));    
    log.n_cons= length( find(cellfun(@(x)[1- isempty(x)],  strfind(log.roicons, wc.firstroiname))) ) ; 
    log.cons = cellfun(@(x)x(length(wc.firstroiname) + length(log.roi_con_connector)+1:end), log.roicons(1:log.n_cons), 'UniformOutput',0);
    log.rois = cellfun(@(x)x(1: length(x)- length(log.cons{1})-length(log.roi_con_connector)), log.roicons( 1:log.n_cons: length(log.roicons)), 'UniformOutput',0); 
    log.n_rois= length( log.rois ); 

    % Add behavioural scores 
    k=size(d_betas,2)+1; 
    d_betas{1,  k}  ='eRT_cell.Sim_valfx';  d_betas(2:end, k) =  num2cell( cell2mat(d_betas(2:end, strcmp(d_betas(1,:), 'eRT_cell.Sim_cN') )) - cell2mat(d_betas(2:end, strcmp(d_betas(1,:), 'eRT_cell.Sim_cR') )) )  ;  k=k+1; 
    d_betas{1,  k}  ='eRT_cell.Dis_valfx';  d_betas(2:end, k) =  num2cell( cell2mat(d_betas(2:end, strcmp(d_betas(1,:), 'eRT_cell.Dis_cN') )) - cell2mat(d_betas(2:end, strcmp(d_betas(1,:), 'eRT_cell.Dis_cR') )) )  ;  k=k+1; 

    % Read out behaviour details from table
    log.behfams_start=find(1-cellfun(@(x)isempty(x), strfind(log.data_names, 'BEH_')))+1;  %
    log.behnames=d_betas(1,log.behfams_start(1):end)';
end
for o=1:1 % Identify FL con types (determines figures, stats, correlations to perform)
    log.con_diffs={};
    if isempty(log.FLtype)
         error('FL type not specified yet!');
    end

    % FL Constrast differences scores (requested)
    if sum(strcmp(log.FLtype, {'co'}))==1
        log.con_diffs={
            'Sim_vfx'    {'SimR' '-' 'SimN'}   % All simple effects 
            'Dis_vfx'    {'DisR' '-' 'DisN'}             
            'Rew_sfx'    {'SimR' '-' 'DisR'}  
            'Neu_sfx'    {'SimN' '-' 'DisN'}  
            };
    elseif sum(strcmp(log.FLtype, {'iohit'}))==1 
        log.con_diffs={
            'Sim_vfx'    {'SR_Hit' '-' 'SN_Hit'}   % All simple effects 
            'Dis_vfx'    {'DR_Hit' '-' 'DN_Hit'}  
            'Rew_sfx'    {'SR_Hit' '-' 'DR_Hit'}  
            'Neu_sfx'    {'SN_Hit' '-' 'DN_Hit'}  
            };
    elseif sum(strcmp(log.FLtype, {'io'}))==1
        log.con_diffs={
            'Sim_vfx'    {['SR_Item' log.io_analysis] '-' ['SN_Item' log.io_analysis]}   % All simple effects 
            'Dis_vfx'    {['DR_Item' log.io_analysis] '-' ['DN_Item' log.io_analysis]}             
            'Rew_sfx'    {['SR_Item' log.io_analysis] '-' ['DR_Item' log.io_analysis]}  
            'Neu_sfx'    {['SN_Item' log.io_analysis] '-' ['DN_Item' log.io_analysis]}  
            };
        if strcmp(log.io_analysis, 'xMem_Hit')==1,  disp('ItemMem pmods chosen');  else disp('Item events (not pmods) chosen');  end 
    elseif sum(strcmp(log.FLtype, {'cm'}))==1
        log.con_diffs={
            'Sim_vfx'    {'SimRxCMem_Hit' '-' 'SimNxCMem_Hit'}   % All simple effects 
            'Dis_vfx'    {'DisRxCMem_Hit' '-' 'DisNxCMem_Hit'}             
            'Rew_sfx'    {'SimRxCMem_Hit' '-' 'DisRxCMem_Hit'}  
            'Neu_sfx'    {'SimNxCMem_Hit' '-' 'DisNxCMem_Hit'}  
            };
        
    elseif sum(strcmp(log.FLtype, {'iom'}))==1
        log.con_diffs={
            'SR_memfx'    {'SR_Hit' '-' 'SR_Miss'}   % Memfx 
            'SN_memfx'    {'SN_Hit' '-' 'SN_Miss'}   
            'DR_memfx'    {'DR_Hit' '-' 'DR_Miss'}   
            'DN_memfx'    {'DN_Hit' '-' 'DN_Miss'}   
            'SimHit_valfx'    {'SR_Hit' '-' 'SN_Hit'}   % Within memhit, compare cells 
            'DisHit_valfx'    {'DR_Hit' '-' 'DN_Hit'}   
            'RewHit_Simfx'    {'SR_Hit' '-' 'DR_Hit'}   
            'NeuHit_Simfx'    {'SN_Hit' '-' 'DN_Hit'}   
            };
    elseif sum(strcmp(log.FLtype, {'co_ppi_cell'}))==1
        log.con_diffs={
            'Sim_vfx'    {'SimR' '-' 'SimN'}   % All simple effects 
            'Dis_vfx'    {'DisR' '-' 'DisN'}             
            'Rew_sfx'    {'SimR' '-' 'DisR'}  
            'Neu_sfx'    {'SimN' '-' 'DisN'}  
            };
    else log.con_diffs=cell(0,2); disp('No simple effects')
    end
    
end
for o1=1:1  % Preprocessing of betas (mean-centre, new beta scores)
    d_newconbetas=cell(log.n_subjs+1,0); log.newcons={};
    for r=1:log.n_rois
        %         % Reverse specific betas!!!!
        %         wr=cell2mat(d_betas(2:end,  find(cellfun(@(x)~isempty(x),  strfind(d_betas(1,:),  log.rois{r})))));
        %         wr=-1 *wr; d_betas(2:end,  find(cellfun(@(x)~isempty(x),  strfind(d_betas(1,:),  log.rois{r}))))=num2cell(wr);
        %         disp(['Reversed-valence betas for ' log.rois{r}])
        
        if request.meancentre
            switch log.FLtype
                case 'co' 
                    wr= cell2mat(  d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'SimR'])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'DisN'])) )    );
                    wr= wr- repmat(mean(wr,2), 1, size(wr,2)); wr=num2cell(wr);
%                     d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{1}])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{end}])) ) = wr;                  
                    d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'SimR'])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'DisN'])) )     = wr;
                    disp('Betas mean centred across co contrasts!')
                case 'cm'
                    if r==1; input('Requested mean centring for cm. ARE YOU SURE? you shoulnt'); end
                    wr= cell2mat(  d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'SimRxCMem_Hit'])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'DisNxCMem_Hit'])) )    );
                    wr= wr- repmat(mean(wr,2), 1, size(wr,2)); wr=num2cell(wr);
%                     d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{1}])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{end}])) ) = wr;                  
                    d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'SimRxCMem_Hit'])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'DisNxCMem_Hit'])) ) =wr;
                    if r==1; disp('Betas mean centred across co contrasts !'), end
                case 'iom'
                    wr= cell2mat(  d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'SR_Hit'])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'DN_Miss'])) )    );
                    wr= wr- repmat(mean(wr,2), 1, size(wr,2)); wr=num2cell(wr);
%                     d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{1}])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{end}])) ) = wr;
                    d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'SR_Hit'])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'DN_Miss'])) )=wr;     
                    if r==1; disp('Betas mean centred across iom contrasts !'), end
                case 'iohit'
                    wr= cell2mat(  d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'SR_Hit'])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'DN_Miss'])) )    );  
                    wr= wr- repmat(mean(wr,2), 1, size(wr,2)); wr=num2cell(wr);
%                     d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{1}])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{end}])) ) = wr;
                    d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'SR_Hit'])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'DN_Miss'])) )=wr;     
                    if r==1; disp('Betas mean centred across io hit & miss contrasts !'), end
                case 'co_ppi_cell'
                    wr= cell2mat(  d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'SimR'])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'DisN'])) )    );
                    wr= wr- repmat(mean(wr,2), 1, size(wr,2)); wr=num2cell(wr);
%                     d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{1}])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{end}])) ) = wr;                  
                    d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'SimR'])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector 'DisN'])) )     = wr;
                    disp('Betas mean centred across all cell ppi contrasts! Are you sure?')
                otherwise, disp('Mean centring across ALL contrasts');
                    wr= cell2mat(  d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{1}])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{end}])) )    );
                    wr= wr- repmat(mean(wr,2), 1, size(wr,2)); wr=num2cell(wr);
                    d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{1}])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{end}])) ) = wr;
                    if r==1; disp('Betas mean centred across ALL contrasts! Probably wrong.'), end
            end
        end
        
        
        % Calculate difference scores
        for d=1:size(log.con_diffs,1)
            if sum(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.con_diffs{d,2}{1}])) + sum(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.con_diffs{d,2}{end}])) ~=2; error(['Error calculating difference scores. Could not find columns headed ' log.rois{r} log.roi_con_connector log.con_diffs{d,2}{1}  '  or  ' log.rois{r} log.roi_con_connector log.con_diffs{d,2}{3} ]);  end
            
            % Find data to manipulate
            wd{1}=cell2mat(  d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.con_diffs{d,2}{1}]))) ) ;
            wd{2}= cell2mat(    d_betas(2:end,  find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.con_diffs{d,2}{3}])))  );
            
            % Compute difference (assume 1 operation only)
            eval(['wd{3} = wd{1} '   log.con_diffs{d,2}{2}  'wd{2};'])
            
            % Add new column to existing data
            d_betas=[d_betas [cellstr([log.rois{r} log.roi_con_connector log.con_diffs{d}] ); num2cell(wd{3})]];
            
            wd=[];
        end
        
        
        % #######################################################
        % Manually specified new scores! ###############################
        % #######################################################
        if sum(strcmp(log.FLtype, {'cm';'co'; 'co_ppi_cell'})) 
            if r==1, disp('Calculating new betas for SimVal'), end
            in.newconname='SR_vOth'; in.cons={'SimR';'SimN'; 'DisR'; 'DisN'};            
            for c=1:length(in.cons)
                in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
            end
            d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname]; if r==1; log.newcons=[log.newcons; in.newconname]; end
            d_newconbetas(2:end, end)=num2cell(   in.b{1} -   mean([in.b{2} in.b{3}  in.b{4}],2)   );  % Manually compute score
        elseif sum(strcmp(log.FLtype, {'io'})) 
            if r==1, disp('Calculating new betas for SimVal'), end            
            in.newconname='SR_vOth'; in.cons={['SR_Item' log.io_analysis];['SN_Item' log.io_analysis]; ['DR_Item' log.io_analysis]; ['DN_Item' log.io_analysis]};
            for c=1:length(in.cons)
                in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
            end
            d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname]; if r==1; log.newcons=[log.newcons; in.newconname]; end
            d_newconbetas(2:end, end)=num2cell(   in.b{1} -   mean([in.b{2} in.b{3}  in.b{4}],2)   );  % Manually compute score
        elseif sum(strcmp(log.FLtype, {'iom'})) 
            if r==1, disp('Calculating new betas for Sim x Val x Memhit'), end
            in.newconname='mfx_SRmOth'; in.cons={'SR_memfx'; 'SN_memfx'; 'DR_memfx'; 'DN_memfx'};
            for c=1:length(in.cons)
                in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
            end
            d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname]; if r==1; log.newcons=[log.newcons; in.newconname]; end
            d_newconbetas(2:end, end)=num2cell(   3*in.b{1} -   mean([in.b{2} in.b{3}  in.b{4}],2)   );  % Manually compute score
            in.newconname='mfx_SRmSN'; in.cons={'SR_memfx'; 'SN_memfx'; 'DR_memfx'; 'DN_memfx'};
            for c=1:length(in.cons)
                in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
            end
            d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname]; if r==1; log.newcons=[log.newcons; in.newconname]; end
            d_newconbetas(2:end, end)=num2cell(   in.b{1} -    in.b{2} );  % Manually compute score
            in.newconname='mfx_SRmDR'; in.cons={'SR_memfx'; 'SN_memfx'; 'DR_memfx'; 'DN_memfx'};
            for c=1:length(in.cons)
                in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
            end
            d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname]; if r==1; log.newcons=[log.newcons; in.newconname]; end
            d_newconbetas(2:end, end)=num2cell(   in.b{1} -    in.b{3} );  % Manually compute score

            % Independent of memory 
            in.newconname='SR'; in.cons={'SR_Hit'; 'SR_Miss'; 'SN_Hit'; 'SN_Miss'; 'DR_Hit'; 'DR_Miss'; 'DN_Hit'; 'DN_Miss'};
            for c=1:length(in.cons)
                in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
            end
            d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname]; if r==1; log.newcons=[log.newcons; in.newconname]; end
            d_newconbetas(2:end, end)=num2cell(  mean([ in.b{1}  in.b{2}],2) );  % Manually compute score
            in.newconname='SN'; in.cons={'SR_Hit'; 'SR_Miss'; 'SN_Hit'; 'SN_Miss'; 'DR_Hit'; 'DR_Miss'; 'DN_Hit'; 'DN_Miss'};
            for c=1:length(in.cons)
                in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
            end
            d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname]; if r==1; log.newcons=[log.newcons; in.newconname]; end
            d_newconbetas(2:end, end)=num2cell(  mean([ in.b{3}  in.b{4}],2) );  % Manually compute score
            in.newconname='DR'; in.cons={'SR_Hit'; 'SR_Miss'; 'SN_Hit'; 'SN_Miss'; 'DR_Hit'; 'DR_Miss'; 'DN_Hit'; 'DN_Miss'};
            for c=1:length(in.cons)
                in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
            end
            d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname]; if r==1; log.newcons=[log.newcons; in.newconname]; end
            d_newconbetas(2:end, end)=num2cell(  mean([ in.b{5}  in.b{6}],2) );  % Manually compute score
            in.newconname='DN'; in.cons={'SR_Hit'; 'SR_Miss'; 'SN_Hit'; 'SN_Miss'; 'DR_Hit'; 'DR_Miss'; 'DN_Hit'; 'DN_Miss'};
            for c=1:length(in.cons)
                in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
            end
            d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname]; if r==1; log.newcons=[log.newcons; in.newconname]; end
            d_newconbetas(2:end, end)=num2cell(  mean([ in.b{7}  in.b{8}],2) );  % Manually compute score
            in.newconname='SRmSN'; in.cons={'SR_Hit'; 'SR_Miss'; 'SN_Hit'; 'SN_Miss'; 'DR_Hit'; 'DR_Miss'; 'DN_Hit'; 'DN_Miss'};
            for c=1:length(in.cons)
                in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
            end
            d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname]; if r==1; log.newcons=[log.newcons; in.newconname]; end
            d_newconbetas(2:end, end)=num2cell(  mean([ in.b{1}  in.b{2}],2)  -  mean([ in.b{3}  in.b{4}],2)  );  % Manually compute score
            in.newconname='SRmDR'; in.cons={'SR_Hit'; 'SR_Miss'; 'SN_Hit'; 'SN_Miss'; 'DR_Hit'; 'DR_Miss'; 'DN_Hit'; 'DN_Miss'};
            for c=1:length(in.cons)
                in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
            end
            d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname]; if r==1; log.newcons=[log.newcons; in.newconname]; end
            d_newconbetas(2:end, end)=num2cell(  mean([ in.b{1}  in.b{2}],2)  -  mean([ in.b{5}  in.b{6}],2)  );  % Manually compute score
        elseif sum(strcmp(log.FLtype, {'iohit'}))==1
            if r==1, disp('Calculating new betas for Sim x Val x Memhit'), end
            in.newconname='hit_SRmOth'; in.cons={'SR_Hit'; 'SN_Hit'; 'DR_Hit'; 'DN_Hit'};
            for c=1:length(in.cons)
                in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
            end
            d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname]; if r==1; log.newcons=[log.newcons; in.newconname]; end
            d_newconbetas(2:end, end)=num2cell(   3*in.b{1} -   mean([in.b{2} in.b{3}  in.b{4}],2)   );  % Manually compute score
            in.newconname='hit_SRmSN'; in.cons={'SR_Hit'; 'SN_Hit'; 'DR_Hit'; 'DN_Hit';};
            for c=1:length(in.cons)
                in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
            end
            d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname]; if r==1; log.newcons=[log.newcons; in.newconname]; end
            d_newconbetas(2:end, end)=num2cell(   in.b{1} -    in.b{2} );  % Manually compute score
            in.newconname='hit_SRmDR'; in.cons={'SR_Hit'; 'SN_Hit'; 'DR_Hit'; 'DN_Hit';};
            for c=1:length(in.cons)
                in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
            end
            d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname]; if r==1; log.newcons=[log.newcons; in.newconname]; end
            d_newconbetas(2:end, end)=num2cell(   in.b{1} -    in.b{3} );  % Manually compute score
        else disp('Simple fx not set up yet!')
        end
                
        % #######################################################
        % #######################################################
        
        wr=[];
    end
        
end
d_betas=[d_betas d_newconbetas];  % Append new betas
log.cons=[log.cons; log.con_diffs(:,1); log.newcons];  log.n_cons=length(log.cons);

% % Subset of subjects ?
% whichsubs=[5 11 24 20 10 9 3 15 12 7 13 17 22]; % HIGH median split by dprime  Dis valfx 
% whichsubs=[9 3 15 12 7 13 17 22]; % dprime DisValfx >0
% log.subjects= log.subjects(whichsubs); log.n_subjs=length(log.subjects);  d_betas = [d_betas(1,:);  d_betas(whichsubs+1, :)];input('Some subjects only! Proceed?') ;


% Z score beh/betas
% zthis='eRT_cell.Sim_cR'; 
% d_betas(2:end, strcmp(d_betas(1,:), zthis)) =  num2cell(   zscore(cell2mat(d_betas(2:end, strcmp(d_betas(1,:), zthis))))   );  disp(['z scored ' zthis ' !!!]' ])



%% Get all stats (r_stats)
%   Col 2: Mean, Std Error, One-sample tstat, df, One-sample p  (row=contrast)
%   Col 3: 2x2 factorial
%   Col 4: Simple-effects ttests

dostats=1;

% Get statistics (switch things on and off here!!!)
if dostats
    r_stats=[log.rois cell(log.n_rois, 5)]; k=1;  printout={}; p=1;
    for r=1:log.n_rois
        
        % Collect betas for all contrasts
        wr.betanames=cellfun(@(x)[log.rois{r}   log.roi_con_connector  x], log.cons, 'UniformOutput',0); wr.d=[];
        for i=1:length(wr.betanames)
            wr.d=[wr.d cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:),  wr.betanames{i}))))];
        end
        
        % (i) Mean + Std error + tstat + tstate df + tstat p
        r_stats{r, 2}=[ mean(wr.d)'     (std(wr.d)/sqrt(size(wr.d,1)))'    ];
        [wr.h wr.p wr.ci wr.stats]=ttest(wr.d);
        r_stats{r, 2}=[r_stats{r, 2} wr.stats.tstat' wr.stats.df' wr.p']; k=1;
    
        switch log.FLtype
            case  'co'
                wr.d=wr.d(:, [find(strcmp(log.cons, 'SimR')) find(strcmp(log.cons, 'SimN')) find(strcmp(log.cons, 'DisR')) find(strcmp(log.cons, 'DisN')) ]);
                [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,1:4) ,[2 2], {'Sim' 'Val'});
                r_stats{r, 3}=wr.anova.R;  % Row=Sim, Val, SxV; Col=F, df1, df2, p
            case  'io'
                wr.d=wr.d(:, [find(strcmp(log.cons, ['SR_Item' log.io_analysis]))  find(strcmp(log.cons, ['SN_Item' log.io_analysis]))  find(strcmp(log.cons, ['DR_Item' log.io_analysis]))  find(strcmp(log.cons, ['DN_Item' log.io_analysis])) ]);
                [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,1:4) ,[2 2], {'Sim' 'Val'});
                r_stats{r, 3}=wr.anova.R;  % Row=Sim, Val, SxV; Col=F, df1, df2, p
            case 'cm' 
                wr.d=wr.d(:, [find(strcmp(log.cons, 'SimRxCMem_Hit')) find(strcmp(log.cons, 'SimNxCMem_Hit')) find(strcmp(log.cons, 'DisRxCMem_Hit')) find(strcmp(log.cons, 'DisNxCMem_Hit')) ]);
                [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,1:4) ,[2 2], {'Sim' 'Val'});
                r_stats{r, 3}=wr.anova.R;  % Row=Sim, Val, SxV; Col=F, df1, df2, p
            case 'iohit' 
                    % Stats for full Sim x Val x Memhit
                    wr.d=wr.d(:, [ find(strcmp(log.cons, 'SR_Hit'))    find(strcmp(log.cons, 'SN_Hit'))   find(strcmp(log.cons, 'DR_Hit'))  find(strcmp(log.cons, 'DN_Hit'))   ]);
                    [wr.anova]=teg_repeated_measures_ANOVA(wr.d,[2 2], {'Sim' 'Val'});
                    r_stats{r, 3}=wr.anova.R;  % Row=Sim, Val, SxV; Col=F, df1, df2, p
            case 'iom' 
                if strcmp(log.iom_analysis, 'SimValMemhit')
                    % Stats for full Sim x Val x Memhit
                    wr.d=wr.d(:, [ find(strcmp(log.cons, 'SR_Hit'))  find(strcmp(log.cons, 'SR_Miss'))  find(strcmp(log.cons, 'SN_Hit'))  find(strcmp(log.cons, 'SN_Miss'))  find(strcmp(log.cons, 'DR_Hit')) find(strcmp(log.cons, 'DR_Miss')) find(strcmp(log.cons, 'DN_Hit'))  find(strcmp(log.cons, 'DN_Miss'))  ]);
                    [wr.anova]=teg_repeated_measures_ANOVA(wr.d,[2 2 2], {'Sim' 'Val' 'Memhit'});
                    r_stats{r, 3}=wr.anova.R;  % Row=Sim, Val, SxV; Col=F, df1, df2, p
                elseif strcmp(log.iom_analysis, 'SimVal')
                    % Stats on memfx, Sim x Val
                    wr.d=wr.d(:, [ find(strcmp(log.cons, 'SR_memfx'))  find(strcmp(log.cons, 'SN_memfx'))   find(strcmp(log.cons, 'DR_memfx')) find(strcmp(log.cons, 'DN_memfx'))   ]);
                    [wr.anova]=teg_repeated_measures_ANOVA(wr.d,[2 2], {'Sim' 'Val'});
                    r_stats{r, 3}=wr.anova.R;  % Row=Sim, Val, SxV; Col=F, df1, df2, p
                end
            case 'co_ppi_cell';
                wr.d=wr.d(:, [find(strcmp(log.cons, 'SimR')) find(strcmp(log.cons, 'SimN')) find(strcmp(log.cons, 'DisR')) find(strcmp(log.cons, 'DisN')) ]);
                [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,1:4) ,[2 2], {'Sim' 'Val'});
                r_stats{r, 3}=wr.anova.R;  % Row=Sim, Val, SxV; Col=F, df1, df2, p
            otherwise disp('stats not set up yet!');
        end
        
        %
        wr=[];
    end
    
    % Print to display ANOVA statistics
    if sum(strcmp(log.FLtype, {'cm';'co';'co_ppi_cell';'io';'iohit'}))
        for r=1:log.n_rois
            printout{p,1}=log.rois{r};
            
            disp([r_stats{r,1} ' ---'])
            c=1; disp( ['ME Sim: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='ME Sim';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            c=2; disp( ['ME Val : F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='ME Val';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            c=3; disp( ['SxV: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='Sim x Val';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            
            p=p+1;
            disp(' ')
        end
    elseif sum(strcmp(log.FLtype, {'iom'}))
        for r=1:log.n_rois
            printout{p,1}=log.rois{r}; disp([r_stats{r,1} ' ---'])
            if strcmp(log.iom_analysis, 'SimValMemhit')
                c=1; disp( ['ME Sim: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
                printout{p,2}='ME Sim';
                printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
                p=p+1;
                
                c=2; disp( ['ME Val : F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
                printout{p,2}='ME Val';
                printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
                p=p+1;
                
                c=3; disp( ['ME Memhit: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
                printout{p,2}='Memhit';
                printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
                p=p+1;
                
                c=4; disp( ['Sim x Val: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
                printout{p,2}='Sim x Val';
                printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
                p=p+1;
                
                c=5; disp( ['Sim x Memhit: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
                printout{p,2}='Sim x Memhit';
                printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
                p=p+1;
                
                c=6; disp( ['Val x Memhit: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
                printout{p,2}='Val x Memhit';
                printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
                p=p+1;
                
                c=7; disp( ['Sim x Val x Memhit: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
                printout{p,2}='Sim x Val x Memhit';
                printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
                p=p+1;
                
            elseif strcmp(log.iom_analysis, 'SimVal')
                c=1; disp( ['ME Sim: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
                printout{p,2}='ME Sim';
                printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
                p=p+1;
                
                c=2; disp( ['ME Val : F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
                printout{p,2}='ME Val';
                printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
                p=p+1;
                
                c=3; disp( ['SxV: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
                printout{p,2}='Sim x Val';
                printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
                p=p+1;
                
            end
            p=p+1;
            disp(' ')
        end
    end
end

% error('done stats')

%% Figures
%   Col 2: Mean, Std Error, One-sample tstat, df, One-sample p  (row=contrast)
%   Col 3: 2x2 factorial
%   Col 4: Simple-effects ttests

dofig=1;
for o1=1:1  % Figure settings
        fontsize=15;
        f.fontname='PT Sans Caption';  % pt serif (caption) ,san serif , pt sans,trebuchet ms
        % f.fontname='Cambria';
        % f.fontname='Arial';
end
if dofig
    close all hidden
    f.plotcols=4;  f.figwidth= 1800; f.figheight=400; f.fontsize=10; f.fontsize_title=25;
    f.subplot_VerHorz=[0.2 0.15]; f.fig_BotTop=[0.15 0.2]; f.fig_LeftRight=[0.1 0.1];
%     f.subplot_VerHorz=[0.1 0.1]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.05 0.05]; % Loads of rois
    figure('Name', 'ROI betas', 'NumberTitle', 'off', 'Position', [100 550 f.figwidth f.figheight], 'Color', 'w');  k=1;
    for r=1: log.n_rois
        set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname)
        
        switch log.FLtype
            case 'co'
                %             if sum(strcmp(log.FLtype, {'cm';'co'})) % [SimxVal ]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={'SimR';'SimN'; 'DisR';'DisN';};
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Sim' 'Dis' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
                %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                %                 ylabel('Parameter estimates (mean)', 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Reward', 'Neutral'), end
                % %                     xticklabel_rotate
                xlim([0.5 2.5])
            case 'cm'
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={  'SimRxCMem_Hit';'SimNxCMem_Hit'; 'DisRxCMem_Hit';'DisNxCMem_Hit';};
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Sim' 'Dis' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
                %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                %                 ylabel('Parameter estimates (mean)', 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Reward', 'Neutral'), end
                % %                     xticklabel_rotate
                xlim([0.5 2.5])
            case 'iom'
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;

                if strcmp(log.iom_analysis, 'SimValMemhit')
                    
                    % Plot full Sim x Val x Memhit
                    wr.whichcons_name={ 'SR_Hit'; 'SR_Miss'; 'SN_Hit'; 'SN_Miss'; 'DR_Hit'; 'DR_Miss'; 'DN_Hit'; 'DN_Miss'};
                    wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                    if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                    barwitherr(  r_stats{r,2}(wr.whichcons,2) ,    r_stats{r,2}(wr.whichcons,1) ); % Mean +/- 1 Std Err
                    set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:8, 'XTickLabel', cellfun(@(x)x(1:4),  wr.whichcons_name, 'UniformOutput',0),  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
                    ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                    title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
                    %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                    %                 if r==log.n_rois,  legend('Reward', 'Neutral'), end
                    %                 xlim([0.5 2.5])
                elseif strcmp(log.iom_analysis, 'SimVal')

                    % Plot memhit, Sim x Val
                    wr.whichcons_name={ 'SR_memfx'; 'SN_memfx';  'DR_memfx'; 'DN_memfx';  };
                    wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                    if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                    barwitherr( reshape( r_stats{r,2}(wr.whichcons,2),2,2) ,    reshape(r_stats{r,2}(wr.whichcons,1),2,2) ); % Mean +/- 1 Std Err
                    set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel',  {'Sim' 'Dis'} ,  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
                    ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                    title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
                    %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                    if r==log.n_rois,  legend('Reward', 'Neutral'), end; xlim([0.5 2.5])
                end
                
            case 'iohit'
                    % Plot memhit, Sim x Val
                    subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                    wr.whichcons_name={ 'SR_Hit'; 'SN_Hit';  'DR_Hit'; 'DN_Hit';  };
                    wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                    if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                    barwitherr( reshape( r_stats{r,2}(wr.whichcons,2),2,2) ,    reshape(r_stats{r,2}(wr.whichcons,1),2,2) ); % Mean +/- 1 Std Err
                    set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel',  {'Sim' 'Dis'} ,  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
                    ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                    title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
                    %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                    if r==log.n_rois,  legend('Reward', 'Neutral'), end; xlim([0.5 2.5])
                    
            case 'io'
                %             if sum(strcmp(log.FLtype, {'cm';'co'})) % [SimxVal ]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={['SR_Item' log.io_analysis]; ['SN_Item' log.io_analysis]; ['DR_Item' log.io_analysis]; ['DN_Item' log.io_analysis]};
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Sim' 'Dis' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
                %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                %                 ylabel('Parameter estimates (mean)', 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Reward', 'Neutral'), end
                % %                     xticklabel_rotate
                xlim([0.5 2.5])
            case 'co_ppi_cell'
                %             if sum(strcmp(log.FLtype, {'cm';'co'})) % [SimxVal ]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={'SimR';'SimN'; 'DisR';'DisN';};
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Sim' 'Dis' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
                %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                %                 ylabel('Parameter estimates (mean)', 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Reward', 'Neutral'), end
                % %                     xticklabel_rotate
                xlim([0.5 2.5])
                
                
            case 'co_ppi_comp'
                %             if sum(strcmp(log.FLtype, {'cm';'co'})) % [SimxVal ]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={'PPI Pos';};
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr(  r_stats{r,2}(wr.whichcons,2) ,     r_stats{r,2}(wr.whichcons,1)   ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:1, 'XTickLabel', {'Sim' 'Dis' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
                %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                %                 ylabel('Parameter estimates (mean)', 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
%                 if r==log.n_rois,  legend('Reward', 'Neutral'), end
                % %                     xticklabel_rotate
                xlim([0.5 2.5])
            otherwise error('Which plots? '); 
        end
        
        
        % Signicant?
        %         [h p ci stats]=ttest(cell2mat(d_betas(2:end, 1+[129 136])));
        %         stats
        %         p
        
        
        
        
%         ylim([-5 2])
%         xlim([0.5 2.5])
        
     
    
end
end


% wr.d=wr.d(:, [find(strcmp(log.cons, ['SR_Item' log.io_analysis]))  find(strcmp(log.cons, ['SN_Item' log.io_analysis]))  find(strcmp(log.cons, ['DR_Item' log.io_analysis]))  find(strcmp(log.cons, ['DN_Item' log.io_analysis])) ]);
%                 [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,1:4) ,[2 2], {'Sim' 'Val'});
%                 r_stats{r, 3}=wr.anova.R;  % Row=Sim, Val, SxV; Col=F, df1, df2, p
%             
% error('Done with stats and fig')



%% Battery correlations (specific)
%       r_batcorr{k}: row= beh, col= roi/con

% Certain ROIs?
request.rois=log.rois;
d_meanbetas=[d_betas(1,2:end)'  num2cell(mean(cell2mat(d_betas(2:end, 2:end)))') num2cell(std(cell2mat(d_betas(2:end, 2:end)))'./sqrt(log.n_subjs))  num2cell( ttest(cell2mat(d_betas(2:end, 2:end)))' )  ]; 
[h p ci st]=ttest(cell2mat(d_betas(2:end, 2:end))); p(p>0.1)=nan; d_meanbetas(:,4)=num2cell(p'); 
openvar d_meanbetas  % mean beta - simple effects
for i=1:size(d_meanbetas,1);
    
    if isnan(d_meanbetas{i,4})==1, 
        d_meanbetas{i,4}='nsf'; 
%         d_meanbetas{i,4}=' '; 
    else d_meanbetas{i,4}=['t('  num2str(st.df(i)) ')= '  num2str(st.tstat(i),3) ', p= '  num2str(p(i),3)];
    end
    
    
%     d_meanbetas{i,4}=['t('  num2str(st.df(i)) ')= '  num2str(st.tstat(i),3) ', p= '  num2str(p(i),3)];
    try d_meanbetas{i,5}= kstest(cell2mat( d_betas(2:end, i+1))   ); end 
end



% titles=d_betas(1,:)';  openvar titles
% error('done with simple fx')

% (1) Behaviour scores to correlate #########################################
request.cor_beh= {
    'eRT_cell.Sim_cR' ; 'eRT_cell.Sim_valfx'; 'dpr_cell.R_simfx';'dpr_cell.Sim_valfx'; 'dpr_Cell.Sim_R'
    'hr_cell.R_simfx';'hr_cell.Sim_valfx'; 'p.NS_all'; 'p.NS_1';  'p.NS_2'; 'p.NS_3'; 'p.NS_4';'p.RD_all';  'p.RD_1';   'p.RD_2';
    'p.RD_3';  'p.RD_4';  }; 


request.cor_beh= {  % Any Dis memfx ? 
%     'dpr_cell.Dis_valfx'; 'dpr_Cell.Dis_R';
%     'hr_cell.Dis_valfx';'hr_Cell.Dis_R';
 'eRT_cell.Dis_valfx';
}; 


for o=1:1 %  Adjust behaviour for correlations
    
 
    
%     k=size(d_betas,2)+1;
%     
%     d_betas{1,k}='per.ct_Reject';
%     d_betas(2:end, k)= num2cell( 2*cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.Reject'),1,'first'))) -  cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.cF_Reject'),1,'first')))     );
%     k=k+1;
%     
%     d_betas{1,k}='per.cF_Accept';
%     d_betas(2:end, k)= num2cell( 1- (cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.cF_Reject'),1,'first'))) +  cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.cF_Explore'),1,'first'))))     );
%     k=k+1;
%     d_betas{1,k}='per.ct_Accept';
%     d_betas(2:end, k)= num2cell( 1- (cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.ct_Reject'),1,'first'))) +  cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.ct_Explore'),1,'first'))))     );
%     k=k+1;
%     
%     d_betas{1,k}='per.Accept_ct-cF';
%     d_betas(2:end, k)= num2cell( cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.ct_Accept'),1,'first'))) -  cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.cF_Accept'),1,'first')))     );
%     k=k+1;
% 
%     request.cor_beh=[request.cor_beh;  {'per.Accept_ct-cF'; 'per.cF_Accept'}]; 
end
request.cor_rois=request.rois;
switch log.FLtype
    case 'cm';   request.cor_con={'SimRxCMem_Hit'; 'Sim_vfx'; 'Rew_sfx';'SR_vOth'}; 
    case 'co';   
        request.cor_con={'SimR'; 'Sim_vfx'; 'Rew_sfx';'SR_vOth'};
%         request.cor_con={'DisR'; 'Dis_vfx'};
    case 'io';   request.cor_con={['SR_Item' log.io_analysis]; 'Sim_vfx'; 'Rew_sfx';'SR_vOth'};
    case 'co_ppi_cell';   
%         request.cor_con={'SimR'; 'Sim_vfx'; 'Rew_sfx';'SR_vOth';};        
        request.cor_con={'SimR'};
    case 'co_ppi_comp';   request.cor_con={'PPI Pos';};
    case 'iom';    request.cor_con={'SR_Hit'; 'SR_memfx'; 'SimHit_valfx'; 'RewHit_Simfx'; 'mfx_SRmOth'; 'mfx_SRmSN'; 'mfx_SRmDR'; 'SR';'SN'; 'SRmSN'; 'SRmDR'};
    case 'iohit';    request.cor_con={'SR_Hit'; 'Sim_vfx'}; 
    otherwise, request.cor_con =  {};
end
rc=1; request.cor_roicons={}; % Compile roicons + execute battery
for r=1: length(request.cor_rois)    
    for c=1:length(request.cor_con)
        request.cor_roicons{rc,1}=[request.cor_rois{r}  log.roi_con_connector    request.cor_con{c}];  rc=rc+1;
    end
end
k=1; r_batcorr{k,1}=[[{' ' } request.cor_roicons'];  [request.cor_beh  cell(length(request.cor_beh),  length(request.cor_roicons))]];
for b=2:size(r_batcorr{k,1},1); for rc=2:size(r_batcorr{k,1},2)
        wbc.beh=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), r_batcorr{k,1}{b,1}), 1, 'first')));
        wbc.roiconbeta= cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), r_batcorr{k,1}{1,rc}))));
        
        if request.kendalls_correlation
            [wbc.r wbc.p]=corr(wbc.beh, wbc.roiconbeta,'type', 'Kendall');  
            if b==1 && rc==1; disp('Kendalls correlation'); end % 
%             wbc.stat=wbc.r;   % Print r or p?5
%             wbc.stat=wbc.p;
            wbc.stat=['ta=' num2str(wbc.r,2) ', p=' num2str(wbc.p,2)];
        else
            [wbc.r wbc.p]=corr(wbc.beh, wbc.roiconbeta);  %
            %         [wbc.r wbc.p]=corr(wbc.beh, wbc.roiconbeta,'type', 'Kendall');  disp('Kendalls correlation'); %
%             wbc.stat=wbc.r;   % Print r or p?5
%             wbc.stat=wbc.p;
            wbc.stat=['r=' num2str(wbc.r,2) ', p=' num2str(wbc.p,2)];
        end
        
        if wbc.p<0.001;  r_batcorr{k,1}{b,rc}= ['*** ' wbc.stat  ];;  % [num2str(wbc.stat,3) ' **'];
        elseif wbc.p<0.01;  r_batcorr{k,1}{b,rc}= ['** ' wbc.stat  ];;  % [num2str(wbc.stat,3) ' **'];
        elseif wbc.p<=0.05;  r_batcorr{k,1}{b,rc}= ['* ' wbc.stat ];  %[num2str(wbc.stat,3) ' *'];
        elseif wbc.p<0.1; r_batcorr{k,1}{b,rc}=wbc.stat ;
        else r_batcorr{k,1}{b,rc}=wbc.stat ;
        end
end,end
openvar r_batcorr{1,1}  % Put probabilities in 2nd col if necessary





% Load ad hoc file for comparison 

for o=1:1  % Load ci3 cm
% [w.n w.t d_nb]= xlsread('F:\1 [Context-Memory] fMRI Data\2 Second level results s4FullCardiacWithDerivAnts\m_ci3_ContextItemNomotor_Hit_Landmarks5\cm_m1_2x2\ROI\ci3 rois\(06-Aug-2015) Extracted betas.xlsx');
% d_ti=d_nb(1,:)';
% nb=  cell2mat(d_nb(2:end , find(strcmp(d_ti, 'HPC_aL_svtp-SimR')))) - cell2mat(d_nb(2:end , find(strcmp(d_ti, 'HPC_aL_svtp-SimN'))));
% d_nb=[d_nb [{'HPC_aL_svtp-Sim_vfx'}; num2cell(nb)]]; 
% nb=  cell2mat(d_nb(2:end , find(strcmp(d_ti, 'HPC_aL_svtp-SimRxCMem_Hit')))) - cell2mat(d_nb(2:end , find(strcmp(d_ti, 'HPC_aL_svtp-SimNxCMem_Hit'))));
% d_nb=[d_nb [{'HPC_aL_svtp-mSim_vfx'}; num2cell(nb)]]; 
% d_ti=d_nb(1,:)';
% 
% % roicon1='HPC_aL_svtp-SimRxCMem_Hit';
% % roicon1='HPC_aL_svtp-SimR';
% % roicon1='HPC_aL_svtp-Sim_vfx';
% roicon1='HPC_aL_svtp-mSim_vfx';

end

% [w.n w.t d_nb]= xlsread('F:\1 [Context-Memory] fMRI Data\2 Second level results s4FullCardiacWithDerivAnts\m_ci3 PPIs\PPI cell\cm co cell betas\cell anat\(25-Feb-2014) Extracted betas sph_SNVTA_R1_cmSRvSN cell.xlsx');

% [w.n w.t d_nb]= xlsread('F:\1 [Context-Memory] fMRI Data\2 Second level results s4FullCardiacWithDerivAnts\m_ci3 PPIs\PPI cell\cm co cell betas\cell anat\(25-Feb-2014) Extracted betas sph_SNVTA_R1_cmSRvSN cell.xlsx');
[w.n w.t d_nb]= xlsread('F:\1 [Context-Memory] fMRI Data\2 Second level results s4FullCardiacWithDerivAnts\m_ci3 PPIs\PPI cell\cm co cell betas\(21-Feb-2014) Extracted betas sph_SNVTA_R1_anat.xlsx');

d_nb(1,:)= cellfun(@(x)strtrim(x), d_nb(1,:), 'UniformOutput',0);
d_ti=d_nb(1,:)';
nb=  cell2mat(d_nb(2:end , find(strcmp(d_ti, 'SimR_HPC2anat_aL')))) - cell2mat(d_nb(2:end , find(strcmp(d_ti, 'SimN_HPC2anat_aL'))));
d_nb=[d_nb [{'HPC2anat_aL-Sim_vfx'}; num2cell(nb)]]; 
nb=  cell2mat(d_nb(2:end , find(strcmp(d_ti, 'SimR_HPC2anat_aR')))) - cell2mat(d_nb(2:end , find(strcmp(d_ti, 'SimN_HPC2anat_aR'))));
d_nb=[d_nb [{'HPC2anat_aR-Sim_vfx'}; num2cell(nb)]]; 
d_ti=d_nb(1,:)';


% roicon1= 'SimR_HPC2anat_aL'; 
% roicon1= 'SimR_HPC2anat_aR'; 
% roicon1= 'HPC2anat_aL-Sim_vfx';
roicon1= 'HPC2anat_aR-Sim_vfx'; 


roicon2='HPC_aL_svtp-SimRxCMem_Hit';
% % roicon2='HPC_aL_svtp-SimR';
% roicon2='HPC_aL_svtp-Sim_vfx';


b1 =  cell2mat(d_nb(2:end , find(strcmp(d_ti, roicon1))));
% b2 =  cell2mat(d_betas(2:end , find(strcmp(d_betas(1,:), [roi2 '-' con2]))));
b2 =  cell2mat(d_betas(2:end , find(strcmp(d_betas(1,:), roicon2))));
[r p]=corr(b1, b2); disp(['r=' num2str(r) '   p=' num2str(p)])
     
close all
f.scatter_dotsize=100;   f.scatter_linewidth=4;   f.FontSize=25; f.fontname='PT Sans Caption';   
figure('color','w','Position', [100 250 800 600])
scatter(b2,b1, f.scatter_dotsize,'LineWidth', 3); h=lsline; set(h,'LineWidth', f.scatter_linewidth)
xlabel(sprintf('Left DG/CA3 \nparameter estimate (univariate)'),'FontSize',f.FontSize, 'FontName', f.fontname)
ylabel(sprintf('Right anterior hippocampus \nparameter estimate (PPI)'),'FontSize',f.FontSize, 'FontName', f.fontname)
title([roicon1 '  ' roicon2 ':  r=' num2str(r) '   p=' num2str(p)])
set(gca,'FontSize',f.FontSize, 'FontName', f.fontname, 'LineWidth', 0.8);





error('Done with beh battery')

% % AD HOC PLOTS ------------------------------
% roicon= 'HPC_aL_tc-cF_Rej-Exp';
% beh='State anxiety scores';
% r={'HPC_aR_sc-cF_EVGain-Loss'}


% % roicon=r{1}; 
wr.rc= cell2mat( d_betas(2:end,  find(strcmp(d_betas(1,:),roicon)) ));  wr.b= cell2mat( d_betas(2:end,  find(strcmp(d_betas(1,:),beh)) ));
figure('color', 'w'), scatter(wr.rc, wr.b), lsline; xlabel(roicon,'FontSize',20), ylabel(sprintf(beh),'FontSize',20)
if      request.kendalls_correlation==0, [r p]= corr(wr.rc, wr.b); title(['r='  num2str(r) ', p=' num2str(p)],'FontSize',20)
else      [r p]= corr(wr.rc, wr.b,'type', 'Kendall'); title(['tau='  num2str(r) ', p=' num2str(p)],'FontSize',20)
end; set(gca,'FontSize',20)

% % -------------------------------------

% corrwhat={'rCA1_aL-cF_Rej-Exp'; 'rDG_aL-cF_Rej-Exp'; 'Trait anxiety scores'};
% 
% 
% d_corrwhat= cellfun(@(x)cell2mat(d_betas(2:end, strcmp(d_betas(1,:),x))), corrwhat, 'UniformOutput',0); d_corrwhat= [d_corrwhat{:}];
% [r p]=corr(d_corrwhat);  
% % [r p]=corr(d_corrwhat,'type', 'Kendall');  
% 
% disp([[{'r'}; corrwhat] [corrwhat'; num2cell(r)]]); disp(' ');disp(' ')
% disp([[{'p'}; corrwhat] [corrwhat'; num2cell(p)]]); disp(' ');disp(' ')

% 
% 
% 
% 
% 
% cons={'cF_Rej_vMargCho'; 'cF_NonRej_vMargCho'; 'ct_Rej_vMargCho'; 'ct_NonRej_vMargCho'};
% cons={'cF_vChosen';'cF_vBestUnchosen'; 'ct_vChosen';'ct_vBestUnchosen'}; 
% for c=1:length(cons)
%     wc.c1= cell2mat(d_betas(2:end, strcmp(d_betas(1,:), ['rCA1_aL-' cons{c}])));
%     wc.c2= cell2mat(d_betas(2:end, strcmp(d_betas(1,:), ['rCA3_aL-' cons{c}])));
%     
%     [h p]= ttest(wc.c1, wc.c2);
%     disp([cons{c} '   : p=' num2str(p)])   
% end
% 

 

 
%% Adhoc random plots/analysis etc

% % Plot imagesc for coactivation
% input('Select rois for imagesc-ing coactivation?'); log.rois=log.rois(1:end-2); log.n_rois=length(log.rois);
% log.rois_names={'BA46';'BA10';'Right striatum';'Superior MFG';'Precuneus'};
% close all hidden; figure; set(gcf,'color','w'); a=cell2mat(  r_batcorr{2,1}(2:end, 2:end));
% imagescnan(a, 'NanColor', [0.9 0.9 0.9]); axis square; colorbar;
% title('Parameter estimates for Explore choices', 'FontSize', 25)
% set(gca, 'FontSize', 15, 'YTick',1:log.n_rois, 'YTickLabel',log.rois_names,'XTick',1:log.n_rois, 'XTickLabel',log.rois_names)
% 
% % Correlations?
% d_exp=d_betas(:,82:end);
% corrcols=[51 81];
% [r p]=corr(cell2mat(d_betas(2:end, corrcols(1))), cell2mat(d_betas(2:end, corrcols(2)))); r,p

%% Targetted/requested correlations + plots (r_corr)

r_corr=cell(0, 5); k=1;% (1) Beta name, (2) Beh name, (3) corr r, (4) corr p, (5) data, (6) corr sig?

request.plotcorr={'o2'}; 

for o1=1:1 % [Choice & c13 rois: HPC TxC Beh inhibition]
    
    if sum(strcmp('o1', request.plotcorr))==1
        
        % LEFT HPC, cF Rej-Exp ##############################
        wc.roiname='HPC_aL_stc'; wc.roiconname='Parameter estimates from Left Inferior\n Hippocampus(Exp, Reject > Explore)';
        wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.Trait';   wc.behname='Trait anxiety scores';  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.State'; wc.behname='State anxiety scores'; 
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        % RIGHT HPC, cF Rej-Exp ##############################
        wc.roiname='HPC_aR_stc'; wc.roiconname='Parameter estimates from Right Inferior\n Hippocampus(Exp, Reject > Explore)';
        wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.behname= 'st.State';   wc.behname='State anxiety scores'; 
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        
        
        % Superior HPC rois 
        
        % LEFT HPC, cF Rej-Exp ##############################
        wc.roiname='HPC_aL_sc'; wc.roiconname='Parameter estimates from Left Superior\n Hippocampus(Exp, Reject > Explore)';
        wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.Trait';   wc.behname='Trait anxiety scores';  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.State'; wc.behname='State anxiety scores'; 
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        % RIGHT HPC, cF Rej-Exp ##############################
        wc.roiname='HPC_aR_sc'; wc.roiconname='Parameter estimates from Right Superior\n Hippocampus(Exp, Reject > Explore)';
        wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.behname= 'st.State';   wc.behname='State anxiety scores'; 
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        
    
    
    end
end
 

% r_corr=r_corr(1:4,:);

% r_corr=r_corr(4:6,:);
% close all hidden
% Automated plot all requested correlations
disp('Targeted/requested correlations for plotting:'); 
f.scatter_dotsize=100;   f.scatter_linewidth=4;   f.FontSize=15; f.fontname='PT Sans Caption';   r_corr{:,1}
f.plotcols=4;  f.figwidth= f.plotcols*400; f.figheight=ceil(size(r_corr,1)/ (f.plotcols))*400;
f.subplot_VerHorz=[0.2 0.1]; f.fig_BotTop=[0.15 0.1]; f.fig_LeftRight=[0.1 0.05]; k=1;
figure('Name', 'Beh - beta ', 'NumberTitle', 'off', 'Position', [100 250 f.figwidth f.figheight], 'Color', 'w');  
for c=1:size(r_corr,1)
    disp(r_corr{k,1})
    subtightplot(ceil(size(r_corr,1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
%     scatter(r_corr{c,5}(:,1),  r_corr{c,5}(:,2)); lsline
    scatter(r_corr{c,5}(:,1),  r_corr{c,5}(:,2), f.scatter_dotsize,'LineWidth', 3); h=lsline; set(h,'LineWidth', f.scatter_linewidth)
    
    
    
    ylabel(r_corr{c,2}, 'FontSize',f.FontSize, 'FontName', f.fontname);  
%     ylabel(sprintf(r_corr{c,2}), 'FontSize',f.FontSize, 'FontName', f.fontname)
%     xlabel(sprintf(['Parameter estimates from\n' r_corr{c,1}]), 'FontSize',15, 'FontName', f.fontname)
    xlabel(sprintf(r_corr{c,1}), 'FontSize',f.FontSize, 'FontName', f.fontname)
    set(gca,'FontSize',f.FontSize, 'FontName', f.fontname, 'LineWidth', 0.8);
    if r_corr{c,6}==1; 
        if request.kendalls_correlation, title(['tau='  num2str(r_corr{c,3},2) ', p=' num2str(r_corr{c,4},2) ' *'], 'FontSize',15, 'FontName', f.fontname)    
        else title(['r='  num2str(r_corr{c,3},2) ', p=' num2str(r_corr{c,4},2) ' *'], 'FontSize',15, 'FontName', f.fontname)
        end
    else
        if request.kendalls_correlation, title(['tau='  num2str(r_corr{c,3},2) ', p=' num2str(r_corr{c,4},2)], 'FontSize',15, 'FontName', f.fontname) 
        else title(['r='  num2str(r_corr{c,3},2) ', p=' num2str(r_corr{c,4},2)], 'FontSize',15, 'FontName', f.fontname)
        end
    end
    
    
% %     Extra adhoc stuff
%     xlim([-1 2])
xlim('auto')
    
    % 
    k=k+1;
end

  
