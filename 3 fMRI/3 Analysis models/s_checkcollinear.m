% s_checkcollinearity
clear all, close all hidden, clc

% Requested analysis
log.specificsubjects={}; 
% log.specificsubjects={'p01_CW'}; 

log.plot_ranger=[];
% log.plot_ranger=[0 0.5];  % Plot range of correlations
 

for o1=1:1 % General settings and specifications


    % Subjects
    where.where='D:\Dropbox\SCRIPPS\1 ContextMem fMRI';  where.data_brain='D:\1 [Context-Memory] fMRI Data\1 MRI data';  where.data_beh=[where.where filesep '2 Behavioural data'];
    addpath(where.where)
    log.w=load([where.data_brain filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    w.onsets=load([where.data_brain filesep 'Onsetslog_m_c1_Contextall_Hit' ]); % Apply subject filter first according to validity of onsets file
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog,log.specificsubjects,vertcat({'Subject' 'onsets'}, w.onsets.subjectsok),'onsets');    
    
end


%% For mixed Context + Item models (where both context mem and item mem are represented)

dothis=1;
if dothis
    
    req.testwhat='cm_iom_out'; 
%     req.testwhat='co_cm_iom_out'; 
%     req.testwhat='co_cm_iom'; 
%     req.testwhat='co_cm_iom_hm'; 
%     req.testwhat='co_iom';  
%     req.testwhat='cop_iom';  
%     req.testwhat='co_cm_iop';
%     req.testwhat='co_cm_io';
    req.testwhat='all'; 
    
    % % WHICH MODEL?
    log.FLthread='2 First level s4FullCardiacWithDerivAnts';
    log.FLmod='m_ci3_ContextItemNomotor_Hit';
%     log.FLmod='m_ci2_ContextonlyItem_Hit';
%     log.FLmod='m_c5_ContextonlyItempresonly';
%     log.FLmod='m_i4_ItemContextevent_Hit'; 
%     log.FLmod='m_c4_ContextallItempresonly_Hit'; % req.testwhat='co_cm_iop';  
%     log.FLmod='m_i3_ItemContextpresonly_Hit';  % req.testwhat='cop_iom';  
%     log.FLmod='m_i6_ItemMempmodContextpresonly_Hit'; %  req.testwhat='co_cm_iomp';  
%     log.FLmod='m_ci6_ContextItempmodWithmotor_Hit'; req.testwhat='all';     cells ={'SimR' 'SimRxCMem_Hit' 'SimN' 'SimNxCMem_Hit' 'DisR' 'DisRxCMem_Hit' 'DisN' 'DisNxCMem_Hit' 'SR_Item' 'SR_ItemxMem_Hit' 'SN_Item' 'SN_ItemxMem_Hit' 'DR_Item' 'DR_ItemxMem_Hit' 'DN_Item' 'DN_ItemxMem_Hit' 'Outcome' 'OutcomexOutcome_shown' 'Motor'};
%         log.FLmod='m_c7_ContexteventItempresonly_Hit'; 


    % Fetch design matrices
    d_mat=cell(log.n_subjs,2);  % Design matrix, Correlation matrix
    for s=1:log.n_subjs
        ws.s=load([where.data_brain filesep log.subjects{s} filesep log.FLthread filesep log.FLmod ' Contrasted\SPM.mat']);
%         ws.s=load([where.data_brain filesep log.subjects{s} filesep log.FLthread filesep log.FLmod ' Estimated\SPM_' log.FLmod '.mat']); disp('loadig wrong thing !')
        

        if s==1; log.RegNames= ws.s.SPM.xX.name(1:3:end)'; log.RegNames = log.RegNames(1:find(strcmp(log.RegNames, 'Sn(1) R1'))-1); log.n_regs= length(log.RegNames);  end
        d_mat{s,1} = ws.s.SPM.xX.X(:, 1:3:end);
        d_mat{s,1} = d_mat{s,1}(:, 1: log.n_regs) ;
%         log.RegNames
%         error
        
        % Which regressors to include? 
        if strcmp(req.testwhat,'cm_iom_out')==1, 
            d_mat{s,1}= d_mat{s,1}(:,  [find(strcmp(log.RegNames, 'Sn(1) SimRxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SR_Hit*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SimNxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SN_Hit*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisRxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DR_Hit*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisNxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DN_Hit*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) Outcome*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) OutcomexOutcome_shown^1*bf(1)')) ]);
        elseif strcmp(req.testwhat,'co_cm_iom_out')==1, 
            d_mat{s,1}= d_mat{s,1}(:,  [find(strcmp(log.RegNames, 'Sn(1) SimR*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimRxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SR_Hit*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimN*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimNxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SN_Hit*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisR*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisRxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DR_Hit*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) DisN*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisNxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DN_Hit*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) Outcome*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) OutcomexOutcome_shown^1*bf(1)')) ]);
        elseif strcmp(req.testwhat,'co_cm_iom')==1, 
%             d_mat{s,1}= d_mat{s,1}(:,  [find(strcmp(log.RegNames, 'Sn(1) SimR*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimRxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SR_Hit*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimN*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimNxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SN_Hit*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisR*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisRxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DR_Hit*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) DisN*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisNxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DN_Hit*bf(1)'))]);
            d_mat{s,1}= d_mat{s,1}(:,  [find(strcmp(log.RegNames, 'Sn(1) SimR*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimRxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SR_Miss*bf(1)'))   find(strcmp(log.RegNames, 'Sn(1) SimN*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimNxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SN_Miss*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) DisR*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisRxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DR_Miss*bf(1)'))   find(strcmp(log.RegNames, 'Sn(1) DisN*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisNxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DN_Miss*bf(1)'))]); disp('looking at mem MISSES'); 
        elseif strcmp(req.testwhat,'co_cm_iomp')==1, 
%             d_mat{s,1}= d_mat{s,1}(:,  [find(strcmp(log.RegNames, 'Sn(1) SimR*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimRxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SR_Hit*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimN*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimNxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SN_Hit*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisR*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisRxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DR_Hit*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) DisN*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisNxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DN_Hit*bf(1)'))]);
            disp('Use all regs!')
        elseif strcmp(req.testwhat,'co_cm_iom_hm')==1, 
            d_mat{s,1}= d_mat{s,1}(:,  [find(strcmp(log.RegNames, 'Sn(1) SimR*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimRxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SR_Hit*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SR_Miss*bf(1)'))   find(strcmp(log.RegNames, 'Sn(1) SimN*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimNxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SN_Hit*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SN_Miss*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) DisR*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisRxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DR_Hit*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) DR_Miss*bf(1)'))   find(strcmp(log.RegNames, 'Sn(1) DisN*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisNxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DN_Hit*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DN_Miss*bf(1)'))]); 
        elseif strcmp(req.testwhat,'co_iom')==1, 
            d_mat{s,1}= d_mat{s,1}(:,  [find(strcmp(log.RegNames, 'Sn(1) SimR*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SR_Hit*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SimN*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SN_Hit*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisR*bf(1)'))   find(strcmp(log.RegNames, 'Sn(1) DR_Hit*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisN*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DN_Hit*bf(1)'))]);
        elseif strcmp(req.testwhat,'cop_iom')==1, 
            d_mat{s,1}= d_mat{s,1}(:,  [find(strcmp(log.RegNames, 'Sn(1) ContextPresonly*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SR_Hit*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SN_Hit*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DR_Hit*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DN_Hit*bf(1)'))]);
        elseif strcmp(req.testwhat,'co_cm_iop')==1,  
%             d_mat{s,1}= d_mat{s,1}(:,  [find(strcmp(log.RegNames, 'Sn(1) SimR*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimRxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SimN*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimNxCMem_Hit^1*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) DisR*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisRxCMem_Hit^1*bf(1)'))   find(strcmp(log.RegNames, 'Sn(1) DisN*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisNxCMem_Hit^1*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) Outcome*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) OutcomexOutcome_shown^1*bf(1)'))      ]);
            d_mat{s,1}= d_mat{s,1}(:,  [find(strcmp(log.RegNames, 'Sn(1) SimR*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimRxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SimN*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimNxCMem_Hit^1*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) DisR*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisRxCMem_Hit^1*bf(1)'))   find(strcmp(log.RegNames, 'Sn(1) DisN*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisNxCMem_Hit^1*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) Items*bf(1)'))   ]);
        elseif strcmp(req.testwhat,'co_cm_io')==1,  
            d_mat{s,1}= d_mat{s,1}(:,  [find(strcmp(log.RegNames, 'Sn(1) SimR*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimRxCMem_Hit^1*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) SR_Item*bf(1)'))   find(strcmp(log.RegNames, 'Sn(1) SimN*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SimNxCMem_Hit^1*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) SN_Item*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) DisR*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisRxCMem_Hit^1*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) DR_Item*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) DisN*bf(1)')) find(strcmp(log.RegNames, 'Sn(1) DisNxCMem_Hit^1*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) DN_Item*bf(1)'))   ]);
        end 
        
        ws=[];
    end
    
    % PLOT all
    if strcmp(req.testwhat,'cm_iom_out')==1,                cells={ 'SRc' 'SRi' 'SNc' 'SNi' 'DRc' 'DRi' 'DNc' 'DNi' 'Out' 'OutPres'}; 
    elseif strcmp(req.testwhat,'co_cm_iom_out')==1,         cells={ 'SRco' 'SRc' 'SRi' 'SNco' 'SNc' 'SNi' 'DRco' 'DRc' 'DRi'  'DNco' 'DNc' 'DNi' 'Out' 'OutPres'}; 
	elseif strcmp(req.testwhat,'co_cm_iom')==1,         cells={ 'SRco' 'SRc' 'SRi' 'SNco' 'SNc' 'SNi' 'DRco' 'DRc' 'DRi'  'DNco' 'DNc' 'DNi'};     
    elseif strcmp(req.testwhat,'co_cm_iomp')==1,          cells={ 'Ctx'  'SRi' 'SRim' 'SNi' 'SNim'     'DRi' 'DRim' 'DNi' 'DNim'     'Out' 'OutPres'};
	elseif strcmp(req.testwhat,'co_cm_iom_hm')==1,      cells={ 'SRco' 'SRc' 'SRih' 'SRim' 'SNco' 'SNc' 'SNih' 'SNim' 'DRco' 'DRc' 'DRih'  'DRim'  'DNco' 'DNc' 'DNih' 'DNim'}; 
	elseif strcmp(req.testwhat,'co_iom')==1,         cells={ 'SRco' 'SRi' 'SNco'  'SNi' 'DRco'  'DRi'  'DNco' 'DNi'}; 
    elseif strcmp(req.testwhat,'cop_iom')==1,         cells={ 'Ctx' 'SRi'  'SNi'  'DRi'  'DNi'};
    elseif strcmp(req.testwhat,'co_cm_iop')==1, cells={'SRco' 'SRc' 'SNco' 'SNc'  'DRco' 'DRc'  'DNco' 'DNc' 'item'}; 
    elseif strcmp(req.testwhat,'co_cm_io')==1,  cells={ 'SRco' 'SRc' 'SRi' 'SNco' 'SNc' 'SNi' 'DRco' 'DRc' 'DRi'  'DNco' 'DNc' 'DNi'}; 
    elseif strcmp(req.testwhat,'co_cm')==1, cells={'SRco' 'SRc' 'SNco' 'SNc'  'DRco' 'DRc'  'DNco' 'DNc' 'item'}; 
    elseif strcmp(req.testwhat,'all')==1,
        disp('cell names not specified!') 
    end
    if exist('cell','var')==0; cells = cellfun(@(x)x(7:end-6), log.RegNames, 'UniformOutput',0);  end 
    
    d_corrs=nan(log.n_subjs,5);  % Mean, min, max, absmean, absmax
    close all hidden, f.plotcols=6;  f.figwidth= 2400; f.figheight=800; f.fontsize=10; f.fontname='PT Sans Caption';  f.subplot_VerHorz=[0.05 0.03]; f.fig_BotTop=[0.05 0.15]; f.fig_LeftRight=[0.05 0.05];
    figure('Name', 'Reg Collinearity', 'NumberTitle', 'off', 'Position', [100 50 f.figwidth f.figheight], 'Color', 'w');  k=1;
    for s=1:log.n_subjs
        
        % Correlate
        [ws.r ws.p]=corr(d_mat{s,1});
        d_mat{s,2} =ws.r ; ws.r(ws.r==1)=nan;
        d_corrs(s,:)= [nanmean(ws.r(:)) min(ws.r(:)) max(ws.r(:)) nanmean(abs(ws.r(:))) max(abs(ws.r(:)))];
        d_outcomecorr(s,1:2)=  [min(ws.r(:)) max(ws.r(:))]; % Greatest anticorr and corr 
%         ws.r(4:end, 1:3)=nan; % Blot out cross-cell correlations (i.e. SR vs SN)
%         ws.r([1:3 7:end], 4:6)=nan;
%         ws.r([1:6 10:end], 7:9)=nan;
%         ws.r(1:9, 10:12)=nan;
%         d_outcomecorr(s,1:2)=  [min(ws.r(:)) max(ws.r(:))]; % Greatest anticorr and corr 
        
        
        % PLOT
        subtightplot(ceil((log.n_subjs)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
        imagescnan(d_mat{s,2}, 'nancolor',[0.1 0.1 0.1]) % , axis square   %,   title(log.subjects{s})
        colorbar
        set(gca, 'FontSize',10, 'XTick', 1:length(cells), 'XTickLabel', cells, 'YTick', 1:length(cells), 'YTickLabel', cells);  % xlim([0.3 2.7])

        
        
        if (k-1)/f.plotcols==floor((k-1)/f.plotcols), colorbar, end
        if ~isempty(log.plot_ranger), caxis(log.plot_ranger), end
    end
    % disp([num2cell(1:length(log.RegNames))' log.RegNames'] )
%     mean(d_corrs)
    
    figure('color','w'), subplot(1,2,1), barwitherr(std(d_corrs)./sqrt(log.n_subjs), mean(d_corrs), 'y'), ylabel('r correlation',  'FontSize', 20), set(gca,  'FontSize', 20, 'XTick', 1:5, 'XTickLabel', {'Mean' 'min' 'max' 'Abs mean' 'Abs Max'})
    d_mat{log.n_subjs+1,2}    =zeros(size(d_mat{1,2}));
    for s=1:log.n_subjs, d_mat{log.n_subjs+1,2}=d_mat{log.n_subjs+1,2}+ d_mat{s,2}    ; end
    d_mat{log.n_subjs+1,2}=d_mat{log.n_subjs+1,2}./log.n_subjs;
    subplot(1,2,2),imagescnan(d_mat{log.n_subjs+1,2}, 'nancolor',[0.1 0.1 0.1]), colorbar, axis square   %,   title(log.subjects{s})
    set(gca, 'FontSize',10, 'XTick', 1:length(cells), 'XTickLabel', cells, 'YTick', 1:length(cells), 'YTickLabel', cells);  % xlim([0.3 2.7])
    xticklabel_rotate
end

% title('Context event Context memory Object (recognized)', 'FontSize',20);
% % axis off
% % subplot(1,2,2)
% % % 
% a={'Context event';'Context memory';'Recognized object'}
% a= [a;a;a;a]
% 
% set(gca, 'FontSize',17, 'XTick', 1:length(cells), 'XTickLabel', cells, 'YTick', 1:length(cells), 'YTickLabel', a);  % xlim([0.3 2.7])
% 

abs(d_outcomecorr)>0.6


mean( d_outcomecorr(:,2)) 
std( d_outcomecorr(:,2)) 




%%

