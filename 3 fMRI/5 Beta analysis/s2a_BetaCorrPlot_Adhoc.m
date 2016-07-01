behname=[]; betname=[]; 
% % AD HOC PLOTS ------------------------------

beh='eRT_cell.Sim_cR';
beh='eRT_cell.Sim_valfx';
% beh='dpr_cell.Sim_valfx'; 
% beh='dpr_cell.R_simfx';


% beh='DG_CA2_3_L1-SimR'; 
beh='eRT_cell.Sim_valfx'; rr={'DG_CA2_3_R1-SimR'}; 
behname='RT speeding in the similar condition (ms)'; betname='PPI Parameter estimates'; 


% ####################
f.scatter_dotsize=100;   f.scatter_linewidth=4;   f.FontSize=25; f.fontname='PT Sans Caption';   
%
roicon=rr{1}; wr.rc= cell2mat( d_betas(2:end,  find(strcmp(d_betas(1,:),roicon)) ));  wr.b= cell2mat( d_betas(2:end,  find(strcmp(d_betas(1,:),beh)) ));
figure('color', 'w')
if isempty(behname)==1, behname=beh; end; if isempty(betname)==1, behname=roicon; end
scatter(wr.rc, wr.b, f.scatter_dotsize,'LineWidth', 3); h=lsline; set(h,'LineWidth', f.scatter_linewidth); xlabel(sprintf(betname),'FontSize',f.FontSize), ylabel(sprintf(behname),'FontSize',f.FontSize)
if      request.kendalls_correlation==0, [r p]= corr(wr.rc, wr.b); title(['r='  num2str(r) ', p=' num2str(p)],'FontSize',f.FontSize)
else      [r p]= corr(wr.rc, wr.b,'type', 'Kendall'); title(['tau='  num2str(r) ', p=' num2str(p)],'FontSize',f.FontSize)
end; set(gca,'FontSize',f.FontSize)


openvar wr.rc,

error('done! :)')
% 
% 
% xlabel(sprintf('Difference in left hippocampal parameter estimates\n  on Non-Reject trials (Exp>Ctrl)'), 'FontSize',f.FontSize,'FontName', 'PT Sans Caption')
% xlabel(sprintf('Difference in right hippocampal parameter \n  estimates (Exp, Non-Reject>Reject)'), 'FontSize',f.FontSize,'FontName', 'PT Sans Caption')
% ylabel(sprintf(beh),'FontSize',f.FontSize)

request.kendalls_correlation=0; 
request.kendalls_correlation=1; 


%% MEDIAN split 


% split_iv={'Probability of following null information\n from Exploring in experimental task'}; 
split_iv={'State anxiety scores'};

wm.ivbeh=cell2mat(d_betas(2:end, strcmp(d_betas(1,:), split_iv))) ; 
sum(wm.ivbeh >median(wm.ivbeh))
wm.subs_high= find(wm.ivbeh > median(wm.ivbeh));   % HIGH score group 
wm.subs_low= find(1- (wm.ivbeh >median(wm.ivbeh)));   % HIGH score group 


wm.subs_low



