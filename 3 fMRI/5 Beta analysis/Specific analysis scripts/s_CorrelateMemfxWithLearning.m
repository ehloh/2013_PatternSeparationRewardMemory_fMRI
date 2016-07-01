% Correlate memory effect with learning & other variables, across all subjects
clear all; close all hidden; clc

where.where='I:\1 fMRI analysis'; 
% where.where='/Volumes/PENNYDISK/1 fMRI analysis';

% Requested analysis
log.specificsubjects={};
%
log.behlearn='(15-8-13) Learning data'; % mat
log.behmem='(23-Sep-2013) Memory scores all'; % txt

for o1=1:1 % General setup 
    
    % Load subjects
    addpath(where.where);
    where.data_beh=[where.where filesep '2 Behavioural data']; where.data_behgroup=[where.data_beh filesep 'Group behaviour files'];
    log.w=load([where.data_beh filesep 'datalog_allsubjects.mat']); log.datalog=log.w.datalog;
    [log.subjects_all log.n_subjs_all log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
    % Include all subjects
    log.subjects=log.subjects_all ;
    log.n_subjs=log.n_subjs_all;
end

%%

% Load all data
w.b=importdata([where.data_behgroup filesep log.behmem '.txt']);
d.mem=[strtrim(w.b.textdata(:,1)) [strtrim(w.b.textdata(1,2:end)); num2cell(w.b.data)]];
w.b=load([where.data_behgroup filesep log.behlearn '.mat']);
d.learn=w.b.ldat;

% Checks
d.learn(2:end,:)=sortrows(d.learn(2:end,:),1);
d.mem(2:end,:)=sortrows(d.mem(2:end,:),1);

% Exclude subjects?
exclude={'p02_MK';'p13_CB'};
for i=1:length(exclude)
    d.mem(find(strcmp(d.mem(:,1), exclude{i})),:)=[];
    d.learn(find(strcmp(d.learn(:,1), exclude{i})),:)=[];
end

%
if mean(strcmp(d.learn(2:end,1), d.mem(2:end,1)))~=1;  error('Subjects mismatch'); end
d.learn(:,1)=[];
d.mem(:,1)=[];

%
log.mem_list=strtrim(d.mem(1,:)');
log.learn_list=strtrim(d.learn(1,:)');
d.mem(1,:)=[];d.learn(1,:)=[];

%% Instructions: which correlations?

log.corr_learn={'eJ_cell.Sim_cR';'eJ_cell.Sim_cN';'eJ_cell.Dis_cR';'eJ_cell.Dis_cN';};
% log.memscores={'Cell.Sim_R';'Cell.Sim_N';'Cell.Dis_R';'Cell.Dis_N';};
% log.learn_list


%
% log.corr_learn={'eJ_cell.SimValfx';'eJ_cell.RewSimfx'};
log.memtype={'md';'mh';'msh';'mr';'mk';'msr';'msk'};
log.memscores={'cell.R_simfx';'cell.Sim_valfx'};
log.corr_mem=cell(length(log.memscores)*length(log.memtype),1); k=1;
for i=1:length(log.memtype)
    for j=1:length(log.memscores)
        log.corr_mem{k}=[log.memtype{i} '.' log.memscores{j}]; k=k+1;
    end
end

%% Perform correlations, r_corr
% [1,2] variables [3] r [4] p [5] sig type

r_corr=cell(length(log.corr_learn)*length(log.corr_mem),5); k=1;
for m=1:length(log.corr_mem)
    for l=1:length(log.corr_learn)
        
        % Checks
        if sum(strcmp(log.mem_list, log.corr_mem{m}))~=1; error('Ill-specified mem score for correlation'); end
        if sum(strcmp(log.learn_list, log.corr_learn{l}))~=1; error('Ill-specified learning score for correlation'); end
        
        % Correlate
        [wc.r wc.p]=corrcoef(cell2mat(d.mem(:,find(strcmp(log.mem_list, log.corr_mem{m})))), cell2mat(d.learn(:,find(strcmp(log.learn_list, log.corr_learn{l})))));
        
            r_corr{k,1}=log.corr_mem{m};
            r_corr{k,2}=log.corr_learn{l};
            r_corr{k,3}=wc.r(2,1);
            r_corr{k,4}=wc.p(2,1);
        if wc.p<0.1    
            % Mark significance
            if wc.p<0.001
                r_corr{k,5}='     * * *';
            elseif wc.p<0.01
                r_corr{k,5}='     * *';
            elseif wc.p<0.05
                r_corr{k,5}='     *';
            end
        end
        
        
            % 
            k=k+1;
    end
end
    
    

openvar r_corr
% absolutely nothing here!



