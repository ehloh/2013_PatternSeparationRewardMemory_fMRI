% evaluate learning 
%   RT effect in conditioning
%   Item judgments in Encoding (Accuracy & RT)
clear all;close all;clc

where.data='I:\2 fMRI behaviour analysis\2 Behavioural Data';
% where.data='/Volumes/MEDIASTAR/1 fMRI/Data';
where.script='I:\2 fMRI behaviour analysis\1 Analysis scripts';
% where.script='/Volumes/MEDIASTAR/2 fMRI behaviour analysis/1 Analysis scripts';
where.printout='I:\2 fMRI behaviour analysis\3 Analysis inputs';
% where.printout='/Volumes/MEDIASTAR/2 fMRI behaviour analysis/3 Analysis inputs';

% Load up
cd(where.data)
w.dir=dir('p*');
% w.dir=dir(where.data);
load([where.script filesep 'i_eval0_learning.mat'])

for i=1:size(instruc,1) % Data headers for results files
    res_cond{1,1+i}=['cRT_' instruc{i,1}];
    res_cond{1,1}='Subject';
    res_encJud{1,1+i}=['eJ_' instruc{i,1}];
    res_encJud{1,1}='Subject';
    res_encRT{1,1+i}=['eRT_' instruc{i,1}];
    res_encRT{1,1}='Subject';
end

for s=1:size(w.dir,1) % insert subject loop stuff here
ws.subjname=w.dir(s).name;
ws.subjname=w.dir(s).name;
disp(['Subject ' num2str(s) ' - ' ws.subjname])

ws.cond_a=load([where.data filesep ws.subjname filesep ws.subjname '_file_1conditioning.mat']);
ws.cond=ws.cond_a.conditioning.trialstats;
ws.encb1=load([where.data filesep ws.subjname filesep ws.subjname '_file_2encodingFMRI_b1.mat']);
ws.encb2=load([where.data filesep ws.subjname filesep ws.subjname '_file_2encodingFMRI_b2.mat']);
ws.enc=vertcat(ws.encb1.encoding.trialstats,ws.encb2.encoding.trialstats);
res_cond{s+1,1}=ws.subjname;
res_encJud{s+1,1}=ws.subjname;
res_encRT{s+1,1}=ws.subjname;

for o1=1:1% Sort into cell samples
% Conditioning
ws.cond=ws.cond(ws.cond(:,8)>0,:); % Get rid of premature responses
ws.c.s1=ws.cond((ws.cond(:,3)==1 & ws.cond(:,4)==1), :); % sample 1= SimR
ws.c.s2=ws.cond((ws.cond(:,3)==1 & ws.cond(:,4)==0), :); % sample 2= SimN
ws.c.s3=ws.cond((ws.cond(:,3)==2 & ws.cond(:,4)==1), :); % sample 3= DisR
ws.c.s4=ws.cond((ws.cond(:,3)==2 & ws.cond(:,4)==0), :); % sample 4= DisN
ws.c.s0=[];
% Encoding
ws.e.s1=ws.enc((ws.enc(:,3)==1 & ws.enc(:,4)==1), :); % sample 1= SimR
ws.e.s2=ws.enc((ws.enc(:,3)==1 & ws.enc(:,4)==0), :); % sample 2= SimN
ws.e.s3=ws.enc((ws.enc(:,3)==2 & ws.enc(:,4)==1), :); % sample 3= DisR
ws.e.s4=ws.enc((ws.enc(:,3)==2 & ws.enc(:,4)==0), :); % sample 4= DisN
ws.e.s0=[];
end

% Data headers
for i=1:size(instruc,1)
    res_cond{1,1+i}=['cRT_' instruc{i,1}];
    res_encJud{1,1+i}=['eJ_' instruc{i,1}];
    res_encRT{1,1+i}=['eRT_' instruc{i,1}];
end

%% Analysis starts

for i=1:size(instruc,1) % Calculate conditioning RT & Enc accuracy 
    % Conditioning
    eval(['wc.c_dat=vertcat(ws.c.s' num2str(instruc{i,2}(1)) ', ws.c.s' num2str(instruc{i,2}(2)) ', ws.c.s' num2str(instruc{i,2}(3)) ', ws.c.s' num2str(instruc{i,2}(4)) ');'])
    res_cond{s+1, i+1}=mean(wc.c_dat(:,8));
    % Encoding accuracy
    eval(['wc.e_dat=vertcat(ws.e.s' num2str(instruc{i,2}(1)) ', ws.e.s' num2str(instruc{i,2}(2)) ', ws.e.s' num2str(instruc{i,2}(3)) ', ws.e.s' num2str(instruc{i,2}(4)) ');'])
    res_encJud{s+1, i+1}=mean([wc.e_dat(:,26)' wc.e_dat(:,27)' wc.e_dat(:,28)']');
    wc=[];
end

% Re-compile encoding RT data, excluding items taht were incorrectly judged
ws.er.s1=vertcat(ws.e.s1(ws.e.s1(:,26)==1, 21),ws.e.s1(ws.e.s1(:,27)==1, 23),ws.e.s1(ws.e.s1(:,28)==1, 25));
ws.er.s2=vertcat(ws.e.s2(ws.e.s2(:,26)==1, 21),ws.e.s2(ws.e.s2(:,27)==1, 23),ws.e.s2(ws.e.s2(:,28)==1, 25));
ws.er.s3=vertcat(ws.e.s3(ws.e.s3(:,26)==1, 21),ws.e.s3(ws.e.s3(:,27)==1, 23),ws.e.s3(ws.e.s3(:,28)==1, 25));
ws.er.s4=vertcat(ws.e.s4(ws.e.s4(:,26)==1, 21),ws.e.s4(ws.e.s4(:,27)==1, 23),ws.e.s4(ws.e.s4(:,28)==1, 25));
ws.er.s0=[];

for i=1:size(instruc,1) % Calculate Enc RT
    eval(['wc.e_dat=vertcat(ws.er.s' num2str(instruc{i,2}(1)) ', ws.er.s' num2str(instruc{i,2}(2)) ', ws.er.s' num2str(instruc{i,2}(3)) ', ws.er.s' num2str(instruc{i,2}(4)) ');'])
    res_encRT{s+1, i+1}=mean([wc.e_dat(:)' wc.e_dat(:)' wc.e_dat(:)']');
    wc=[];
end

ws=[];
end

%%  Statistical analysis?



%% Print out 
res=horzcat(res_cond, res_encJud(:,2:size(res_encJud,2)),res_encRT(:,2:size(res_encRT,2)));
printreport=print2txt(where.printout, ['(' date ') Learning results'], res);
printreport