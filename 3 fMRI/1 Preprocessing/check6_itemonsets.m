clear all; close all hidden; clc

subj='p01_CW';
memtype='Roc'; % 'Roc' 'Hit' 'Rem'
model='m_ci5_ContextItempmod';


where.databrain='C:\Users\eloh\Desktop\1 [Context-Memory] fMRI Data\1 MRI data';
where.databeh='I:\1 fMRI analysis\2 Behavioural data';
where.encscript='I:\3 fMRI study documentation\1 Data acquisition scripts';
% where.databeh='/Volumes/PENNYDISK/1 fMRI analysis/2 Behavioural data/';
% where.databrain='/Volumes/SANDISK/1 CONTEXT Brain data/1 MRI data';
% where.encscript='/Volumes/PENNYDISK/3 fMRI study documentation/1 Data acquisition scripts'; 
w.s=load([where.encscript filesep 'Stimuli' filesep 'stimlist.mat']); itemlist=w.s.itemlist; w.m=load([where.databeh filesep  subj filesep subj '_file_3memorytest.mat']);

%% Display (check against cogent log): each presented item, Sim x Val x Mem category

a=w.m.memtest.trialstats(w.m.memtest.trialstats(:,3)==1,:);
for i=1:size(a,1)
    showitem{i,1}=a(i,6);
    showitem{i,2}=a(i,2);
    showitem{i,2}=itemlist{a(i,2),1};
    if a(i,3)==1
        showitem{i,3}='Sim';
    else 
        showitem{i,3}='Dis';
    end
    if a(i,4)==1
        showitem{i,3}=[showitem{i,3} ' Reward'];
    else 
        showitem{i,3}=[showitem{i,3} ' Neutral'];
    end
    
    if a(i,12)==1
        showitem{i,4}='Hit';
    else 
        showitem{i,4}='Miss';
    end
    
    %
%     item{i,1}=a(i,6);
    item(i,2)=a(i,2); % 2=Item num
    item(i,3)=a(i,12); % Old/New
    item(i,4)=a(i,14); % Rem/Know
    item(i,5)=a(i,16); % Sure/Guess
    switch [num2str(item(i,3)) ' ' num2str(item(i,4)) ' ' num2str(item(i,5))]
        case '2 999 1'  % Sure New (1)
            item(i,6)=1;
        case '2 999 2'  % Guess New (2)
            item(i,6)=2;
        case '1 2 2'  % Guess Know (3)
            item(i,6)=3;
        case '1 2 1'  % Sure Know (4)
            item(i,6)=4;
        case '1 1 2'  % Guess Rem (5)
            item(i,6)=5;
        case '1 1 1'  % Sure Rem (6)
            item(i,6)=6;
        otherwise
            error([num2str(item(i,3)) ' ' num2str(item(i,4)) ' ' num2str(item(i,5))])
    end
end
showitem=sortrows(showitem,1);
    
% error('Stop here')


%% Compile all stats in fMRI-timing format

e1=load([where.databeh filesep subj filesep subj '_file_2encodingFMRI_b1.mat']); start{1}=e1.encoding.times.start; e1=e1.encoding.trialstats;
e2=load([where.databeh filesep subj filesep subj '_file_2encodingFMRI_b2.mat']);  start{2}=e2.encoding.times.start; e2=e2.encoding.trialstats;

for o1=1:1 % dataspecs
ocol.SimDis=3; % original data set
ocol.RewNeu=4;
ocol.ItemStim{1}=11;
ocol.ItemStim{2}=12;
ocol.ItemStim{3}=13;
ocol.Context_onset=32;
ocol.Context_end=40;
col.SimDis=1; % New spex (re-read)
col.RewNeu=2;
col.SimValcell=15;
col.Context_memscore=3;
col.Context_onset=4;
col.Context_duration=5;
col.ItemStim{1}=6;
col.ItemStim{2}=7;
col.ItemStim{3}=8;
col.ItemRoc{1}=9;
col.ItemRoc{2}=10;
col.ItemRoc{3}=11;
col.ItemMem{1}=12;
col.ItemMem{2}=13;
col.ItemMem{3}=14;
end

% Correct timings
f=spm_select('List', [where.databrain filesep subj '1 Preprocessed' filesep 'Func_b1'], '^fMQ');
if isempty(f)==1; nscans=250; disp('Fake block 1 n scans used'); else; nscans=size(f,1);end
e1(:, [30 31 32 33 35 37 39 40])=e1(:, [30 31 32 33 35 37 39 40])-start{1}/1000;
e2(:, [30 31 32 33 35 37 39 40])=e2(:, [30 31 32 33 35 37 39 40])-start{2}/1000+nscans*3;
e=vertcat(e1,e2);

% Write
d=nan*zeros(size(e,1),1);
d(:,col.SimDis)=e(:,ocol.SimDis);
d(:,col.RewNeu)=e(:,ocol.RewNeu);
d(:,col.Context_onset)=e(:,ocol.Context_onset);
d(:,col.Context_duration)=e(:,ocol.Context_end)-e(:,ocol.Context_onset);
d(d(:,col.SimDis)==1 & d(:,col.RewNeu)==1,col.SimValcell)=1; % Cell
d(d(:,col.SimDis)==1 & d(:,col.RewNeu)==0,col.SimValcell)=2;
d(d(:,col.SimDis)==2 & d(:,col.RewNeu)==1,col.SimValcell)=3;
d(d(:,col.SimDis)==2 & d(:,col.RewNeu)==0,col.SimValcell)=4;
for i=1:3 % Mark memory
    d(:,col.ItemStim{i})=e(:,ocol.ItemStim{i});
    
    for t=1:size(d,1) % Write Roc status of the item
        d(t, col.ItemRoc{i})=item(find(item(:,2)==d(t,col.ItemStim{i})), 6);
    end
end
switch memtype
    case 'Roc'
        d(:,col.Context_memscore)=d(:,col.ItemRoc{1})+d(:,col.ItemRoc{2})+d(:,col.ItemRoc{3});
    case 'Hit'
        d(:,col.Context_memscore)=(d(:,col.ItemRoc{1})>2.5)+(d(:,col.ItemRoc{2})>2.5)+(d(:,col.ItemRoc{3})>2.5);
end

% Vectors
% c_onset=d(:, col.Context_onset);
% c_dur=d(:, col.Context_duration);
% c_mempmod=d(:, col.Context_memscore);
for c=1:4
    dd=d(d(:,col.SimValcell)==c,:);
    c_onset{c}=dd(:, col.Context_onset);
    c_dur{c}=dd(:, col.Context_duration);
    c_mempmod{c}=dd(:, col.Context_memscore);
end

%% Load onsets for comparison

v=load([where.databrain filesep subj filesep '2 First level' filesep subj '_onsets_' model '_' memtype '.mat']);
error('DONE - from here, manual comparison')

% Manual compare - evalute following, changing r: 
r=4; 
disp('memscore pmods')
for r=1:4
    p=v.pmod(r).param{1}; % context scores
m=c_mempmod{r};
disp(sum(abs(p-m))) % this should ==0
end

disp('Onsets')
r=4; 
for r=1:4
    o=v.onsets{r}; % onsets
o1=c_onset{r};
disp(sum(abs(o-o1))) % this should ==0
end

disp('Durations')
r=3; 
for r=1:4
    d1=v.durations{r}; % duration
d2=c_dur{r};
disp(sum(abs(d1-d2))) % this should ==0
end

%% Load SPM.mat file, to compare

s=load([where.databrain filesep subj filesep '2 First level' filesep model '_' memtype ' Contrasted' filesep 'SPM.mat']); s=s.SPM;

r=1; disp('Check onsets')
for r=1:4
o=s.Sess.U(r).ons;
o1=c_onset{r};
disp(sum(abs(o-o1)))
end

r=1; disp('Check durations')
for r=1:4
o=s.Sess.U(r).dur;
o1=c_dur{r};
disp(sum(abs(o-o1)))
end

r=1; disp('Check mem pmods')
for r=1:4
p=s.Sess.U(r).P.P;
p1=c_mempmod{r};
disp(sum(abs(p-p1)))
end





