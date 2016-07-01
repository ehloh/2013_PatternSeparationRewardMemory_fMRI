



%% Load data for a manual check
where.data_beh=['H:\1 fMRI analysis' filesep '2 Behavioural'];
load([where.data_beh filesep 'p02_MK' filesep 'p02_MK_file_2encodingfMRI_b2.mat'])

mm=memtest.trialstats;
m(:,1)=mm(:,2);
m(:,3)=mm(:,12);
m(:,4)=mm(:,14);
m(:,5)=mm(:,16);
%
a=subjdata{1,3};
SimR=a.onsets{1}; %
SR_Hit=a.onsets{5};
SR_Miss=a.onsets{6};
SimN=a.onsets{2};  %
SN_Hit=a.onsets{7}; 
SN_Miss=a.onsets{8}; 
DisR=a.onsets{3}; %
DR_Hit=a.onsets{9};
DR_Miss=a.onsets{10};
DisN=a.onsets{4};  %
DN_Hit=a.onsets{11}; 
DN_Miss=a.onsets{12};

% MANUALLY INPUT THESE!
start=114155/1000;
start2=63042/1000;
offset=762;
% if keypresses, k must be divided by 1000. Otherwise, don't. 


%% Paste these items during hard check

% BLOCK 1
k=1;
while k~=999
    disp(['Item onset: ' num2str(-start+k)]); disp('------------'); k=input('Time for next item:  ');
end

% BLOCK 2
k=1;
while k~=999
    disp(['Item onset: ' num2str(-start2+offset+k)]); disp('------------'); k=input('Time for next item:  ');
end

% Check dat
subjdata{1,2}.items=sortrows(subjdata{1,2}.items,[-1 2 3 17 4]);
b=subjdata{1,2}.items;
a=subjdata{1,3}; a.names=a.names';


