% Check the inputs to the 2nd level model (ingoing contrast, cell
% allocations)
clear all; close all hidden; clc

subj='p01_CW';
nsubs=26;
memtype='Hit'; % 'Roc' 'Hit' 'Rem'
model='m_ci5_ContextItempmod';
% model='m_ci3_ContextItemNomotor';
secondlevel='cm_m1_2x2';

where.databrain='C:\Users\eloh\Desktop\1 [Context-Memory] fMRI Data';


%% Load 2nd level SPM.mat file, to check cell allocations (Assume factorial design)

m=load([where.databrain filesep '2 Second level analysis' filesep 'results     ' model '_' memtype '       ' secondlevel filesep 'SPM.mat']); m=m.SPM;

% Ingoing contrasts - Manual check here, within-cell contrast no.s
% consistent across subjects? Take contrast #, go to check 1st level
c_images=m.xY.P;
c_facAlloc=m.xX.X;
c_cellnames=m.xX.name; disp('Ingoing contrasts -  cell names: ');disp(c_cellnames)
% possibly also relevant: m.xVi.I - cell allocations? m.xY= ingoing
% scans/contrast images?

% Names of in-going contrasts (what exactly is each cell? 1,1=SR)
disp(' '); disp('What are the Names of the contrasts going into each cell?')
for i=1:4
    disp(['Cell ' num2str(i) '---'])
    disp(m.xY.VY(1+nsubs*(i-1)).descrip )
end

%% Load a subject's 1st level
%   check that the contrast images going into the 2nd level have correct contrast weights for regrssors

s=load([where.databrain filesep '1 MRI data' filesep subj filesep '2 First level' filesep model '_' memtype ' Contrasted' filesep 'SPM.mat']); s=s.SPM;
s_regnames=s.xX.name';
s_allcons=cell(size(s.xCon,2),1);
for i=1:size(s.xCon,2); s_allcons{i}=s.xCon(i).name; end
disp('List of contrasts from the single subject''s specified model:')
for i=1:size(s.xCon,2); disp(['con' num2str(i) '    ' s_allcons{i}]); end

% Read from above, which contrast no. to check?
connum=8;
disp('################################################')
disp(['Requested contrast num: ' num2str(connum)])
disp('(1) Name of this contrast (in subject''s specified model):  ')
disp(s.xCon(connum).name)
disp('(2) Regressors (in subject''s specified model) that are weighted in this contrast:')
disp(find(s.xCon(connum).c))
disp('(3) The name of this weighted regressor (in the subject''s specified model) is:')
disp(s_regnames(find(s.xCon(connum).c)))
disp('--------------------------------------------------------------------------------------------')




