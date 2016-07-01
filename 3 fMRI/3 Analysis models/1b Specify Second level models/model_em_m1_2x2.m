function [ matlabbatch contrastfiles secondlevelfolder] = model_em_m1_2x2(where_brain, secondlevelfolder, subjectlist, analysistype, firstlevelmodel, covariates, memtype,memdata)
% [ matlabbatch contrastfiles secondlevelfolder] = model_em_m1_2x2(where_brain, secondlevelfolder, subjectlist, analysistype, firstlevelmodel, covariates, memtype,memdata)
% Second level model, for basic 2x2 Factorial (Sim x Val)
%
% Options for memtype:      'Hit', 'Surehit', 'Rem', 'Surerem'
%
% Note: This model is for event-memory parametric modulators, i.e. memory
% pmods for items and memory pmods for contexts are combined at 1st level
% contrasts, to form first-level event-memory contrasts. These event-memory
% contrasts are fed into 2nd level. 
% -------------------------------------------------------------------------------

% Execute to debug:   where_brain=where.data_brain; secondlevelfolder=log.secondlevelfolder; subjectlist=log.subjects;  analysistype=log.AnalysisType; firstlevelmodel=log.firstlevelmodel; secondlevelmodel=log.secondlevelmodel;  memtype=w.memtype; covariates=log.covariates;

%% (1) Details for this model: 

% (a) Factorial design (Sim x Val) - Details of the 'factorialcells'  
%           Col 1=Name of cell's contrast, Col 2=Cell's status in each factor, 
%           Col 3= Relative contrast number (ignore), Col 4= Contrast file name
factorialcells={'SimR' [1 1];
                        'SimN' [1 2];
                        'DisR' [2 1];
                        'DisN' [2 2];};
disp(' '); disp('Factorial cells for requested 2nd-level model:'); for i=1:size(factorialcells,1); disp(['            ' factorialcells{i,1} '   '  num2str(factorialcells{i,2})]); end; input('OK?   ')

%% (2) Set up Covariates (if requested)

if ~isempty(covariates)

    ncells=4; % 2x2 Factorial
    eval(['[cov] = addcov_' covariates '(subjectlist, ncells, memdata);'])
    matlabbatch{1}.spm.stats.factorial_design=cov.covbatch;
    
    % Name + Organization of model
    secondlevelfolder =[secondlevelfolder ' ' cov.name]; mkdir(secondlevelfolder)
    
else
    
    % No covariates
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    mkdir(secondlevelfolder)
end
                          
%% (3) Specify 2nd level model (Full factorial)

matlabbatch{1}.spm.stats.factorial_design.dir = {secondlevelfolder};
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Similarity';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'Valence';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1; 1;1];
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Specify cells to Factorial Design
for i=1:size(factorialcells,1)
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).levels=factorialcells{i,2};
end

% Identify 1st-level contrast files for each subjects    
disp('Specifying contrast files for 2nd level ############# ')
contrastfiles=cell(length(subjectlist),1);
for s=1:length(subjectlist) 
    disp(['Subject ' num2str(s) '  -  ' subjectlist{s}])
    ws.where1st=[where_brain filesep subjectlist{s} filesep '2 First level' analysistype filesep firstlevelmodel ' Contrasted' filesep];
    ws.s=load([ws.where1st 'SPM.mat']); ws.SPM=ws.s.SPM;% Load contrasted SPM file
    
    % Identify & Asisgn the correct contrasts for this cell
    ws.cellcontrastnums=factorialcells; contrastfiles{s}=cell(4,1);
    for i=1:length(factorialcells)
        disp(['Identifying contrast no. ' num2str(i) '  -  ' factorialcells{i,1}])
        if strcmp(memtype,'Roc')==1
            [ws.contrastnames] = f_findcontrastname(ws.SPM, {[factorialcells{i,1} 'xEMem_' memtype] });
        else
            [ws.contrastnames] = f_findcontrastname(ws.SPM, {[factorialcells{i,1} 'xEMem_' memtype{1}] });
        end
        if isempty(ws.contrastnames)==1
            error(['Error: Could not find the correct contrast file for specified cell   -- ' subjectlist{s} '  cell: ' factorialcells{i,1}])
        end
        contrastfiles{s}{i}= [ws.where1st ws.contrastnames{1} ,',1']; % Assign
    end
    %
    ws=[];
end
disp('Specifying contrast files for 2nd level: DONE ############# ')

% Compile contrasts, between-subjects
for i=1:4
    for s=1:length(subjectlist)
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans{s,1}=contrastfiles{s}{i};
    end
end


end

