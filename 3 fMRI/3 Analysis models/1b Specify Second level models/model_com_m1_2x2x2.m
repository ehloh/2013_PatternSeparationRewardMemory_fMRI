function [ matlabbatch contrastfiles secondlevelfolder] = model_com_m1_2x2x2(where_brain, secondlevelfolder, subjectlist, analysistype, firstlevelmodel, covariates, memtype,memdata)
% [ matlabbatch contrastfiles secondlevelfolder] = model_iom_m1_2x2x2(where_brain, secondlevelfolder, subjectlist, analysistype, firstlevelmodel, covariates, memtype,memdata)
% Second level model, for basic 2x2x2 Factorial (Sim x Val x Mem)
%
% Options for memtype:      'Hit', 'Surehit', 'Rem', 'Surerem'
%
% -------------------------------------------------------------------------------

% Execute to debug: where_brain=where.data_brain; subjectlist=log.subjects; firstlevelmodel=log.firstlevelmodel; secondlevelmodel=log.secondlevelmodel;  memtype=w.memtype; covariates=log.covariates

%% (1) Details for this model: Factorial design (Sim x Val x Item Mem)

% Details of the 'factorialcells'  
%           Col 1=Name of cell's contrast, Col 2=Cell's status in each factor, 
%           Col 3= Relative contrast number (ignore), Col 4= Contrast file name
factorialcells={'SimR_L' [1 1 1 ];           % Contrasts names must match the names specified in the 1st level
                        'SimR_H' [1 1 2]; 
                        'SimN_L' [1 2 1];  
                        'SimN_H' [1 2 2]; 
                        'DisR_L' [2 1 1]; 
                        'DisR_H' [2 1 2];  
                        'DisN_L' [2 2 1]; 
                        'DisN_H' [2 2 2];}; 

% for i=1:size(factorialcells,1) % Specify names according to Type of Memory 
%     if strfind(factorialcells{i,1}, '_m1')
%         factorialcells{i,1}=[factorialcells{i,1}(1:length(factorialcells{i,1})-2) memtype{1}];
%     elseif strfind(factorialcells{i,1}, '_m0')
%         factorialcells{i,1}=[factorialcells{i,1}(1:length(factorialcells{i,1})-2) memtype{2}];
%     else
%         error('Could not find specified item memory statuses')
%     end
% end

mkdir(secondlevelfolder)
disp(' '); disp('Factorial cells for requested 2nd-level model:');  for i=1:size(factorialcells,1); disp(['            ' factorialcells{i,1} '   '  num2str(factorialcells{i,2})]); end; input('OK?   ')
             
%% (2) Set up Covariates (if requested)

if ~isempty(covariates)

    ncells=8; % 2x2 Factorial
    eval(['[cov] = addcov_' covariates '(subjectlist, ncells, memdata);'])
    matlabbatch{1}.spm.stats.factorial_design=cov.covbatch;
    
    % Name + Organization of model
    secondlevelfolder =[secondlevelfolder ' ' cov.name]; mkdir(secondlevelfolder)
    
else
    
    % No covariates
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
end

%% (3) Specify model for Factorial analysis

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
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).name = ['MemLevel' memtype{1}];
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).ancova = 0;
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
contrastfiles=cell(length(subjectlist),1);
for s=1:length(subjectlist) 
    ws.where1st=[where_brain filesep subjectlist{s} filesep '2 First level' analysistype filesep firstlevelmodel ' Contrasted' filesep];
    ws.s=load([ws.where1st 'SPM.mat']); ws.SPM=ws.s.SPM;% Load contrasted SPM file
    
    % Identify the correct contrasts for this cell
    ws.cellcontrastnums=factorialcells;
    for i=1:length(factorialcells)
        ws.f=0; j=1;
        while ws.f==0
            if strcmp(ws.SPM.xCon(j).name, factorialcells{i})==1
                ws.cellcontrastnums{i,3}=j;
                ws.cellcontrastnums{i,4}=ws.SPM.xCon(j).Vcon.fname;
                ws.f=1;
            elseif j==size(ws.SPM.xCon,2)
                error(['Error: Could not find the correct contrast file for specified cell   -- ' subjectlist{s} '  cell: ' factorialcells{i}])
            else
                j=j+1;
            end
        end
    end
    
    % Assign contrasts to the correct factorial cells
    contrastfiles{s}=cell(8,1);
    for i=1:8
        contrastfiles{s}{i}=[ws.where1st  ws.cellcontrastnums{i,4}];
    end
    
    %
    ws=[];
end

% Compile contrasts, between-subjects
for i=1:8
    for s=1:length(subjectlist)
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans{s,1}=contrastfiles{s}{i};
    end
end

end

