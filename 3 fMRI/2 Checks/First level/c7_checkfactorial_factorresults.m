function [ details] = c7_checkfactorial_factorresults(spm, ContextOrItem, memtype)
% [ details] = c7_checkfactorial_factorresults(spm, ContextOrItem, memtype)
% Find contrast weights corresponding to factor contrasts
%
%       spm:                    SPM.mat variable for model (after 1st level contrasts)
%       ContextorItem:     Are we assessing Context or Item regressors?
%       memtype:            Hit/Surehit/Rem/Surerem  - Roc not written yet
%
%       details:                 Cell 
%                                       Col 1: Name of Factor/comparison
%                                       Col 2: Regressor weights corresponding to this comparison
%                                       Col 3: Contrast number
%
% ----------------------------------------------------------------------------------

% Execute: spm=wm.spm; ContextOrItem=2; memtype=log.memtype;

for o1=1:1 % Set up 
    
    % Which cells to combine for factors?
    %   Col 1: Name of factor/comparison
    %   Col 2: Cell, containing Cell names & their corresponding weights
    %               Col 1: Cell name (must correspond to contrast name)
    %               Col 2: Weights for this cell (these weights are then applied to
    %                             the regressor-weights that correspond to this
    which_context={'Visual'         {'SimR' 1; 'SimN' 1; 'DisR' 1; 'DisN' 1};
                            'Similarity'    {'SimR' 1; 'SimN' 1; 'DisR'  -1; 'DisN' -1};
                            'Valence'       {'SimR'  1; 'DisR' 1; 'SimN'  -1; 'DisN' -1};
                            };
    
    which_item={'Visual'            {'SR_m1' 1; 'SR_m0' 1; 'SN_m1' 1; 'SN_m0' 1;    'DR_m1' 1; 'DR_m0' 1; 'DN_m1' 1; 'DN_m0' 1};
                        'Similarity'        {'SR_m1' 1; 'SR_m0' 1; 'SN_m1' 1; 'SN_m0' 1;    'DR_m1' -1; 'DR_m0' -1; 'DN_m1' -1; 'DN_m0' -1};
                        'Valence'           {'SR_m1' 1; 'SR_m0' 1; 'SN_m1' -1; 'SN_m0' -1;    'DR_m1' 1; 'DR_m0' 1; 'DN_m1' -1; 'DN_m0' -1};
                        'Memory'            {'SR_m1' 1; 'SR_m0' -1; 'SN_m1' 1; 'SN_m0' -1;    'DR_m1' 1; 'DR_m0' -1; 'DN_m1' 1; 'DN_m0' -1};
                        };
    
    switch ContextOrItem
        case 1 % Context
            which=which_context;
        case 2 % Item
            which=which_item;
        otherwise
            error('Context=1, Item=2')
    end
    
    % If item regressors are being assessed, change condition names for memory
    if ContextOrItem==2
        switch memtype
            case 'Hit'
                memyes='Hit'; memno='Miss';
            case 'Surehit'
                memyes='Surehit'; memno='NotSurehit';
            case 'Rem'
                memyes='Rem'; memno='NotRem';
            case 'Surerem'
                memyes='Surerem'; memno='NotSurerem';
            otherwise; error(['Could not find requested memory type  -- ' memtype])
        end
        for w=1:size(which,1) % Re-naming
            for c=1:size(which{w,2},1)
                if strcmp(which{w,2}{c,1}(length(which{w,2}{c,1})),'1')==1
                    which{w,2}{c,1}=[which{w,2}{c,1}(1:length(which{w,2}{c,1})-2) memyes];
                elseif strcmp(which{w,2}{c,1}(length(which{w,2}{c,1})),'0')==1
                    which{w,2}{c,1}=[which{w,2}{c,1}(1:length(which{w,2}{c,1})-2) memno];
                end
            end
        end
    end
    
end

%% What are we working with?

% Details of all available contrasts
connames=cell(size(spm.xCon,2),1);weights=cell(size(spm.xCon,2),1);
for i=1:size(spm.xCon,2)
    connames{i}=spm.xCon(i).name;
    weights{i}=spm.xCon(i).c;
end

%% Which regressors to weight, for this factor?

details=cell(size(which,1),2);
for w=1:size(which,1)
    details{w,1}=which{w,1};
    regweights=nan*zeros(size(which{w,2},1), size(weights{1},1));
    
    % Identify conditions/regressors to be included in this comparison
    for c=1:size(which{w,2},1)
        connum=find(strcmp(connames, which{w,2}{c,1})==1);
        if length(connum)~=1; error('Error: No. of cells to be combined (in factorial check) ~=1. Requested contrast may not exist'); end
%         disp([which{w,2}{c,1}  '   --   ' connames{connum}])
        regweights(c,:)=(weights{connum}').*which{w,2}{c,2};
    end
    
    details{w,2}=sum(regweights,1);
    if sum(abs(details{w,2})>1)~=0
        input('Warning: In assigning condition-cells to factors, factorial-cell may be uneven or doubly assigned. Hit Enter to ignore');
    elseif sum(abs(details{w,2}))==0 
        error('Error: Some factorial cells are empty (no conditions/regressors assigned)')
    end
end

end

