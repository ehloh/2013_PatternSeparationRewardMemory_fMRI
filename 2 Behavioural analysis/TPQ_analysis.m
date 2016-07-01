% Extract TPQ scores
%   Scores are saved to subject's file, and printed to txt file

clear all
clc

load i_Blueprint;
scriptloc='H:\2 fMRI behaviour analysis\';
dataloc='C:\Users\eloh\Desktop\0 [Context-Memory] fMRI Data\Behavioural';
% n_subjects=22; % How many subjects have data?
exptversion='2.5 fMRI';
cd(scriptloc)

for o1=1:1 % Set up txt files (headers) 
cd('3 Analysis inputs')
w.fid=fopen(['(' date ') 1.3 TPQ.txt'],'wt');
fprintf(w.fid,['ExptVersion' '\t' 'Subject' '\t'  ]);
fprintf(w.fid, ['p.NS_s1' '\t' 'p.NS_s2' '\t' 'p.NS_s3' '\t' 'p.NS_s4' '\t']);
fprintf(w.fid, ['p.RD_s1' '\t' 'p.RD_s2' '\t' 'p.RD_s3' '\t' 'p.RD_s4']);
fprintf(w.fid,'\n');
end

%%

cd(dataloc)
w.dir=dir;

   for i1=3:length(w.dir)
        subjname=w.dir(i1).name;
%         subjname='p120_KL'; 
        disp(subjname)
        cd(subjname)
%%

% Templates
bp = Blueprint(2:101,1);
Res.NS=[0 0 0 0]'; % Novelty Seeking Subscale with NS1, NS2, NS3 and NS4 scores in the respective columns
Res.HA=[0 0 0 0]'; % Harm Avoidance Subscale with HA1, HA2, HA3 and HA4 scores in the respective columns
Res.RD=[0 0 0 0]'; % Reward Dependance Subscale with RD1, RD2, RD3 and RD4 scores in the respective columns
try % Is TPQ data available?
    load ([subjname '_data_TPQ.mat'])    
    tpq = TPQdata(2:101, 4);
    % Counters    
    for ItemNumb=1:size(tpq,1); 
          Item = tpq(ItemNumb,1);
          bpItem = bp(ItemNumb,1);

          if Item == bpItem;
              add=1;
          else
              add=0;
          end

       if ItemNumb == 2 || ItemNumb == 4 || ItemNumb == 43 || ItemNumb == 9 || ItemNumb == 11 || ItemNumb == 40 || ItemNumb == 85 || ItemNumb == 93 || ItemNumb == 96;
              Res.NS(1) = Res.NS(1)+ add;
          elseif ItemNumb == 30 || ItemNumb == 48 || ItemNumb == 50 || ItemNumb == 46 || ItemNumb == 55 || ItemNumb == 56 || ItemNumb == 81 || ItemNumb == 99;
              Res.NS(2) = Res.NS(2)+ add;
          elseif ItemNumb == 70 || ItemNumb == 72 || ItemNumb == 32 || ItemNumb == 66 || ItemNumb == 76 || ItemNumb == 78 || ItemNumb == 87;
              Res.NS(3) = Res.NS(3)+ add;
          elseif ItemNumb == 13 || ItemNumb == 22 || ItemNumb == 24 || ItemNumb == 28 || ItemNumb == 60 || ItemNumb == 62 || ItemNumb ==16 || ItemNumb == 21 || ItemNumb == 35 || ItemNumb == 65;
              Res.NS(4) = Res.NS(4)+ add;

          elseif ItemNumb == 5 || ItemNumb == 10 || ItemNumb == 14 || ItemNumb == 1 || ItemNumb == 8 || ItemNumb == 82 || ItemNumb == 84 || ItemNumb == 91 || ItemNumb == 95 || ItemNumb == 98;
              Res.HA(1) = Res.HA(1)+ add;
          elseif ItemNumb == 18 || ItemNumb == 19 || ItemNumb == 23 || ItemNumb == 26 || ItemNumb == 29 || ItemNumb == 47 || ItemNumb == 51;
              Res.HA(2) = Res.HA(2)+ add;
          elseif ItemNumb == 33 || ItemNumb ==37 || ItemNumb ==38 || ItemNumb ==42 || ItemNumb ==44 || ItemNumb ==89 || ItemNumb ==100;
              Res.HA(3) = Res.HA(3)+ add;
          elseif ItemNumb == 49 || ItemNumb ==54 || ItemNumb ==57 || ItemNumb ==68 || ItemNumb ==69 || ItemNumb ==73 || ItemNumb ==59 || ItemNumb == 63 || ItemNumb == 75 || ItemNumb == 80;
              Res.HA(4) = Res.HA(4)+ add;

          elseif ItemNumb == 27 || ItemNumb ==31 || ItemNumb ==34 || ItemNumb ==83 || ItemNumb ==94;
              Res.RD(1) = Res.RD(1)+ add;
          elseif ItemNumb == 39 || ItemNumb ==41 || ItemNumb ==77 || ItemNumb ==92 || ItemNumb ==97 || ItemNumb ==45 || ItemNumb ==52 || ItemNumb ==53 || ItemNumb ==79;
              Res.RD(2) = Res.RD(2)+ add;
          elseif ItemNumb == 3 || ItemNumb ==6 || ItemNumb ==7 || ItemNumb ==64 || ItemNumb ==67 || ItemNumb ==74 || ItemNumb ==12 || ItemNumb ==15 || ItemNumb ==86 || ItemNumb ==88 || ItemNumb ==90;
              Res.RD(3) = Res.RD(3)+ add;
          elseif ItemNumb == 17 || ItemNumb ==20 || ItemNumb ==25 || ItemNumb ==36 || ItemNumb ==58;
              Res.RD(4) = Res.RD(4)+ add;

          else 
              Res.NS(1) = Res.NS(1);
        end
    end
catch % If TPQ data is not available
    Res.NS=[999 999 999 999]';
    Res.HA=[999 999 999 999]';
    Res.RD=[999 999 999 999]';
end
%%

    for o5=1:1 % Finish off for this subject
        % Write numbers to txt file
%         fprintf(w.fid,[exptversion '\t' subjname '\t']);
%         fprintf(w.fid,'%0.5g \t',[Res.NS(1) Res.NS(2) Res.NS(3) Res.NS(4)]);
%         fprintf(w.fid,'%0.5g \t',[Res.RD(1) Res.RD(2) Res.RD(3) Res.RD(4)]);
%         fprintf(w.fid,'\n');
        % Save to memorytest
        load ([subjname '_file_3memorytest.mat'])
        memtest.TPQ=Res;
        save([subjname '_file_3memorytest.mat'], 'memtest')
    end
%     save(['p' subjname '_TPQResults'], 'Res');    
    clear Res
    cd ..    
   end

    
%
cd(dataloc)
fclose(w.fid);