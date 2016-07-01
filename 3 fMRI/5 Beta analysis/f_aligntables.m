function [ table headers  subs n_subs excl] = f_aligntables(table1, table2)
% [ table headers  subs n_subs excl] = f_aligntables(table1 table2)
% Align tables 1 & 2 by subject (i.e. so each subject's data 
%   from table 1 and 2 are aligned). 1st column assumed subjects, 1st row
%   assumed header titles. Subjects excluded from table 1 and table 2 are
%   in variable 'excl'
%
% ---------------------------------------------------------------------------

% Fake inputs:  clear all; clc;  w1=importdata('/Users/EleanorL/Dropbox/WorkPC/Explore beta/c3 choice betas/c3 ROI battery/(03-Nov-2013) Extracted betas.txt'); table1=[w1.textdata(:,1) [w1.textdata(1,2:end); num2cell(w1.data)]]; [w2.num w2.txt table2]=xlsread('/Users/EleanorL/Dropbox/SANDISK/5 Explore fMRI/1 Behavioural data/Group behaviour files/Behavioural profile for correlation.xlsx'); table2(24,:)=table2(2,:); table2{2,1}='bum'; table2(2,:)=[];

% Headers
head1=table1(1,2:end)'; head2=table2(1,2:end)';
headers=[head1;  head2];

% Subject list - remove blank spaces
table1(:,1)=cellfun(@deblank, table1(:,1), 'UniformOutput',0);
table2(:,1)=cellfun(@deblank, table2(:,1), 'UniformOutput',0);

% Compile aligned table (data only)
k=1;table=cell(0,0); excl.table1={}; excl.table2={}; subs={};
for s=2:size(table1,1)
    if length(find(strcmp(table2(:,1),table1{s,1})))==1
        subs=[subs; table1(s,1)];
        table(k,1:length(head1))=table1(s,2:end);
        table(k,length(head1)+1:length(head1)+length(head2))=table2(find(strcmp(table2(:,1),table1{s,1})),2:end);
        k=k+1;
    else
        excl.table1=[excl.table1; table1(s,1)];
    end
end
n_subs=length(subs);

% Exclusion list for table 2
for s=2:size(table1,1)
    if length(find(strcmp(table1(:,1),table2{s,1})))~=1
        excl.table2=[excl.table2; table2(s,1)];
    end
end

% Combine aligned table with subject list and headers
table=[[{'Subjects'}; subs] [headers';  table]];


end