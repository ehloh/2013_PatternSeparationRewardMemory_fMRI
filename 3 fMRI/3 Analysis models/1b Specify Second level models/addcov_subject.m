function [covbatch] = addcov_subject(subjectlist, ncells, memeffecttype)
% [covbatch] = addcov_subject(subjectlist, ncells, memeffecttype)
% Include subjects as covariates
%
% ------------------------------------------------------------------------------
for s=1:length(subjectlist)
    sub=zeros(1,length(subjectlist)); sub(s)=1;
   covbatch.cov(s).c = repmat(sub,[1 ncells])';
   covbatch.cov(s).cname = ['sub_' subjectlist{s}];
   covbatch.cov(s).iCFI = 1;
   covbatch.cov(s).iCC = 1;
end

end

