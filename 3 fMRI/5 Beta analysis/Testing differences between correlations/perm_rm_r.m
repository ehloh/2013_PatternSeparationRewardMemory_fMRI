function p = perm_rm_r(X,Y,iter,plot)
% function p = perm_rm_r(X,Y,iter,plot)
% Takes predictors (m by 2 matrix) (e.g. 2 scores per subject). 
% These values are shuffled iter number of times (default: 10000). 
% Y is the dependent variable, i.e. the one both columns in X are trying to
% predict.
% The distribution for the difference in correlation with Y is plot and calculated how
% extreme this value is, which is characterized by p-value. Plot indicates whether or not you
% want histogram. NaN observations are thrown out with their paired observation

% rng('shuffle') %seed number generator at unpredictable location

if size(X,2) ~= 2
    error('only takes 2 predictors per subject/row')
end
if size(X,1) ~= size(Y,1)
    error('Predictors X and data Y need to have same number of rows (observations)')
end
if nargin == 2
    iter = 10000;
    plot = true;
elseif nargin == 3
    plot = true;
elseif nargin == 4
    if isempty(iter)
        iter = 10000;
    end
end
   
%set any row that has a NaN to 1
whereNans = isnan(X(:,1)) | isnan(X(:,2)) | isnan(Y);
X = X(~whereNans,:);
Y = Y(~whereNans,:);
% calculate correlation between predictor and dependent variable.
r = corrcoef([X Y]);
% calculate difference in correlation
test_score = r(1,3)-r(2,3);
clear r

permuted_scores = nan(iter,1);
for i = 1:iter
    swap = randn(size(X,1),1)>0; %logical indicating which to swap
    permuted_data(swap,1:2) = [X(swap,2) X(swap,1)];
    permuted_data(~swap,1:2) = [X(~swap,1) X(~swap,2)];
    r = corrcoef([permuted_data Y]);
    % calculate difference in correlation
    permuted_scores(i,1) = r(1,3)-r(2,3);
    clear r
end

if plot
    figure
    hist(permuted_scores,100);
    hold on
    lims = get(gca,'YLim');
    line([test_score test_score],[lims(1) lims(2)]);
    hold off
end
%test how unlikely the difference in r-score is
p = sum(abs(permuted_scores)>=abs(test_score))/iter;
end
