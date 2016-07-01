function[ h, p, SteigerZ ] = BenStuff_CorrCompSteigerZ( r12,r13, r23, n, alpha)
%[ h, p, SteigerZ ] = BenStuff_CorrCompSteigerZ( r12,r13, r23, n, alpha) 
%   Compares dependent (matched) r's for overlapping variables (two predictors same criterion)
%   h: hypothesis: 1 if difference is significant, else 0
%   p value of difference
%
%   --- note that h and p are based on z-distribution and rest on
%   sufficiently large n! ---
%
%   SteigerZ: z-value of difference
%
%   r12: corr between first predictor and criterion
%   r13: corr between second predictor and criterion
%   r23: corr between preictors
%   n= number of subjects (pairings for each predictor)
%   alpha: alpha level (defaults to .05)

%   Equations nicked from ...
%   http://jeromyanglim.blogspot.com/2009/09/significance-tests-on-correlations.html (23/2/2012)
%
%   keen on overkill? original publication: Steiger, 1980
%   
%   Fisher's Transformation is used
%   en.wikipedia.org/wiki/Fisher_transformation
%   http://www.jstor.org/stable/2331838 - Fisher, 1915 - to really do your
%   head in...
%   
%
%   found a bug? please let me know!
%   benjamindehaas@gmail.com 2/2012

if nargin<5
    alpha=.05;
end

rmSq=(r12^2+r13^2)./2;
f=(1-r23)/2*(1-rmSq);
h=(1-(f*rmSq))/(1-rmSq);

Z12=atanh(r12);%Fisher transformation
Z13=atanh(r13);

SteigerZ=(Z12-Z13)*(sqrt(n-3)/(sqrt(2*(1-r23)*h)));
CriticalZ=(-erfinv(alpha - 1)) .* sqrt(2); %e.g. 1.96 for alpha=.05 (n is not assumed to be tiny here!)

if abs(SteigerZ)>CriticalZ
    h=1;
else
    h=0;
end

p=normcdf(-abs(SteigerZ),0,1)*2;%two-tailed


end

