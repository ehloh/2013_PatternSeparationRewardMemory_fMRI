function[ h, p, Z ] = BenStuff_CorrComp( r1, r2, n1, n2, alpha)
%[ h, p, Z ] = BenStuff_CorrComp( r1, r2, n1, n2, alpha)
%   Compare two independent r's
%   h: hypothesis: 1 if difference is significant, else 0
%   Z: z-value of difference
%   p value of difference
%
%   r1: corr in first sample
%   r2: corr in second sample
%   n1= number of subjects for r1
%   n2= number of subjects for r2
%   alpha: alpha level (defaults to .05)

%   Equations nicked from ...
%   http://www.jessicagrahn.com/uploads/6/0/8/5/6085172/comparecorrcoeff.doc
%   cf. Fisher, 1921: http://digital.library.adelaide.edu.au/dspace/handle/2440/15169
%   
%   Fisher's Transformation is used
%   en.wikipedia.org/wiki/Fisher_transformation
%   http://www.jstor.org/stable/2331838 - Fisher, 1915 - to really do your
%   head in...
%   
%
%   found a bug? please let me know!
%   benjamindehaas@gmail.com 10/2012

if nargin<5
    alpha=.05;
end


Z1=atanh(r1);%Fisher transformation
Z2=atanh(r2);

Z=(Z1-Z2)/(sqrt((1/(n1-3))+(1/(n2-3))));
CriticalZ=(-erfinv(alpha - 1)) .* sqrt(2); %e.g. 1.96 for alpha=.05 (n is not assumed to be tiny here!)

if abs(Z)>CriticalZ
    h=1;
else
    h=0;
end

p=normcdf(-abs(Z),0,1)*2;%two-tailed


end

