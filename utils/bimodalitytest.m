% R. Pfister, K. Schwarz, M. Janczyk, R. Dale, J. Freeman. Good things peak in pairs: a note on the bimodality coefficient. Frontiers in Psychology, Vol. 4, Art. 700, Oct. 2013. 
function [isbimodal, coef] = bimodalitytest(xpdf)
n = numel(xpdf);

skew = skewness(xpdf, 0); % skew
kurt = kurtosis(x, 0) - 3; % kurtosis, directly copied from cited paper

coefnum = skew^2 + 1;
coefden_fracnum = (n-1)^2;
coefden_fracden = (n-2)*(n-3);
coefden = kurt + 3*coefden_fracnum/coefden_fracden;

coef = coefnum / coefden;

isbimodal = coef > 5/9;
end