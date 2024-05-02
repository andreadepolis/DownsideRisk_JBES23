function [ecdf] = pdf2ecdf(pdf, support)

p    = pdf./trapz(support, pdf) ;
pxi  = linspace(min(support), max(support), size(support, 1))' ;
pi   = interp1(support, p, pxi, 'linear') ;
ecdf = cumtrapz(pxi, pi) ;
