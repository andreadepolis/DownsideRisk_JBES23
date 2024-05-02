function stackedArea(xax, Y, clr, alp, hv)
% Plot stacked areas.
% Input:
%   xax:   (Tx1) vector of x-axis; can be numeric, datenum or datetime.
%   Y:     (Txn) matrix of bands (e.g., percentiles). If n is odd, the
%           central value is considere to be the median.
%   clr:   bands color; (nx3) RGB triplet.
%   alp:   bands color transparency. Can be specified as a (nx1) vector of
%           transparency levels, or as a single value. In the latter case,
%           the color gradient is achieved by overlaying bands with the
%           same transparency.
% 
% Andrea De Polis, 2022. a.de-polis@warwick.ac.uk

if nargin < 5

    hv = 'on' ;

    if nargin < 4

        alp = 1 ;

    end

end

n = size(Y, 2) ;

% Positive components. 
fh(1) = figure(get(gcf, 'Number')) ; 
hP    = area(xax, Y.*(Y>0), 'Facecolor', 'flat', 'Facealpha', 1, 'linestyle', 'none', 'EdgeColor', 'none', 'linewidth', .01) ; 
ahP   = gca ; 

% Negative components. 
fh(2) = figure ; 
hN    = area(xax, Y.*(Y<0), 'Facecolor', 'flat', 'Facealpha', 1, 'linestyle', 'none', 'EdgeColor', 'none', 'linewidth', .01) ; 

alp = repmat(alp, 1, 2*n) ;

for ii = 1 : n

    hP(ii).FaceColor = clr(ii, :) ;
    hP(ii).FaceAlpha = alp(ii) ;
    hP(ii).HandleVisibility = hv ;

    hN(ii).FaceColor = clr(ii, :) ;
    hN(ii).FaceAlpha = alp(ii) ;
    hN(ii).HandleVisibility = hv ;
    
end
   
set([hP hN], 'parent', ahP) ; 
close(get(gcf, 'Number'))

end