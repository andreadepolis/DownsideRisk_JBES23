function plotbands(xax, bands, clr, alp)
% Plot shaded bands.
% Input:
%   xax:   (Tx1) vector of x-axis; can be numeric, datenum or datetime.
%   bands: (Txn) matrix of bands (e.g., percentiles). If n is odd, the
%           central value is considere to be the median.
%   clr:   bands color; can be color name (e.g., 'b') of RGB triplet.
%   alp:   bands color transparency. Can be specified as a (nx1) vector of
%           transparency levels, or as a single value. In the latter case,
%           the color gradient is achieved by overlaying bands with the
%           same transparency.
% 
% Andrea De Polis, 2022. a.de-polis@warwick.ac.uk

[r, c] = size(bands) ;

if r > c

    T = r ;
    n = c ;

    bands = bands' ;

else

    T = c ;
    n = r ;
    
end

hn = floor(n/2) ;

if nargin < 4

    alp = repmat(.2, hn, 1) ;

    if nargin < 3 

        clr = 'b' ;

    end

elseif length(alp) ~= hn

    alp = repmat(alp(1), hn, 1) ;

end

if size(xax, 1) == T

    xax = xax' ;

end

nanID = all(isfinite(xax), 1) & all(isfinite(bands), 1) ;

xax   = xax(:, nanID) ;
bands = bands(:, nanID) ;

if size(xax, 1) ~= hn

    xax = repmat(xax, n, 1) ;

end

isVisible = 'off' ;
med = nan(T, 1) ;
if mod(n - 1, 2) == 0

    med = bands(n - hn, :) ;
    bands(n - hn, :) = [] ;

    isVisible = 'on' ;

end

hold on
for ii = 1 : hn

    xx = [xax(ii, :) fliplr(xax(ii, :))] ;
    yy = [bands(ii, :) fliplr(bands(ii + hn, :))] ;

    fill(xx, yy, clr, 'facealpha', alp(ii), 'edgecolor', 'none', 'handlevisibility', 'off') ;
 
end
plot(xax, med, 'color', clr, 'handlevisibility', isVisible)
hold off