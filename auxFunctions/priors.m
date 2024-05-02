function [logpriors, p] = priors(param, ylag, nolong) 

% Flags
if nargin < 3
      
    % One-component filter
    nolong = 0 ;
    
    if nargin < 2 
        
        % Add lags of y as predictors
        ylag = 0 ;
        
    end
    
end

au = 1 - 1e-4 ;

% 1) cycle roots
ar2roots = param(1 : 2) ;

p1 = normpdf(ar2roots(1), au, .2) ;

p2 = normpdf(sum(ar2roots), au, .2) * (sum(ar2roots) <= 1 ) ; 

% 2) Stationary variance
p3 = normpdf(param(3), .9, .2) * (param(3) < .95) ;

% 3) Stationary asymmetry
p4 = normpdf(param(4), .9, .2) * (param(4) < .95) ;

% 4) Student-t dof
p5 =  gampdf(1/param(end), 2, 10) * (param(end) < .34) ;

% 5) Score loading priors
if nolong
    
    scloadL = [] ;
    scloadS = param(5 : 7) ;
    
else

    scloadL = param(5 : 2 : 10) ;
    scloadS = param(6 : 2 : 10) ;
    
end

a1 = 4 ;
b1 = .3 ;
b2 = .1 ;

p6 = prod(gampdf(1./scloadS, a1, 1/b1)./(scloadS.^2))*prod(gampdf(1./scloadL, a1, 1/b2)./(scloadL.^2)) ;

% 6) Beta ridge priors (independent)
betas = param(11 : end - 5) ;

lsd = .01 ;
ssd = .2 ;

if ylag
    
    lsd = .01 ;
    
end

p7 = prod([normpdf(betas(1 : 3 : end), 0, lsd); normpdf(betas(2 : 3 : end), 0, ssd); normpdf(betas(3 : 3 : end), 0, ssd)]) ;

if isempty(p7); p7 = 1; end

% 7) Hessian smoothing parameter
p8 = unifpdf(param(end - 1), 0, .5) ;

% 8) Initial values/constant
p9 = prod(normpdf(param(end - 4 : end - 2), [2.5 1.5 -.05]', .5)) ; % var should be .5

% Joint prior distribution
logpriors = log(p1 * p2 * p3 * p4 * p5 * p6 * p7 * p8 * p9) ;

p = [p1 p2 p3 p4 p5 p6 p7] ;

end
