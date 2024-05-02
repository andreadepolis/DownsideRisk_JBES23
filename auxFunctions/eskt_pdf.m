function f = eskt_pdf(x, mu, sigma, epsilon, eta)
% Probability density function (PDF) of the epsilon-Skew-t of Gomez, Bolfarine and Torres (2007).
% 
% Input:  x      : support of the e-Skew-t;
%         mu     : location of the e-Skew-t; 
%         sigma  : scale of the e-Skew-t;
%         epsilon: shape of the e-Skew-t;
%         eta    : inverse of the degrees of freedom, nu = 1/eta.
% 
% Output: f: PDF of the e-Skew-t.
% 
% Andrea De Polis (a.de-polis@warwick.ac.uk), 2022.
x = x - mu ;

c = exp(gammaln((eta + 1)./(2*eta)) - gammaln(1./(2*eta)))./sqrt(pi./eta) ;

f = c./sigma.*(1 + eta*((x)./(sigma.*(1 - sign(x).*epsilon))).^2).^(-(1 + eta)/(2*eta)) ;

