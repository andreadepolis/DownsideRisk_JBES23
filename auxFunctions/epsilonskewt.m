function d = epsilonskewt(x, mu, sigma, epsilon, eta)

x = x - mu ;

c = gamma((eta + 1)/(2*eta))/(gamma(1/(2*eta))*sqrt(pi/eta)) ;

d = c./sigma.*(1 + eta*((x)./(sigma.*(1 - sign(x).*epsilon))).^2).^(-(1 + eta)/(2*eta)) ;

