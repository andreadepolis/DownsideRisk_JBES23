function [LL, TVP, checkres, checkfrc, parms] = sktFilter(parms, y, x, h, fdraws, c2, sumll)

% INPUT:  parms  = model's static parameters ;
%           y    = vector of data ;
%           x    = matrix of exogenous regressors ;
%           h    = bootcast spteps
%           d    = direct forecast (4)
%         inival = TVP initial values ;
%            frc = in-sample(0) or out-of-sample(1)
%         fdraws = number of drwas
% 
% OUTPUT: LL       = model's loglikelihood ;
%         TVP      = series of TVP ;
%         checkres = selection amtrics
%         checkfrc = forecast details

if ~exist('c2', 'var')
    
    c2 = 1 ;
    
end

if ~exist('sumll', 'var')
    
    sumll = 1 ;
    
end

% Sizes
T = max(size(y)) ;
v = 7 ; %max(size(inival)) ;  
K = size(x, 2) ;

% Preallocations 
[sig21, sig22, crho1, crho2, ezf, ezl, vrf, vrl, mu, aux, s2] = deal(zeros(T + h + 1, 1)) ;
[N, J, I, JIJ, S] = deal(zeros(T + h, 3)) ;
ll                = zeros(T, 1) ;
TVP               = zeros(T + h + 1, v) ;
[ecf, den]        = deal(nan(fdraws, h + (h == 0))) ;
yF                = nan(7, h + (h == 0)) ;
yf                = nan(1, h + (h == 0)) ;
checkres          = struct ;
checkfrc          = struct ;

c  = .95 ;
a  = 50 ;
z  = linspace(-a, a, fdraws) ;
dz = 2*a/fdraws ;

% Selection matrices -> f_t = W + A*f_{t-1} + B*S_{t-1} + C*X_{t-1}

A = zeros(v, v) ; 

A(1, 1) = 1 ; %parms(4) ;
A(2, 2) = parms(1) ; 
A(2, 3) = parms(2) ;
A(3, 2) = 1 ;
A(4, 4) = 1 ; % parms(7) ;
A(5, 5) = parms(3) ;
A(6, 6) = 1 ; %parms(7) ;
A(7, 7) = parms(4) ;

B = zeros(v, 3) ;

if c2

    B(1, 1) = parms(5) ;  % k_beta
    B(2, 1) = parms(6) ;  % k_c
    B(4, 2) = parms(7) ;  % k_sig2
    B(5, 2) = parms(8) ;  % k_sig2
    B(6, 3) = parms(9) ;  % k_vrho
    B(7, 3) = parms(10) ; % k_vrho
        
else

    B(1, 1) = 0 ;%parms(5) ;  % k_beta
    B(2, 1) = parms(5) ;  % k_c
    B(4, 2) = 0 ;%parms(7) ;  % k_sig2
    B(5, 2) = parms(6) ;  % k_sig2
    B(6, 3) = 0 ;%parms(9) ;  % k_vrho
    B(7, 3) = parms(7) ; % k_vrho
    
end

X = zeros(T + h, 1) ;
C = zeros(v, 1) ;

if ~isempty(x)

    C   = zeros(v, K) ;
    
    C(2, :) = parms(11 : 3 : end - 4 - 3*c2) ;
    C(5, :) = parms(12 : 3 : end - 3 - 3*c2) ;
    C(7, :) = parms(13 : 3 : end - 2 - 3*c2) ;

    X = [x; x(end, :).*ones(h, 1)] ; 

end

smooth = parms(end - 1) ; 
eta    = parms(end) ;

vdf = 1/(1 - 2*eta) ; % Student-t 

% Loglikelihood constant
cost = log(gamma((eta + 1)/(2*eta))) - log(gamma(1/(2*eta))) - 0.5*(log(1/eta) + log(pi)) ;

% Inital values
TVP(1 : h + 1, 1) = parms(end - 4) ;%inival(1) ;
TVP(1 : h + 1, 2) = 0 ;
TVP(1 : h + 1, 3) = 0 ;
TVP(1 : h + 1, 4) = parms(end - 3) ;%inival(4) ;
TVP(1 : h + 1, 5) = 0 ;
TVP(1 : h + 1, 6) = parms(end - 2) ;%inival(6) ;
TVP(1 : h + 1, 7) = 0 ;

% Flag to control for bad likelihood values
regularization = 0 ;

j  = randi([h + 1, T - h]) ;
kk = 1 ;

for t = h + 1 : T + h
    
    mu(t) = TVP(t, 1) + TVP(t, 2) ;

    sig21(t) = exp(2*TVP(t, 4)) ;
    sig22(t) = exp(2*TVP(t, 5)) ;
    sig2     = sig21(t) * sig22(t) ;

    crho1(t) = tanh(TVP(t, 6)) ;
    aux(t)   = tanh(TVP(t, 7)) ;
    crho     = c*(crho1(t) + aux(t))/(1 + crho1(t)*aux(t)) ;
    crho2(t) = crho - crho1(t) ;
    
    if t <= T
    
        % Epsilon
        eps  = y(t) - mu(t) ;
        
        eps2  = eps^2 ;
        zeta2 = eps2/sig2 ; 

        % Likelihood
        ll(t) = - 0.5*log(sig2) - ((1 + eta)/(2*eta))*log(1 + (eta*eps2)/((1 - sign(eps)*crho)^2*sig2)) ;

        % Scaled score
        w = (1 + eta)/((1 - sign(eps)*crho)^2 + eta*zeta2) ;

        N(t, 1) = (w*eps)/sig2 ;
        N(t, 2) = (w*zeta2-1)/(2*sig2) ;
        N(t, 3) = -(sign(eps)/(1-sign(eps)*crho))*w*zeta2 ;

        J(t, 1) = 1 ;
        J(t, 2) = 2*sig2 ;
        J(t, 3) = (c^2 - crho^2)/c ;

        I(t, 1) = (1+eta)/((1+3*eta)*(1-crho^2)*sig2) ;
        I(t, 2) = 1/(2*(1+3*eta)*sig2^2) ;
        I(t, 3) = (3*(1+eta))/((1+3*eta)*(1-crho^2)) ;

        JIJ(t, :) = (J(t, :).*I(t, :).*J(t, :)).^.5 ;

        % Smoothing hessian
        if t ~= 1
    
            JIJ(t, :) = smooth * JIJ(t, :) + (1 - smooth) * JIJ(t - 1, :) ;
    
        end

        S(t, :) = JIJ(t, :).\J(t, :).*N(t, :) ;
        
    else
        
        den(:, kk) = epsilonskewt(z, mu(t), sqrt(sig2), crho, eta) ;
       
        if nansum(den(:, kk))*dz > .8
        
            ecf(:, kk) = pdf2ecdf(den(:, kk), z') ;

            [~, pos] = min(abs([.5 .05 .25 .32 .68 .75 .95] - ecf(:, kk))) ;

            yF(:, kk) = z(pos) ;
            yf(1, kk) = randsample(z, 1, true, den(:, kk)) ;

        else
            
            regularization = inf ;
            
        end
          
        % Bootcast (Koopman et al., 2019)
        S(t, :) = JIJ(j, :).\J(j, :).*N(j, :) ;

        kk = kk + 1 ;
         j =  j + 1 ;
        
    end

    % Updating
    TVP(t + 1, :) = A*TVP(t, :)' + B*S(t, :)' + C*X(t, :)' ;
    
    % Compute E[X] and V(X) hom components
    ezf(t) = (4*exp(cost)*crho*sqrt(sig2))/(1 - eta) ;
    ezl(t) = (4*exp(cost)*crho1(t)*sqrt(sig21(t)))/(1 - eta) ;
    vrf(t) = (3*crho^2)/(1 - 2*eta) - (ezf(t)/sqrt(sig2))^2 ;
    vrl(t) = (3*crho1(t)^2)/(1 - 2*eta) - (ezl(t)/sqrt(sig21(t)))^2 ;
    
    s2(t) = sig2 ;
    
end

llrec = - ll(12 : end)/T + regularization ;

if sumll
    
    LL = - cost - sum(ll(12 : end))/T + regularization ;
    
else
   
    LL = (-cost*(T/(T-11)) - ll(12 : end))/T + regularization ;

end
    
if (isnan(LL)); LL = inf ; end

% E[X] and V(X)

Exf = mu - ezf ;
Exl = TVP(:, 1) - ezl ;
Vxf = s2.*(vdf + vrf) ;
Vxl = sig21.*(vdf + vrl) ;

nu = 1/eta ;

checkres.seleA = A ;
checkres.seleB = B ;
checkres.seleC = C ;
checkres.mu    = mu ;
checkres.Expvf = ezf ;
checkres.Expvl = ezl ;
checkres.Varnf = vrf ;
checkres.Varnl = vrl ;
checkres.Vdof  = vdf ;
checkres.Exf   = Exf ;
checkres.Exl   = Exl ;
checkres.Vxf   = Vxf ;
checkres.Vxl   = Vxl ;
checkres.tvp   = TVP ;
checkres.scs   = S ;
checkres.gnu   = 4*nu*exp(cost)/(nu - 1) ;
checkres.hnu   = 3*nu/(nu - 2) -  checkres.gnu^2 ;
checkres.llrec = llrec ;
checkres.cost  = -cost ;

TVP(:, [4 5 6 7 8]) = [sig21 sig22 crho1 crho2 aux] ;

checkfrc.ECDF = ecf ; 
checkfrc.EPDF = den ;
checkfrc.fden =  yF ;
checkfrc.forc =  yf ;

end
