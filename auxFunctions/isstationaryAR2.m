function res = isstationaryAR2(phi1, phi2)

X = [phi1 phi2; 1 0] ;

eigv = eig(X) ;
meig = max(abs(eigv)) ;

res = meig <= .989 ;
