function y = scaledInvX2pdf(x, tau2, v)
% scaled inverse chi-square distribution
c = ((tau2*v/2)^(v/2))/gamma(v/2);
ep = exp((-v*tau2)./(2*x));
dn = x.^(1+v/2);
y = c*ep./dn;

end