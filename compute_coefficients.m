function [a_i, b_i,k,w0] = compute_coefficients(t_p)
syms t
f0 = 0.1/t_p; p0 = 1/f0; w0 = 2*pi*f0;  I = [-p0/2; p0/2];
a = I(1);  b = I(2); P = 1; k = 10000; i = linspace(1,k,k);
f = @(t) P.*exp(-2*log(2).*((t-2*t_p)./t_p).^2);
ai = @(t) (2./p0).*(P.*exp(-2*log(2).*((t-2*t_p)./t_p).^2).*cos(2.*pi.*i.*t./p0));
bi = @(t) (2./p0).*(P.*exp(-2*log(2).*((t-2*t_p)./t_p).^2).*sin(2.*pi.*i.*t./p0));
a0 = @(t) (1./p0).*(P.*exp(-2*log(2).*((t-2*t_p)./t_p).^2));
a_0 = integral(a0,0,b);
a_i = integral(ai,a,b,'ArrayValued',true);
b_i = integral(bi,a,b,'ArrayValued',true);