function [pts d] = HalfBellShape(Minx, Maxx, n, sigma2)

tx = linspace(-5,5,2*n);

d = exp(-(tx.*tx)/(2*sigma2))/(sqrt(2*pi*sigma2));

d = d((n+1):end);
pts = Minx + cumsum(d);

pts = Minx + (Maxx - Minx)*cumtrapz(d)/max(cumtrapz(d));