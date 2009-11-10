function [pts d] = BellShape(Minx, Maxx, n, sigma2)

tx = linspace(-5,5,n);

d = exp(-(tx.*tx)/(2*sigma2))/(sqrt(2*pi*sigma2));
pts = Minx + cumsum(d);

pts = Minx + (Maxx - Minx)*cumtrapz(d)/max(cumtrapz(d));