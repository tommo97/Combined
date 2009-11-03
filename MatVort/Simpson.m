function I = Simpson(f,h)
n=length(f)-1;
I=0;
switch n
    case 1
        disp('Too few accissas');
    case 2
        I = h/3*(f(1) + 4*f(2) + f(3));
    case 3
        I = 3/8*h*(f(1) + 3*f(2) + 3*f(3) + f(4));
    otherwise
        if  2*floor(n/2)~=n,
            I = 3/8*h*(f(n -2) + 3*f(n -1)  ...
                + 3*f(n) + f(n+1));
            m=n -3;
        else
            m=n;
        end
        I = I+ (h/3)*( f(1)+ 4*sum(f(2:2:m)) + f(m+1));
        if m>2
            I = I+ (h/3)*2*sum(f(3:2:m));
        end
end
