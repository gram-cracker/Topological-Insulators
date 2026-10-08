
y0 = 0;
nn = 101;
y = linspace(-5,5,nn);
ly=1;
k0 = 2;

sy = zeros(1,nn);

for n = 1:nn
    sy(n) = fun_sy( y(n), y0, ly, k0 );
end

plot(y,imag(sy));