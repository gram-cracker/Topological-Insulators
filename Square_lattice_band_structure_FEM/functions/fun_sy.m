function sy = fun_sy( y, y0, ly, k0 )
%FUN_SY returns sy value in the PML layer


sigma0 = 12;
pow = 2;

sigma_y = 0;

if (y<y0)
    sigma_y = sigma0 * ((y0-y)/ly)^pow;
end

sy = 1 + 1i*sigma_y/k0;
% sy = 1 + 1i*sigma_y;