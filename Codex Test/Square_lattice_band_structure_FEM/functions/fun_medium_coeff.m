function val = fun_medium_coeff(x,y,c,c0,r1,r2)
% medium coefficient at (x,y)

shiftx = 0.0; shifty=0.0;

c0_x = 0.5; c0_y = 0.5;
c1_x = 0+shiftx;   c1_y = 0+shifty;
c2_x = 1-shiftx;   c2_y = 0+shifty;
c3_x = 1-shiftx;   c3_y = 1-shifty;
c4_x = 0+shiftx;   c4_y = 1-shifty;

% flag0 = ((x-0.5).^2 + (y-0.5).^2 - r1^2) <= 0;
% flag1 = (x.^2 + y.^2 - r2^2) <= 0;
% flag2 = ((x-1).^2 + y.^2 - r2^2) <=0 ;
% flag3 = ((x-1).^2 + (y-1).^2 - r2^2) <= 0;
% flag4 = (x.^2 + (y-1).^2 - r2^2) <= 0;

flag0 = ((x-c0_x).^2 + (y-c0_y).^2 - r1^2) <= 0;
flag1 = ((x-c1_x).^2 + (y-c1_y).^2 - r2^2) <= 0;
flag2 = ((x-c2_x).^2 + (y-c2_y).^2 - r2^2) <= 0;
flag3 = ((x-c3_x).^2 + (y-c3_y).^2 - r2^2) <= 0;
flag4 = ((x-c4_x).^2 + (y-c4_y).^2 - r2^2) <= 0;

flag = flag0 + flag1 + flag2 + flag3 + flag4;

val = (flag>0) * c + (flag == 0) * c0;

end

