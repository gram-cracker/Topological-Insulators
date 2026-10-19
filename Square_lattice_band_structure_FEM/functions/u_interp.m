function [ u_xy, dux_xy, duy_xy ] = u_interp( u, p, t, x, y )
%U_INTERP returns the fem solution at (x,y)

% find the triangle that (x,y) belongs to
[ t_xy ] = fun_xy_loc( x, y, p, t );

np = length(x);
u_xy = zeros(size(x));
dux_xy = zeros(size(x));
duy_xy = zeros(size(x));

for n=1:np
    
    ne = t_xy(n);
    nodes_loc = t(ne,:);
    xi = p(nodes_loc(1),1); yi = p(nodes_loc(1),2);
    xj = p(nodes_loc(2),1); yj = p(nodes_loc(2),2);
    xk = p(nodes_loc(3),1); yk = p(nodes_loc(3),2);
 
    area_e=det([1 xi yi; 1 xj yj; 1 xk yk])/2;
    
    % u=ui*Ni(x,y)+uj*Nj(x,y)+uk*Nk(x,y)
    % where Ni(x,y)=(ai+bi*x+ci*y)/2area_e; 
    ai=xj*yk-xk*yj; bi=yj-yk; ci=xk-xj;
    aj=xk*yi-xi*yk; bj=yk-yi; cj=xi-xk;
    ak=xi*yj-xj*yi; bk=yi-yj; ck=xj-xi;
           
    Ni = (ai+bi*x(n)+ci*y(n))/(2*area_e);
    Nj = (aj+bj*x(n)+cj*y(n))/(2*area_e);
    Nk = (ak+bk*x(n)+ck*y(n))/(2*area_e);     
    
    u_xy(n) = [ Ni, Nj, Nk ] * u(nodes_loc);
    dux_xy(n) = 1/(2*area_e) * [ bi, bj, bk ] * u(nodes_loc);
    duy_xy(n) = 1/(2*area_e) * [ ci, cj, ck ] * u(nodes_loc);

end

end