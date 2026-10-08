function [ t_xy ] = fun_xy_loc( x, y, p, t )
%FUN_XY_LOC returns the triangle that xy belongs to

threshold = 0;

np = length(x);
t_xy = zeros(size(x));

% [xj, yj ]: coordinates of the three vertices for each triangle
x1 = p(t(:,1),1);  y1 = p(t(:,1),2);
x2 = p(t(:,2),1);  y2 = p(t(:,2),2);
x3 = p(t(:,3),1);  y3 = p(t(:,3),2);

xmin =  min(min(x1,x2),x3);  xmax =  max(max(x1,x2),x3);
ymin =  min(min(y1,y2),y3);  ymax =  max(max(y1,y2),y3);

for j=1:np
    
    idx_t = find(x(j)>=xmin & x(j)<=xmax & y(j)>=ymin & y(j)<=ymax);

    for k = 1:length(idx_t)
        ne = idx_t(k);
        nodes_loc = t(ne,:);
        xy_loc = p(nodes_loc,:);
        
        A1 = [ 1 x(j) y(j);   1 xy_loc(2,:);  1 xy_loc(3,:) ];
        A2 = [ 1 xy_loc(1,:); 1 x(j) y(j);    1 xy_loc(3,:) ];
        A3 = [ 1 xy_loc(1,:); 1 xy_loc(2,:);   1 x(j) y(j)  ];
        
        % if ( det(A1)>=0 & (det(A2)>=0) & (det(A3)>=0) )
        if ( det(A1)>=threshold & (det(A2)>=threshold) & (det(A3)>=threshold) )            
            t_xy(j) = ne;            
            % plot( [xy_loc(:,1); xy_loc(1,1)], [xy_loc(:,2);xy_loc(1,2)] );
            % hold on; plot(x(j),y(j),'-o');
        end                
        
    end

end

end

