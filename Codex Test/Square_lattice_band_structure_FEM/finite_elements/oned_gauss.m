function [r,w] = oned_gauss(rule)
%-----------------------------------------------------------------------
%  oned_gauss.m - calculate Gauss integration points on (-1,1)
%
%  Copyright (c) 2001, Jeff Borggaard, Virginia Tech
%  Version: 1.0
%
%  Usage:    [r,w] = oned_gauss(rule)
%
%  Variables:     rule  
%                        Number of Gauss points:
%                 r
%                        Gauss points located between (-1,1)      
%                 w
%                        Gauss weights corresponding to r
%-----------------------------------------------------------------------

  r = zeros(rule,1);
  w = zeros(rule,1);

  if rule == 1
    r(1) = 0;
    w(1) = 2;
  elseif rule == 2
    r(1) =-1.0d0 / sqrt(3.0d0);
    r(2) =-r(1);
    w(1) = 1.0;
    w(2) = 1.0;
  elseif rule == 3
    r(1) =-sqrt(3.0d0/5.0d0);
    r(2) = 0.0;
    r(3) =-r(1);
    w(1) = 5.0d0 / 9.0d0;
    w(2) = 8.0d0 / 9.0d0;
    w(3) = w(1);
  elseif rule == 4
    r(1) =-sqrt((3.0d0+2.0*sqrt(6.0d0/5.0d0))/7.0d0);
    r(2) =-sqrt((3.0d0-2.0*sqrt(6.0d0/5.0d0))/7.0d0);
    r(3) =-r(2);
    r(4) =-r(1);
    w(1) = 0.5d0 - 1.0d0 / ( 6.0d0 * sqrt(6.0d0/5.0d0) );
    w(2) = 0.5d0 + 1.0d0 / ( 6.0d0 * sqrt(6.0d0/5.0d0) );
    w(3) = w(2);
    w(4) = w(1);
  elseif rule == 5
    r(1) =-sqrt(5.0d0+4.0d0*sqrt(5.0d0/14.0d0)) / 3.0d0;
    r(2) =-sqrt(5.0d0-4.0d0*sqrt(5.0d0/14.0d0)) / 3.0d0;
    r(3) = 0.0d0;
    r(4) =-r(2);
    r(5) =-r(1);
    w(1) = 161.0d0/450.0d0-13.0d0/(180.d0*sqrt(5.0d0/14.0d0));
    w(2) = 161.0d0/450.0d0+13.0d0/(180.d0*sqrt(5.0d0/14.0d0));
    w(3) = 128.0d0/225.0d0;
    w(4) = w(2);
    w(5) = w(1);
  elseif rule == 6
    r(1) = -0.2386191861;
    r(2) = -0.6612093865;
    r(3) = -0.9324695142;
    r(4) = - r(1);
    r(5) = - r(2);
    r(6) = - r(3);
    w(1) = .4679139346;
    w(2) = .3607615730;
    w(3) = .1713244924;
    w(4) = w(1);
    w(5) = w(2);
    w(6) = w(3);
  else
    error('Quadrature rule not supported')
    keyboard
  end

