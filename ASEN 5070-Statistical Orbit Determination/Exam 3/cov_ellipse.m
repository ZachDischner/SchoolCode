function cov_ellipse(P,x0,y0,C)
% Written by Zach Dischner. Wrapper for ellipse.m that parses a 2x2
% covariance matrix into proper input parameters. 


ra  =sqrt(P(1,1));
rb  = sqrt(P(2,2));
ang = P(1,2)/(ra*rb);
ellipse(ra,rb,ang,x0,y0,C);



end