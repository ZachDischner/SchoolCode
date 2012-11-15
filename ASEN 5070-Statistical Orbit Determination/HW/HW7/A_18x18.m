function A = A_18x18(Area,Cd,H,J2,R_E,m,r,r0,rho0,theta_dot,uE,x,xdot,y,ydot,z,zdot)
%A_18X18
%    A = A_18X18(AREA,CD,H,J2,R_E,M,R,R0,RHO0,THETA_DOT,UE,X,XDOT,Y,YDOT,Z,ZDOT)

%    This function was generated by the Symbolic Math Toolbox version 5.9.
%    25-Oct-2012 13:25:54

t2 = 1.0./r.^2;
t3 = r.^2;
t4 = 1.0./t3.^(3.0./2.0);
t5 = R_E.^2;
t6 = z.^2;
t7 = t2.*t6.*5.0;
t8 = t7-1.0;
t9 = J2.*t2.*t5.*t8.*(3.0./2.0);
t10 = t9-1.0;
t11 = theta_dot.*y;
t12 = t11+xdot;
t20 = theta_dot.*x;
t13 = -t20+ydot;
t14 = 1.0./H;
t15 = 1.0./m;
t16 = r-r0;
t17 = t14.*t16;
t18 = exp(t17);
t19 = t12.^2;
t21 = t13.^2;
t22 = zdot.^2;
t23 = t19+t21+t22;
t24 = 1.0./r.^4;
t25 = 1.0./r.^6;
t26 = 1.0./t3.^(5.0./2.0);
t27 = sqrt(t23);
t28 = 1.0./sqrt(t23);
t29 = 1.0./sqrt(t3);
t30 = t3.^2;
t31 = t30.^2;
t32 = 1.0./t3.^(9.0./2.0);
t33 = x.^2;
t34 = theta_dot.^2;
t35 = 1.0./t3.^(7.0./2.0);
t36 = y.^2;
t37 = J2.*t5.*t8.*t24.*x.*3.0;
t38 = J2.*t5.*t6.*t25.*x.*1.5e1;
t39 = t37+t38;
t40 = t4.*t10.*uE;
t41 = J2.*t5.*t8.*t24.*y.*3.0;
t42 = J2.*t5.*t6.*t25.*y.*1.5e1;
t43 = t41+t42;
t44 = Area.*Cd.*rho0.*t12.*t13.*t15.*t18.*t28.*theta_dot.*(1.0./2.0);
t45 = xdot.^2;
t46 = ydot.^2;
t47 = t30.*2.0;
t48 = J2.*t3.*t5.*3.0;
t57 = J2.*t5.*t6.*1.5e1;
t49 = t47+t48-t57;
t50 = t6.*-4.0+t33+t36;
t51 = t7-3.0;
t52 = J2.*t2.*t5.*t51.*(3.0./2.0);
t53 = t52-1.0;
t54 = t33.*t34;
t55 = t34.*t36;
t56 = theta_dot.*xdot.*y.*2.0;
A = reshape([0.0,0.0,0.0,t40+t44-t10.*t26.*t33.*uE.*3.0-t4.*t39.*uE.*x-Area.*Cd.*rho0.*t12.*t14.*t15.*t18.*t27.*t29.*x.*(1.0./2.0),-t4.*t39.*uE.*y-t10.*t26.*uE.*x.*y.*3.0+Area.*Cd.*rho0.*t15.*t18.*t27.*theta_dot.*(1.0./2.0)+Area.*Cd.*rho0.*t15.*t18.*t21.*t28.*theta_dot.*(1.0./2.0)-Area.*Cd.*rho0.*t13.*t14.*t15.*t18.*t27.*t29.*x.*(1.0./2.0),-t4.*uE.*z.*(t38+J2.*t5.*t24.*t51.*x.*3.0)-t26.*t53.*uE.*x.*z.*3.0+Area.*Cd.*rho0.*t13.*t15.*t18.*t28.*theta_dot.*zdot.*(1.0./2.0)-Area.*Cd.*rho0.*t14.*t15.*t18.*t27.*t29.*x.*zdot.*(1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t4.*t43.*uE.*x-t10.*t26.*uE.*x.*y.*3.0-Area.*Cd.*rho0.*t15.*t18.*t27.*theta_dot.*(1.0./2.0)-Area.*Cd.*rho0.*t15.*t18.*t19.*t28.*theta_dot.*(1.0./2.0)-Area.*Cd.*rho0.*t12.*t14.*t15.*t18.*t27.*t29.*y.*(1.0./2.0),t40-t44-t10.*t26.*t36.*uE.*3.0-t4.*t43.*uE.*y-Area.*Cd.*rho0.*t13.*t14.*t15.*t18.*t27.*t29.*y.*(1.0./2.0),-t4.*uE.*z.*(t42+J2.*t5.*t24.*t51.*y.*3.0)-t26.*t53.*uE.*y.*z.*3.0-Area.*Cd.*rho0.*t12.*t15.*t18.*t28.*theta_dot.*zdot.*(1.0./2.0)-Area.*Cd.*rho0.*t14.*t15.*t18.*t27.*t29.*y.*zdot.*(1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t15.*t32.*z.*(m.*t30.*uE.*x.*6.0+J2.*m.*t3.*t5.*uE.*x.*4.5e1-J2.*m.*t5.*t6.*uE.*x.*1.05e2).*(1.0./2.0)-t14.*t15.*t32.*z.*(Area.*Cd.*rho0.*t18.*t27.*t31.*xdot+Area.*Cd.*rho0.*t18.*t27.*t31.*theta_dot.*y).*(1.0./2.0),t15.*t32.*z.*(m.*t30.*uE.*y.*6.0+J2.*m.*t3.*t5.*uE.*y.*4.5e1-J2.*m.*t5.*t6.*uE.*y.*1.05e2).*(1.0./2.0)-t14.*t15.*t32.*z.*(Area.*Cd.*rho0.*t18.*t27.*t31.*ydot-Area.*Cd.*rho0.*t18.*t27.*t31.*theta_dot.*x).*(1.0./2.0),t15.*t32.*(m.*t3.*t30.*uE.*2.0-m.*t6.*t30.*uE.*6.0+J2.*m.*t5.*t30.*uE.*9.0+J2.*m.*t5.*t6.^2.*uE.*1.05e2-J2.*m.*t3.*t5.*t6.*uE.*9.0e1).*(-1.0./2.0)-Area.*Cd.*rho0.*t14.*t15.*t18.*t27.*t31.*t32.*z.*zdot.*(1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,Area.*Cd.*rho0.*t15.*t18.*t28.*(t22+t45.*2.0+t46+t54+t34.*t36.*2.0-theta_dot.*x.*ydot.*2.0+theta_dot.*xdot.*y.*4.0).*(-1.0./2.0),Area.*Cd.*rho0.*t12.*t13.*t15.*t18.*t28.*(-1.0./2.0),Area.*Cd.*rho0.*t12.*t15.*t18.*t28.*zdot.*(-1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,Area.*Cd.*rho0.*t12.*t13.*t15.*t18.*t28.*(-1.0./2.0),Area.*Cd.*rho0.*t15.*t18.*t28.*(t22+t45+t46.*2.0+t55+t56+t33.*t34.*2.0-theta_dot.*x.*ydot.*4.0).*(-1.0./2.0),Area.*Cd.*rho0.*t13.*t15.*t18.*t28.*zdot.*(-1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,Area.*Cd.*rho0.*t12.*t15.*t18.*t28.*zdot.*(-1.0./2.0),Area.*Cd.*rho0.*t13.*t15.*t18.*t28.*zdot.*(-1.0./2.0),Area.*Cd.*rho0.*t15.*t18.*t28.*(t22.*2.0+t45+t46+t54+t55+t56-theta_dot.*x.*ydot.*2.0).*(-1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t35.*t49.*x.*(-1.0./2.0),t35.*t49.*y.*(-1.0./2.0),t35.*z.*(t47-t57+J2.*t3.*t5.*9.0).*(-1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t35.*t50.*uE.*x.*(-3.0./2.0),t5.*t35.*t50.*uE.*y.*(-3.0./2.0),t5.*t35.*uE.*z.*(t6.*-2.0+t33.*3.0+t36.*3.0).*(-3.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,Area.*rho0.*t12.*t15.*t18.*t27.*(-1.0./2.0),Area.*rho0.*t13.*t15.*t18.*t27.*(-1.0./2.0),Area.*rho0.*t15.*t18.*t27.*zdot.*(-1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[18, 18]);