% This calculates the strain energy density function. 
clear
syms x y z s ro tw st s1 s2 s3 s4 b1 a1 a2 a3 a4 a5 a6 op off l L vo vo1 vo2 a11 a12 a13 a21 a22 a23 a31 a32 a33 pai real
%% primitives
%rotate & twist
frotate=[cos(ro+tw*z)*x-sin(ro+tw*z)*y;
    sin(ro+tw*z)*x+cos(ro+tw*z)*y;
    z];
%stretch
fstretch=[1/(z*(z-l)*st+1)^0.5*x;
    1/(z*(z-l)*st+1)^0.5*y;
    st/3*z^3-st*l/2*z^2+z];
%shear
fshear=[x+s1*z;y+s2*z;z];
%bend
fbend=[int(sin(a1*L/pi*(1-cos(pi*s/L))+a2*L/(2*pi)*(1-cos(2*pi*s/L))+a3*L/(3*pi)*(1-cos(3*pi*s/L))+a4*L/(4*pi)*(1-cos(4*pi*s/L))+a5*L/(5*pi)*(1-cos(5*pi*s/L))+a6*L/(6*pi)*(1-cos(6*pi*s/L))),s,0,z)+(1-(1-2*(a1*sin(pi*z/L) + a2*sin(2*pi*z/L) + a3*sin(3*pi*z/L) + a4*sin(4*pi*z/L) + a5*sin(5*pi*z/L) + a6*sin(6*pi*z/L))*x)^0.5)/(a1*sin(pi*z/L) + a2*sin(2*pi*z/L) + a3*sin(3*pi*z/L) + a4*sin(4*pi*z/L) + a5*sin(5*pi*z/L) + a6*sin(6*pi*z/L))*cos(a1*L/pi*(1-cos(pi*z/L))+a2*L/(2*pi)*(1-cos(2*pi*z/L))+a3*L/(3*pi)*(1-cos(3*pi*z/L))+a4*L/(4*pi)*(1-cos(4*pi*z/L))+a5*L/(5*pi)*(1-cos(5*pi*z/L))+a6*L/(6*pi)*(1-cos(6*pi*z/L)));
       y;
       int(cos(a1*L/pi*(1-cos(pi*s/L))+a2*L/(2*pi)*(1-cos(2*pi*s/L))+a3*L/(3*pi)*(1-cos(3*pi*s/L))+a4*L/(4*pi)*(1-cos(4*pi*s/L))+a5*L/(5*pi)*(1-cos(5*pi*s/L))+a6*L/(6*pi)*(1-cos(6*pi*s/L))),s,0,z)-(1-(1-2*(a1*sin(pi*z/L) + a2*sin(2*pi*z/L) + a3*sin(3*pi*z/L) + a4*sin(4*pi*z/L) + a5*sin(5*pi*z/L) + a6*sin(6*pi*z/L))*x)^0.5)/(a1*sin(pi*z/L) + a2*sin(2*pi*z/L) + a3*sin(3*pi*z/L) + a4*sin(4*pi*z/L) + a5*sin(5*pi*z/L) + a6*sin(6*pi*z/L))*sin(a1*L/pi*(1-cos(pi*z/L))+a2*L/(2*pi)*(1-cos(2*pi*z/L))+a3*L/(3*pi)*(1-cos(3*pi*z/L))+a4*L/(4*pi)*(1-cos(4*pi*z/L))+a5*L/(5*pi)*(1-cos(5*pi*z/L))+a6*L/(6*pi)*(1-cos(6*pi*z/L)))];
%rotate back
fbrotate=[cos(-ro)*x-sin(-ro)*y;
    sin(-ro)*x+cos(-ro)*y;
    z];

%% composition
f1 = subs(fstretch,[x,y,z],[frotate(1),frotate(2),frotate(3)]);
f2 = subs(fshear,[x,y,z],[f1(1),f1(2),f1(3)]);
f3 = subs(fbend,[x,y,z],[f2(1) f2(2) f2(3)]);

%% strain energy
J = jacobian(f3,[x y z]); %Jacobian
B = J*J.'; %deformation gradient
%principle invariants of tensor B
I1=trace(B);
I2=0.5*(trace(B)*trace(B)-trace(B*B));
I3=det(B);
%Ecoflex00-30 Mooney-Rivlin coefficients kPa
alpha=5.6;
beta=6.3;
%strain energy density function
phi=alpha*(I1-3)+beta*(I2-3);