% This computes the mode parameters p that satisfies the desired end position and orientation using the Lie group method (Kim J S, Chirikjian G S. Conformational analysis of stiff chiral polymers with end-constraints[J]. Molecular Simulation, 2006, 32(14): 1139-1154.). This is for the block and film case.
clear
tic
w=9;%number of modes
p=zeros(w,1)+[0;0;0;0;0;0.001;0;0;0];%initial guess of mode parameters
% desired orientation
alpha_fdesired=0.5;
beta_fdesired =0.5;
gamma_fdesired=0.5;
% rotation matrix
A_fdesired=[cos(gamma_fdesired),-sin(gamma_fdesired),0;sin(gamma_fdesired),cos(gamma_fdesired),0;0,0,1]*...
    [cos(beta_fdesired),0,sin(beta_fdesired);0,1,0;-sin(beta_fdesired),0,cos(beta_fdesired)]*...
    [1,0,0;0,cos(alpha_fdesired),-sin(alpha_fdesired);0,sin(alpha_fdesired),cos(alpha_fdesired)];
%desired position
a_fdesired=[0.2;0.5;1.1];
% homogeneous transformation matrix
gd=[A_fdesired,a_fdesired;0,0,0,1];
W=eye(w);%weight matrix in jacobian pseudoinverse

max=10;% number of steps
deltatk=0.1; % step size
epsilon=10^-3;
for step=1:max
    tk=(step-1)*deltatk;

g0=get_g(p);
A=g0(1:3,1:3);
aa=g0(1:3,4);
E=epsilon*eye(w);
ginv=[A',-A'*aa;0,0,0,1];
for i=1:w
    gds=(get_g(p+E(:,i))-g0)/epsilon;
    ggd=ginv*gds;
    Jacobian(:,i)=[ggd(3,2);ggd(1,3);ggd(2,1);ggd(1,4);ggd(2,4);ggd(3,4)];
end

gpt=[A*expm(tk*logm(A'*A_fdesired)),aa+tk*(a_fdesired-aa);0 0 0 1];
dgpt=[A*logm(A'*A_fdesired)*expm(tk*logm(A'*A_fdesired)),a_fdesired-aa;0 0 0 0];
gptinv=[(A*expm(tk*logm(A'*A_fdesired)))',-(A*expm(tk*logm(A'*A_fdesired)))*(aa+tk*(a_fdesired-aa));0 0 0 1];
Adfirst=ginv*gpt;
Ad=[Adfirst(1:3,1:3),zeros(3);[0,-Adfirst(3,4),Adfirst(2,4);Adfirst(3,4),0,-Adfirst(1,4);-Adfirst(2,4),Adfirst(1,4),0],Adfirst(1:3,1:3)];
vfirst=gptinv*dgpt;
vsecond=logm(ginv*gpt);
detadt=W\Jacobian'/(Jacobian*inv(W)*Jacobian'+epsilon*eye(6))*Ad*[-vfirst(2,3);vfirst(1,3);-vfirst(1,2);vfirst(1,4);vfirst(2,4);vfirst(3,4)];
pc=W\Jacobian'/(Jacobian*inv(W)*Jacobian'+epsilon*eye(6))*Ad*[-vsecond(2,3);vsecond(1,3);-vsecond(1,2);vsecond(1,4);vsecond(2,4);vsecond(3,4)];
p=p+deltatk*detadt+pc;
end
% g0=get_g(p);
% A=g0(1:3,1:3);
% aa=g0(1:3,4);
% A_fdesired-A
% a_fdesired-aa%check error
toc

function [g]=get_g(p)
%This function calculates the homogeneous transformation matrix
ro=p(1);
tw=p(2);
st=p(3);
s1=p(4);
s2=p(5);
a1=p(6);
a2=p(7);
a3=p(8);
a4=p(9);

l=1;%original length of the block
L=st/3*l^3-st*l/2*l^2+l;% length of the block after stretch
func = @(s) cos(a1*L/pi*(1-cos(pi*s/L))+a2*L/(2*pi)*(1-cos(2*pi*s/L))+a3*L/(3*pi)*(1-cos(3*pi*s/L))+a4*L/(4*pi)*(1-cos(4*pi*s/L)));
funs = @(s) sin(a1*L/pi*(1-cos(pi*s/L))+a2*L/(2*pi)*(1-cos(2*pi*s/L))+a3*L/(3*pi)*(1-cos(3*pi*s/L))+a4*L/(4*pi)*(1-cos(4*pi*s/L)));

% compute the orientation of the top plane and the position of the center point of the top plane
z=1-10^-10;
ifuns = integral(funs, 0, (st*z^3)/3 - (l*st*z^2)/2 + z);
ifunc = integral(func, 0, (st*z^3)/3 - (l*st*z^2)/2 + z);
x=0;
y=0;
X=cos(ro)*(ifuns - (cos((L*a1*(cos((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/pi + (L*a2*(cos((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(2*pi) + (L*a3*(cos((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(3*pi) + (L*a4*(cos((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(4*pi))*((1 - ((x*cos(ro + tw*z) - y*sin(ro + tw*z))/(1 - st*z*(l - z))^(1/2) + s1*((st*z^3)/3 - (l*st*z^2)/2 + z))*(2*a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L)))^(1/2) - 1))/(a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L))) + sin(ro)*((y*cos(ro + tw*z) + x*sin(ro + tw*z))/(1 - st*z*(l - z))^(1/2) + s2*((st*z^3)/3 - (l*st*z^2)/2 + z));
Y=cos(ro)*((y*cos(ro + tw*z) + x*sin(ro + tw*z))/(1 - st*z*(l - z))^(1/2) + s2*((st*z^3)/3 - (l*st*z^2)/2 + z)) - sin(ro)*(ifuns - (cos((L*a1*(cos((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/pi + (L*a2*(cos((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(2*pi) + (L*a3*(cos((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(3*pi) + (L*a4*(cos((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(4*pi))*((1 - ((x*cos(ro + tw*z) - y*sin(ro + tw*z))/(1 - st*z*(l - z))^(1/2) + s1*((st*z^3)/3 - (l*st*z^2)/2 + z))*(2*a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L)))^(1/2) - 1))/(a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L)));
Z=ifunc - (sin((L*a1*(cos((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/pi + (L*a2*(cos((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(2*pi) + (L*a3*(cos((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(3*pi) + (L*a4*(cos((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(4*pi))*((1 - ((x*cos(ro + tw*z) - y*sin(ro + tw*z))/(1 - st*z*(l - z))^(1/2) + s1*((st*z^3)/3 - (l*st*z^2)/2 + z))*(2*a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L)))^(1/2) - 1))/(a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L));

x=0.05;
y=0;
X1=cos(ro)*(ifuns - (cos((L*a1*(cos((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/pi + (L*a2*(cos((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(2*pi) + (L*a3*(cos((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(3*pi) + (L*a4*(cos((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(4*pi))*((1 - ((x*cos(ro + tw*z) - y*sin(ro + tw*z))/(1 - st*z*(l - z))^(1/2) + s1*((st*z^3)/3 - (l*st*z^2)/2 + z))*(2*a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L)))^(1/2) - 1))/(a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L))) + sin(ro)*((y*cos(ro + tw*z) + x*sin(ro + tw*z))/(1 - st*z*(l - z))^(1/2) + s2*((st*z^3)/3 - (l*st*z^2)/2 + z));
Y1=cos(ro)*((y*cos(ro + tw*z) + x*sin(ro + tw*z))/(1 - st*z*(l - z))^(1/2) + s2*((st*z^3)/3 - (l*st*z^2)/2 + z)) - sin(ro)*(ifuns - (cos((L*a1*(cos((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/pi + (L*a2*(cos((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(2*pi) + (L*a3*(cos((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(3*pi) + (L*a4*(cos((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(4*pi))*((1 - ((x*cos(ro + tw*z) - y*sin(ro + tw*z))/(1 - st*z*(l - z))^(1/2) + s1*((st*z^3)/3 - (l*st*z^2)/2 + z))*(2*a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L)))^(1/2) - 1))/(a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L)));
Z1=ifunc - (sin((L*a1*(cos((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/pi + (L*a2*(cos((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(2*pi) + (L*a3*(cos((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(3*pi) + (L*a4*(cos((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(4*pi))*((1 - ((x*cos(ro + tw*z) - y*sin(ro + tw*z))/(1 - st*z*(l - z))^(1/2) + s1*((st*z^3)/3 - (l*st*z^2)/2 + z))*(2*a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L)))^(1/2) - 1))/(a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L));

x=0;
y=0.05;
X2=cos(ro)*(ifuns - (cos((L*a1*(cos((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/pi + (L*a2*(cos((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(2*pi) + (L*a3*(cos((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(3*pi) + (L*a4*(cos((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(4*pi))*((1 - ((x*cos(ro + tw*z) - y*sin(ro + tw*z))/(1 - st*z*(l - z))^(1/2) + s1*((st*z^3)/3 - (l*st*z^2)/2 + z))*(2*a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L)))^(1/2) - 1))/(a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L))) + sin(ro)*((y*cos(ro + tw*z) + x*sin(ro + tw*z))/(1 - st*z*(l - z))^(1/2) + s2*((st*z^3)/3 - (l*st*z^2)/2 + z));
Y2=cos(ro)*((y*cos(ro + tw*z) + x*sin(ro + tw*z))/(1 - st*z*(l - z))^(1/2) + s2*((st*z^3)/3 - (l*st*z^2)/2 + z)) - sin(ro)*(ifuns - (cos((L*a1*(cos((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/pi + (L*a2*(cos((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(2*pi) + (L*a3*(cos((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(3*pi) + (L*a4*(cos((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(4*pi))*((1 - ((x*cos(ro + tw*z) - y*sin(ro + tw*z))/(1 - st*z*(l - z))^(1/2) + s1*((st*z^3)/3 - (l*st*z^2)/2 + z))*(2*a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L)))^(1/2) - 1))/(a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L)));
Z2=ifunc - (sin((L*a1*(cos((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/pi + (L*a2*(cos((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(2*pi) + (L*a3*(cos((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(3*pi) + (L*a4*(cos((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) - 1))/(4*pi))*((1 - ((x*cos(ro + tw*z) - y*sin(ro + tw*z))/(1 - st*z*(l - z))^(1/2) + s1*((st*z^3)/3 - (l*st*z^2)/2 + z))*(2*a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + 2*a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L)))^(1/2) - 1))/(a1*sin((pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a2*sin((2*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a3*sin((3*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L) + a4*sin((4*pi*((st*z^3)/3 - (l*st*z^2)/2 + z))/L));

E1=[X1-X;Y1-Y;Z1-Z]/norm([X1-X;Y1-Y;Z1-Z]);
E2=[X2-X;Y2-Y;Z2-Z]/norm([X2-X;Y2-Y;Z2-Z]);
E3=cross(E1,E2);
R=[E1,E2,E3];

length=[X;Y;Z];
g=[R,length;0 0 0 1];% homogeneous transformation matrix
end