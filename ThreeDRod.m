% 3D rod, Lie group approach in "Kim J S, Chirikjian G S. Conformational analysis of stiff chiral polymers with end-constraints[J]. Molecular Simulation, 2006, 32(14): 1139-1154.".
clear
tic % start timing
W=eye(6); % weight matrix
a_fdesired=[-0.3;-0.1;0.6]; % desired position
alpha_fdesired=0.5; % desired orientation
beta_fdesired=0.5;
gamma_fdesired=0.5;
A_fdesired=[cos(gamma_fdesired),-sin(gamma_fdesired),0;sin(gamma_fdesired),cos(gamma_fdesired),0;0,0,1]*...
    [cos(beta_fdesired),0,sin(beta_fdesired);0,1,0;-sin(beta_fdesired),0,cos(beta_fdesired)]*...
    [1,0,0;0,cos(alpha_fdesired),-sin(alpha_fdesired);0,sin(alpha_fdesired),cos(alpha_fdesired)]; % rotation matrix
gd=[A_fdesired,a_fdesired;0,0,0,1]; 
eta=[0;0;0;0;0;0]; % original angular velocity and lagrange multipliers
A0=eye(3); % original orientation
max=20; %number of steps
deltatk=0.05; % 1/max
epsilon=10^-6;

for step=1:max
    tk=(step-1)*deltatk;
    
tspan=[0:0.005:1];%range of s
y0=[eta(1);eta(2);eta(3);A0(1,1);A0(1,2);A0(1,3);A0(2,1);A0(2,2);A0(2,3);A0(3,1);A0(3,2);A0(3,3)];%initial values of variables
[~,y]=ode45(@(t,y) odefuni(t,y,eta(4:6)),tspan,y0); %solve ode
aa(:,1)=[trapz(0:0.005:1,y(1:201,6));trapz(0:0.005:1,y(1:201,9));trapz(0:0.005:1,y(1:201,12))];


y0=[eta(1)+epsilon;eta(2);eta(3);A0(1,1);A0(1,2);A0(1,3);A0(2,1);A0(2,2);A0(2,3);A0(3,1);A0(3,2);A0(3,3)];%initial values of variables
[~,y1]=ode45(@(t,y) odefuni(t,y,eta(4:6)),tspan,y0); %solve ode
aa1(:,1)=[trapz(0:0.005:1,y1(1:201,6));trapz(0:0.005:1,y1(1:201,9));trapz(0:0.005:1,y1(1:201,12))];
y0=[eta(1);eta(2)+epsilon;eta(3);A0(1,1);A0(1,2);A0(1,3);A0(2,1);A0(2,2);A0(2,3);A0(3,1);A0(3,2);A0(3,3)];%initial values of variables
[~,y2]=ode45(@(t,y) odefuni(t,y,eta(4:6)),tspan,y0); %solve ode
aa2(:,1)=[trapz(0:0.005:1,y2(1:201,6));trapz(0:0.005:1,y2(1:201,9));trapz(0:0.005:1,y2(1:201,12))];
y0=[eta(1);eta(2);eta(3)+epsilon;A0(1,1);A0(1,2);A0(1,3);A0(2,1);A0(2,2);A0(2,3);A0(3,1);A0(3,2);A0(3,3)];%initial values of variables
[~,y3]=ode45(@(t,y) odefuni(t,y,eta(4:6)),tspan,y0); %solve ode
aa3(:,1)=[trapz(0:0.005:1,y3(1:201,6));trapz(0:0.005:1,y3(1:201,9));trapz(0:0.005:1,y3(1:201,12))];
y0=[eta(1);eta(2);eta(3);A0(1,1);A0(1,2);A0(1,3);A0(2,1);A0(2,2);A0(2,3);A0(3,1);A0(3,2);A0(3,3)];%initial values of variables
[~,y4]=ode45(@(t,y) odefuni(t,y,[eta(4)+epsilon;eta(5);eta(6)]),tspan,y0); %solve ode
aa4(:,1)=[trapz(0:0.005:1,y4(1:201,6));trapz(0:0.005:1,y4(1:201,9));trapz(0:0.005:1,y4(1:201,12))];
y0=[eta(1);eta(2);eta(3);A0(1,1);A0(1,2);A0(1,3);A0(2,1);A0(2,2);A0(2,3);A0(3,1);A0(3,2);A0(3,3)];%initial values of variables
[~,y5]=ode45(@(t,y) odefuni(t,y,[eta(4);eta(5)+epsilon;eta(6)]),tspan,y0); %solve ode
aa5(:,1)=[trapz(0:0.005:1,y5(1:201,6));trapz(0:0.005:1,y5(1:201,9));trapz(0:0.005:1,y5(1:201,12))];
y0=[eta(1);eta(2);eta(3);A0(1,1);A0(1,2);A0(1,3);A0(2,1);A0(2,2);A0(2,3);A0(3,1);A0(3,2);A0(3,3)];%initial values of variables
[~,y6]=ode45(@(t,y) odefuni(t,y,[eta(4);eta(5);eta(6)+epsilon]),tspan,y0); %solve ode
aa6(:,1)=[trapz(0:0.005:1,y6(1:201,6));trapz(0:0.005:1,y6(1:201,9));trapz(0:0.005:1,y6(1:201,12))];

A=[y(end,4),y(end,5),y(end,6);y(end,7),y(end,8),y(end,9);y(end,10),y(end,11),y(end,12)];%Rotation matrix
g=[A,aa;[0 0 0],1];%aa: translation vector
ginv=[[y(end,4),y(end,5),y(end,6);y(end,7),y(end,8),y(end,9);y(end,10),y(end,11),y(end,12)]',-[y(end,4),y(end,5),y(end,6);y(end,7),y(end,8),y(end,9);y(end,10),y(end,11),y(end,12)]'*aa;[0 0 0],1];% inverse of g
dg1=1/epsilon*[[y1(end,4),y1(end,5),y1(end,6);y1(end,7),y1(end,8),y1(end,9);y1(end,10),y1(end,11),y1(end,12)]-[y(end,4),y(end,5),y(end,6);y(end,7),y(end,8),y(end,9);y(end,10),y(end,11),y(end,12)],aa1-aa;0 0 0 0];
dg2=1/epsilon*[[y2(end,4),y2(end,5),y2(end,6);y2(end,7),y2(end,8),y2(end,9);y2(end,10),y2(end,11),y2(end,12)]-[y(end,4),y(end,5),y(end,6);y(end,7),y(end,8),y(end,9);y(end,10),y(end,11),y(end,12)],aa2-aa;0 0 0 0];
dg3=1/epsilon*[[y3(end,4),y3(end,5),y3(end,6);y3(end,7),y3(end,8),y3(end,9);y3(end,10),y3(end,11),y3(end,12)]-[y(end,4),y(end,5),y(end,6);y(end,7),y(end,8),y(end,9);y(end,10),y(end,11),y(end,12)],aa3-aa;0 0 0 0];
dg4=1/epsilon*[[y4(end,4),y4(end,5),y4(end,6);y4(end,7),y4(end,8),y4(end,9);y4(end,10),y4(end,11),y4(end,12)]-[y(end,4),y(end,5),y(end,6);y(end,7),y(end,8),y(end,9);y(end,10),y(end,11),y(end,12)],aa4-aa;0 0 0 0];
dg5=1/epsilon*[[y5(end,4),y5(end,5),y5(end,6);y5(end,7),y5(end,8),y5(end,9);y5(end,10),y5(end,11),y5(end,12)]-[y(end,4),y(end,5),y(end,6);y(end,7),y(end,8),y(end,9);y(end,10),y(end,11),y(end,12)],aa5-aa;0 0 0 0];
dg6=1/epsilon*[[y6(end,4),y6(end,5),y6(end,6);y6(end,7),y6(end,8),y6(end,9);y6(end,10),y6(end,11),y6(end,12)]-[y(end,4),y(end,5),y(end,6);y(end,7),y(end,8),y(end,9);y(end,10),y(end,11),y(end,12)],aa6-aa;0 0 0 0];

j1=ginv*dg1;
j2=ginv*dg2;
j3=ginv*dg3;
j4=ginv*dg4;
j5=ginv*dg5;
j6=ginv*dg6;

Jacobian=[-j1(2,3),-j2(2,3),-j3(2,3),-j4(2,3),-j5(2,3),-j6(2,3);...
    j1(1,3),j2(1,3),j3(1,3),j4(1,3),j5(1,3),j6(1,3);...
    -j1(1,2),-j2(1,2),-j3(1,2),-j4(1,2),-j5(1,2),-j6(1,2);...
    j1(1,4),j2(1,4),j3(1,4),j4(1,4),j5(1,4),j6(1,4);...
    j1(2,4),j2(2,4),j3(2,4),j4(2,4),j5(2,4),j6(2,4);...
    j1(3,4),j2(3,4),j3(3,4),j4(3,4),j5(3,4),j6(3,4)];

gpt=[A*expm(tk*logm(A'*A_fdesired)),aa+tk*(a_fdesired-aa);0 0 0 1];
dgpt=[A*logm(A'*A_fdesired)*expm(tk*logm(A'*A_fdesired)),a_fdesired-aa;0 0 0 0];
gptinv=[(A*expm(tk*logm(A'*A_fdesired)))',-(A*expm(tk*logm(A'*A_fdesired)))*(aa+tk*(a_fdesired-aa));0 0 0 1];
Adfirst=ginv*gpt;
Ad=[Adfirst(1:3,1:3),zeros(3);[0,-Adfirst(3,4),Adfirst(2,4);Adfirst(3,4),0,-Adfirst(1,4);-Adfirst(2,4),Adfirst(1,4),0]*Adfirst(1:3,1:3),Adfirst(1:3,1:3)];
vfirst=gptinv*dgpt;
vsecond=logm(ginv*gpt);
detadt=inv(W)*Jacobian'*inv(Jacobian*inv(W)*Jacobian'+epsilon*eye(6))*Ad*[-vfirst(2,3);vfirst(1,3);-vfirst(1,2);vfirst(1,4);vfirst(2,4);vfirst(3,4)];
etac=inv(W)*Jacobian'*inv(Jacobian*inv(W)*Jacobian'+epsilon*eye(6))*Ad*[-vsecond(2,3);vsecond(1,3);-vsecond(1,2);vsecond(1,4);vsecond(2,4);vsecond(3,4)];
eta=eta+deltatk*detadt+etac;
end
toc

% plot result
row=201;
for i=2:row
range=linspace(0,i/row,i);
a(:,i)=[trapz(range,y(1:i,6));trapz(range,y(1:i,9));trapz(range,y(1:i,12))];
end
plot3(a(1,:),a(2,:),a(3,:),'r')
axis([-1 1 -1 1 -1 1])
hold on


function dydt=odefuni(t,y,lambda)
B=eye(3);
b=[0;0;0];
dw=inv(B)*(-cross([y(1);y(2);y(3)],B*[y(1);y(2);y(3)]-b)+[-lambda(1)*y(5)-lambda(2)*y(8)-lambda(3)*y(11);lambda(1)*y(4)+lambda(2)*y(7)+lambda(3)*y(10);0]);
dw1=dw(1);
dw2=dw(2);
dw3=dw(3);
da11=y(3)*y(5)-y(2)*y(6);
da12=-y(3)*y(4)+y(1)*y(6);
da13=y(2)*y(4)-y(1)*y(5);
da21=y(3)*y(8)-y(2)*y(9);
da22=-y(3)*y(7)+y(1)*y(9);
da23=y(2)*y(7)-y(1)*y(8);
da31=y(3)*y(11)-y(2)*y(12);
da32=-y(3)*y(10)+y(1)*y(12);
da33=y(2)*y(10)-y(1)*y(11);
dydt=[dw1;dw2;dw3;da11;da12;da13;da21;da22;da23;da31;da32;da33];
end
%}