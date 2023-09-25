%Inflate a chamber. 
clear
tic
Vd=80; % desired input volume
w=5; % number of modes
max=20; % number of steps
p=zeros(w,1); % original mode parameters
epsilon=10^-5; 
E=epsilon*eye(w);
%weight matrix in jacobian pseudoinverse. It is computed using code GetWeightMatrix.m
W=[0.132729571256575,0.00321559961454290,0.00779796061466052,0.0106157358317054,0.0110695448411330;0.00321559961454290,0.216377548448813,0.0216686218949095,0.0308517375523768,0.0452669760350615;0.00779796061466052,0.0216686218949095,0.376464687350446,0.0557203179586730,0.0747717294861613;0.0106157358317054,0.0308517375523768,0.0557203179586730,0.624584764745024,0.106193166394618;0.0110695448411330,0.0452669760350615,0.0747717294861613,0.106193166394618,0.952060687096547];

for step=1:max
phi0=get_phi(p,Vd);
delta_phi=-phi0/(max-step+1)^0.1;
for i=1:w
    Jacobian(1,i)=(get_phi(p+E(:,i),Vd)-phi0)/epsilon;
end
ji=inv(W)*Jacobian'*inv(Jacobian*inv(W)*Jacobian');
p=p+ji*delta_phi;
end
toc

function [phi]=get_phi(p,Vd)
% This function computes the error between current volume and desired volume
a1=p(1);
a2=p(2);
a3=p(3);
a4=p(4);
a5=p(5);
L=10; % height of chamber
r0=1.6; % radius of chamber
Vo=r0^2*pi*L; % original volume of chamber
Vt=Vo+Vd; % Desired volume of chamber after deformation
fun = @(s) pi*(r0^2+a1*sin(pi*s/L)+a2*sin(3*pi*s/L)+a3*sin(5*pi*s/L)+a4*sin(7*pi*s/L)+a5*sin(9*pi*s/L));
volume=integral(fun,0,L);% Volume of chamber
phi=(Vt-volume)^2; % error between current volume and desired volume
end