%Rod with gravity. Projected gradient method
clear
tic

epsilon=10^-5;
max=40;% number of iterations
yf_desired=0.45;
xf_desired=0.70;
orif_desired=0;
n=6;%number of modes
length=1;% arclength of curve
a=0.0001*[1*(rand-0.5);1*(rand-0.5);1*(rand-0.5);1*(rand-0.5);1*(rand-0.5);1*(rand-0.5)];%initial guess of mode parameters

forward_x = @(a1,a2,a3,a4,a5,a6,x) cos(a1*length/pi*(1-cos(pi*x/length))+a2*length/(2*pi)*(1-cos(2*pi*x/length))+a3*length/(3*pi)*(1-cos(3*pi*x/length))+a4*length/(4*pi)*(1-cos(4*pi*x/length))+a5*length/(5*pi)*(1-cos(5*pi*x/length))+a6*length/(6*pi)*(1-cos(6*pi*x/length)));
forward_y = @(a1,a2,a3,a4,a5,a6,x) sin(a1*length/pi*(1-cos(pi*x/length))+a2*length/(2*pi)*(1-cos(2*pi*x/length))+a3*length/(3*pi)*(1-cos(3*pi*x/length))+a4*length/(4*pi)*(1-cos(4*pi*x/length))+a5*length/(5*pi)*(1-cos(5*pi*x/length))+a6*length/(6*pi)*(1-cos(6*pi*x/length)));
for iteration=1:max
Jacobian=jacob(n,a);
[forwardkin] = forward(a);%this computes the current end point position
x_current=forwardkin(1);
y_current=forwardkin(2);
ori_current=a(1)*length/pi*2+a(3)*length/(3*pi)*2+a(5)*length/(5*pi)*2;

energy0=get_energy(a);
E=eye(n);
for i=1:n
    g(i)=(get_energy(a+epsilon*E(:,i))-energy0)/epsilon;%gradient of energy
end

%get desired end point position for current iteration
x_desired(iteration)=x_current+(xf_desired-x_current)/(max-iteration+1);
y_desired(iteration)=y_current+(yf_desired-y_current)/(max-iteration+1);
ori_desired(iteration)=ori_current+(orif_desired-ori_current)/(max-iteration+1);

lambdas=0.00001;
invjacob=Jacobian'*inv(Jacobian*Jacobian'+lambdas*eye(3));
Z=eye(n)-invjacob*Jacobian;
    
c=[x_desired(iteration)-x_current;y_desired(iteration)-y_current;ori_desired(iteration)-ori_current];
middle=invjacob*c;
alpha=1;
last=alpha*Z*g';
a(:)=a(:)+middle-last;

end
toc

%plot result
i=0;
for s=0:0.001:length
    i=i+1;
    y1(i)=integral(@(x) forward_x(a(1),a(2),a(3),a(4),a(5),a(6),x),0,s);
    y2(i)=integral(@(x) forward_y(a(1),a(2),a(3),a(4),a(5),a(6),x),0,s);            
end
plot(1*y2,1*y1,'b','Linewidth',1.0,'color',[0 0.6 0])
axis([-1 1 -1 1])
hold on



function [forwardkin] = forward(a)
length=1;
forward_x = @(a,x) cos(a(1)*length/pi*(1-cos(pi*x/length))+a(2)*length/(2*pi)*(1-cos(2*pi*x/length))+a(3)*length/(3*pi)*(1-cos(3*pi*x/length))+a(4)*length/(4*pi)*(1-cos(4*pi*x/length))+a(5)*length/(5*pi)*(1-cos(5*pi*x/length))+a(6)*length/(6*pi)*(1-cos(6*pi*x/length)));
forward_y = @(a,x) sin(a(1)*length/pi*(1-cos(pi*x/length))+a(2)*length/(2*pi)*(1-cos(2*pi*x/length))+a(3)*length/(3*pi)*(1-cos(3*pi*x/length))+a(4)*length/(4*pi)*(1-cos(4*pi*x/length))+a(5)*length/(5*pi)*(1-cos(5*pi*x/length))+a(6)*length/(6*pi)*(1-cos(6*pi*x/length)));
forwardkin=[integral(@(x) forward_x(a,x),0,length);integral(@(x) forward_y(a,x),0,length)];
end

function [Jacobian] = jacob(n,a)
length=1;
epsilon=10^-5;
forwardkin0 = forward(a);
for count=1:n
    position=zeros(n,1);
    position(count)=1;
    forwardkin1 = forward(a+epsilon*position);
    Jacobian1(:,count)=1/epsilon*(forwardkin1-forwardkin0);
end
    Jacobian2=[1*length/pi*2,0,1*length/(3*pi)*2,0,1*length/(5*pi)*2,0];
    Jacobian=[Jacobian1;Jacobian2];
end


function [energy]=get_energy(a)
a1=a(1);
a2=a(2);
a3=a(3);
a4=a(4);
a5=a(5);
a6=a(6);
%
r0=[0.00 0.01 0.02];
thetao=0:pi:3/2*pi;
zo=[0.05:0.1:0.95];
roo=repelem(r0,length(thetao));
R=repmat(roo,1,length(zo));
thetaoo=repmat(thetao,1,length(r0));
THETA=repmat(thetaoo,1,length(zo));
z=repelem(zo,1,length(r0)*length(thetao));
x=R.*cos(THETA);
y=R.*sin(THETA);
volume=(pi*0.03^2)/4*(R+0.005)/0.045*0.1*125;
 
phi=10^5*10^-3*sum(volume.*(- (28.*cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).^2)./(5.*(2.*a1.*x.*sin(z.*pi) + 2.*a2.*x.*sin(2.*z.*pi) + 2.*a3.*x.*sin(3.*z.*pi) + 2.*a4.*x.*sin(4.*z.*pi) + 2.*a5.*x.*sin(5.*z.*pi) + 2.*a6.*x.*sin(6.*z.*pi) - 1)) - (63.*(cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).^2./(2.*a1.*x.*sin(z.*pi) + 2.*a2.*x.*sin(2.*z.*pi) + 2.*a3.*x.*sin(3.*z.*pi) + 2.*a4.*x.*sin(4.*z.*pi) + 2.*a5.*x.*sin(5.*z.*pi) + 2.*a6.*x.*sin(6.*z.*pi) - 1) - (sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)) + sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1) + (pi.*cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./(a1.*sin(pi.*z) + a2.*sin(2.*pi.*z) + a3.*sin(3.*pi.*z) + a4.*sin(4.*pi.*z) + a5.*sin(5.*pi.*z) + a6.*sin(6.*pi.*z)).^2 + (x.*pi.*cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./((a1.*sin(z.*pi) + a2.*sin(2.*z.*pi) + a3.*sin(3.*z.*pi) + a4.*sin(4.*z.*pi) + a5.*sin(5.*z.*pi) + a6.*sin(6.*z.*pi)).*(1 - 2.*a2.*x.*sin(2.*pi.*z) - 2.*a3.*x.*sin(3.*pi.*z) - 2.*a4.*x.*sin(4.*pi.*z) - 2.*a5.*x.*sin(5.*pi.*z) - 2.*a6.*x.*sin(6.*pi.*z) - 2.*a1.*x.*sin(pi.*z)).^(1./2))).^2).^2)./20 - (28.*sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).^2)./(5.*(2.*a1.*x.*sin(z.*pi) + 2.*a2.*x.*sin(2.*z.*pi) + 2.*a3.*x.*sin(3.*z.*pi) + 2.*a4.*x.*sin(4.*z.*pi) + 2.*a5.*x.*sin(5.*z.*pi) + 2.*a6.*x.*sin(6.*z.*pi) - 1)) - (63.*((sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)) + sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1) + (pi.*cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./(a1.*sin(pi.*z) + a2.*sin(2.*pi.*z) + a3.*sin(3.*pi.*z) + a4.*sin(4.*pi.*z) + a5.*sin(5.*pi.*z) + a6.*sin(6.*pi.*z)).^2 + (x.*pi.*cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./((a1.*sin(z.*pi) + a2.*sin(2.*z.*pi) + a3.*sin(3.*z.*pi) + a4.*sin(4.*z.*pi) + a5.*sin(5.*z.*pi) + a6.*sin(6.*z.*pi)).*(1 - 2.*a2.*x.*sin(2.*pi.*z) - 2.*a3.*x.*sin(3.*pi.*z) - 2.*a4.*x.*sin(4.*pi.*z) - 2.*a5.*x.*sin(5.*pi.*z) - 2.*a6.*x.*sin(6.*pi.*z) - 2.*a1.*x.*sin(pi.*z)).^(1./2))).*(cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)) + cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1) - (pi.*sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./(a1.*sin(pi.*z) + a2.*sin(2.*pi.*z) + a3.*sin(3.*pi.*z) + a4.*sin(4.*pi.*z) + a5.*sin(5.*pi.*z) + a6.*sin(6.*pi.*z)).^2 - (x.*pi.*sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./((a1.*sin(z.*pi) + a2.*sin(2.*z.*pi) + a3.*sin(3.*z.*pi) + a4.*sin(4.*z.*pi) + a5.*sin(5.*z.*pi) + a6.*sin(6.*z.*pi)).*(1 - 2.*a2.*x.*sin(2.*pi.*z) - 2.*a3.*x.*sin(3.*pi.*z) - 2.*a4.*x.*sin(4.*pi.*z) - 2.*a5.*x.*sin(5.*pi.*z) - 2.*a6.*x.*sin(6.*pi.*z) - 2.*a1.*x.*sin(pi.*z)).^(1./2))) + (cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)))./(2.*a1.*x.*sin(z.*pi) + 2.*a2.*x.*sin(2.*z.*pi) + 2.*a3.*x.*sin(3.*z.*pi) + 2.*a4.*x.*sin(4.*z.*pi) + 2.*a5.*x.*sin(5.*z.*pi) + 2.*a6.*x.*sin(6.*z.*pi) - 1)).^2)./10 - (63.*(sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).^2./(2.*a1.*x.*sin(z.*pi) + 2.*a2.*x.*sin(2.*z.*pi) + 2.*a3.*x.*sin(3.*z.*pi) + 2.*a4.*x.*sin(4.*z.*pi) + 2.*a5.*x.*sin(5.*z.*pi) + 2.*a6.*x.*sin(6.*z.*pi) - 1) - (cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)) + cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1) - (pi.*sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./(a1.*sin(pi.*z) + a2.*sin(2.*pi.*z) + a3.*sin(3.*pi.*z) + a4.*sin(4.*pi.*z) + a5.*sin(5.*pi.*z) + a6.*sin(6.*pi.*z)).^2 - (x.*pi.*sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./((a1.*sin(z.*pi) + a2.*sin(2.*z.*pi) + a3.*sin(3.*z.*pi) + a4.*sin(4.*z.*pi) + a5.*sin(5.*z.*pi) + a6.*sin(6.*z.*pi)).*(1 - 2.*a2.*x.*sin(2.*pi.*z) - 2.*a3.*x.*sin(3.*pi.*z) - 2.*a4.*x.*sin(4.*pi.*z) - 2.*a5.*x.*sin(5.*pi.*z) - 2.*a6.*x.*sin(6.*pi.*z) - 2.*a1.*x.*sin(pi.*z)).^(1./2))).^2).^2)./20 + (28.*(sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)) + sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1) + (pi.*cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./(a1.*sin(pi.*z) + a2.*sin(2.*pi.*z) + a3.*sin(3.*pi.*z) + a4.*sin(4.*pi.*z) + a5.*sin(5.*pi.*z) + a6.*sin(6.*pi.*z)).^2 + (x.*pi.*cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./((a1.*sin(z.*pi) + a2.*sin(2.*z.*pi) + a3.*sin(3.*z.*pi) + a4.*sin(4.*z.*pi) + a5.*sin(5.*z.*pi) + a6.*sin(6.*z.*pi)).*(1 - 2.*a2.*x.*sin(2.*pi.*z) - 2.*a3.*x.*sin(3.*pi.*z) - 2.*a4.*x.*sin(4.*pi.*z) - 2.*a5.*x.*sin(5.*pi.*z) - 2.*a6.*x.*sin(6.*pi.*z) - 2.*a1.*x.*sin(pi.*z)).^(1./2))).^2)./5 + (28.*(cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)) + cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1) - (pi.*sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./(a1.*sin(pi.*z) + a2.*sin(2.*pi.*z) + a3.*sin(3.*pi.*z) + a4.*sin(4.*pi.*z) + a5.*sin(5.*pi.*z) + a6.*sin(6.*pi.*z)).^2 - (x.*pi.*sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./((a1.*sin(z.*pi) + a2.*sin(2.*z.*pi) + a3.*sin(3.*z.*pi) + a4.*sin(4.*z.*pi) + a5.*sin(5.*z.*pi) + a6.*sin(6.*z.*pi)).*(1 - 2.*a2.*x.*sin(2.*pi.*z) - 2.*a3.*x.*sin(3.*pi.*z) - 2.*a4.*x.*sin(4.*pi.*z) - 2.*a5.*x.*sin(5.*pi.*z) - 2.*a6.*x.*sin(6.*pi.*z) - 2.*a1.*x.*sin(pi.*z)).^(1./2))).^2)./5 + (63.*(- cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).^2./(2.*a1.*x.*sin(z.*pi) + 2.*a2.*x.*sin(2.*z.*pi) + 2.*a3.*x.*sin(3.*z.*pi) + 2.*a4.*x.*sin(4.*z.*pi) + 2.*a5.*x.*sin(5.*z.*pi) + 2.*a6.*x.*sin(6.*z.*pi) - 1) - sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).^2./(2.*a1.*x.*sin(z.*pi) + 2.*a2.*x.*sin(2.*z.*pi) + 2.*a3.*x.*sin(3.*z.*pi) + 2.*a4.*x.*sin(4.*z.*pi) + 2.*a5.*x.*sin(5.*z.*pi) + 2.*a6.*x.*sin(6.*z.*pi) - 1) + (sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)) + sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1) + (pi.*cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./(a1.*sin(pi.*z) + a2.*sin(2.*pi.*z) + a3.*sin(3.*pi.*z) + a4.*sin(4.*pi.*z) + a5.*sin(5.*pi.*z) + a6.*sin(6.*pi.*z)).^2 + (x.*pi.*cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./((a1.*sin(z.*pi) + a2.*sin(2.*z.*pi) + a3.*sin(3.*z.*pi) + a4.*sin(4.*z.*pi) + a5.*sin(5.*z.*pi) + a6.*sin(6.*z.*pi)).*(1 - 2.*a2.*x.*sin(2.*pi.*z) - 2.*a3.*x.*sin(3.*pi.*z) - 2.*a4.*x.*sin(4.*pi.*z) - 2.*a5.*x.*sin(5.*pi.*z) - 2.*a6.*x.*sin(6.*pi.*z) - 2.*a1.*x.*sin(pi.*z)).^(1./2))).^2 + (cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)) + cos((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1) - (pi.*sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*((- 2.*a1.*x.*sin(z.*pi) - 2.*a2.*x.*sin(2.*z.*pi) - 2.*a3.*x.*sin(3.*z.*pi) - 2.*a4.*x.*sin(4.*z.*pi) - 2.*a5.*x.*sin(5.*z.*pi) - 2.*a6.*x.*sin(6.*z.*pi) + 1).^(1./2) - 1).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./(a1.*sin(pi.*z) + a2.*sin(2.*pi.*z) + a3.*sin(3.*pi.*z) + a4.*sin(4.*pi.*z) + a5.*sin(5.*pi.*z) + a6.*sin(6.*pi.*z)).^2 - (x.*pi.*sin((60.*a1 + 30.*a2 + 20.*a3 + 15.*a4 + 12.*a5 + 10.*a6 - 60.*a1.*cos(z.*pi) - 30.*a2.*cos(2.*z.*pi) - 20.*a3.*cos(3.*z.*pi) - 15.*a4.*cos(4.*z.*pi) - 12.*a5.*cos(5.*z.*pi) - 10.*a6.*cos(6.*z.*pi))./(60.*pi)).*(a1.*cos(z.*pi) + 2.*a2.*cos(2.*z.*pi) + 3.*a3.*cos(3.*z.*pi) + 4.*a4.*cos(4.*z.*pi) + 5.*a5.*cos(5.*z.*pi) + 6.*a6.*cos(6.*z.*pi)))./((a1.*sin(z.*pi) + a2.*sin(2.*z.*pi) + a3.*sin(3.*z.*pi) + a4.*sin(4.*z.*pi) + a5.*sin(5.*z.*pi) + a6.*sin(6.*z.*pi)).*(1 - 2.*a2.*x.*sin(2.*pi.*z) - 2.*a3.*x.*sin(3.*pi.*z) - 2.*a4.*x.*sin(4.*pi.*z) - 2.*a5.*x.*sin(5.*pi.*z) - 2.*a6.*x.*sin(6.*pi.*z) - 2.*a1.*x.*sin(pi.*z)).^(1./2))).^2 + 1).^2)./20 - 133./4),2);
L=1;
gp=0;
forward_y = @(xx) cos(a(1)*L/pi*(1-cos(pi*xx/L))+a(2)*L/(2*pi)*(1-cos(2*pi*xx/L))+a(3)*L/(3*pi)*(1-cos(3*pi*xx/L))+a(4)*L/(4*pi)*(1-cos(4*pi*xx/L))+a(5)*L/(5*pi)*(1-cos(5*pi*xx/L))+a(6)*L/(6*pi)*(1-cos(6*pi*xx/L)));
for i=1:10
    zoo=i*0.1;
    znew=integral(@(xx) forward_y(xx),0,zoo);
    gp=gp+znew*5*0.15^2*pi*1.07*9.81;
end

energy=phi+gp;% total energy = strain energy phi + gravitational potential energy gp
end