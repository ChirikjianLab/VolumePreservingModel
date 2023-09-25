% 2D rod, projected gradient approach.
clear
tic
xf_desired=0.8; %desired x coordinate
yf_desired=0.2; %desired y coordinate
orif_desired=0.5; %desired orientation
n=6; %number of modes
length=1; %arclength of backbone curve
a=0.001*[1*(rand-0.5);1*(rand-0.5);1*(rand-0.5);1*(rand-0.5);1*(rand-0.5);1*(rand-0.5)]; %initial guess of mode parameters (close to zero)
max=50; %number of iterations

forward_x = @(a1,a2,a3,a4,a5,a6,x) cos(a1*length/pi*(1-cos(pi*x/length))+a2*length/(2*pi)*(1-cos(2*pi*x/length))+a3*length/(3*pi)*(1-cos(3*pi*x/length))+a4*length/(4*pi)*(1-cos(4*pi*x/length))+a5*length/(5*pi)*(1-cos(5*pi*x/length))+a6*length/(6*pi)*(1-cos(6*pi*x/length)));
forward_y = @(a1,a2,a3,a4,a5,a6,x) sin(a1*length/pi*(1-cos(pi*x/length))+a2*length/(2*pi)*(1-cos(2*pi*x/length))+a3*length/(3*pi)*(1-cos(3*pi*x/length))+a4*length/(4*pi)*(1-cos(4*pi*x/length))+a5*length/(5*pi)*(1-cos(5*pi*x/length))+a6*length/(6*pi)*(1-cos(6*pi*x/length)));

for iteration=1:max
Jacobian=jacob(n,a);
[forwardkin] = forward(a);%get the coordinates of the end of backbone curve
x_current=forwardkin(1);
y_current=forwardkin(2);
ori_current=a(1)*length/pi*(1-cos(pi*length/length))+a(2)*length/(2*pi)*(1-cos(2*pi*length/length))+a(3)*length/(3*pi)*(1-cos(3*pi*length/length))+a(4)*length/(4*pi)*(1-cos(4*pi*length/length))+a(5)*length/(5*pi)*(1-cos(5*pi*length/length))+a(6)*length/(6*pi)*(1-cos(6*pi*length/length));

%get desired end point position for current iteration
x_desired(iteration)=x_current+(xf_desired-x_current)/(max-iteration+1);
y_desired(iteration)=y_current+(yf_desired-y_current)/(max-iteration+1);
ori_desired(iteration)=ori_current+(orif_desired-ori_current)/(max-iteration+1);

lambdas=0.00001; %a very small number to avoid singularity when doing matrix inverse
invjacob=Jacobian'*inv(Jacobian*Jacobian'+lambdas*eye(3));
Z=eye(n)-invjacob*Jacobian;
    
c=[x_desired(iteration)-x_current;y_desired(iteration)-y_current;ori_desired(iteration)-ori_current];
middle=invjacob*c;
alpha=0.01;
last=alpha*Z*[a(1);...
    a(2);...
    a(3);...
    a(4);...
    a(5);...
    a(6)]*length;

a(:)=a(:)+middle-last;
end
toc

% plot the result
i=0;
for s=0:0.001:length
    i=i+1;
    y1(i)=integral(@(x) forward_x(a(1),a(2),a(3),a(4),a(5),a(6),x),0,s);
    y2(i)=integral(@(x) forward_y(a(1),a(2),a(3),a(4),a(5),a(6),x),0,s);            
end
plot(1*y2,1*y1,'b','Linewidth',1.0,'color',[0 0.6 0])
axis([-1 1 -1 1])
hold on
scatter(1*yf_desired,1*xf_desired,'r','filled');
hold off




function [forwardkin] = forward(a)
% this function computes the end point coordinates for mode parameters a
length=1;
forward_x = @(a,x) cos(a(1)*length/pi*(1-cos(pi*x/length))+a(2)*length/(2*pi)*(1-cos(2*pi*x/length))+a(3)*length/(3*pi)*(1-cos(3*pi*x/length))+a(4)*length/(4*pi)*(1-cos(4*pi*x/length))+a(5)*length/(5*pi)*(1-cos(5*pi*x/length))+a(6)*length/(6*pi)*(1-cos(6*pi*x/length)));
forward_y = @(a,x) sin(a(1)*length/pi*(1-cos(pi*x/length))+a(2)*length/(2*pi)*(1-cos(2*pi*x/length))+a(3)*length/(3*pi)*(1-cos(3*pi*x/length))+a(4)*length/(4*pi)*(1-cos(4*pi*x/length))+a(5)*length/(5*pi)*(1-cos(5*pi*x/length))+a(6)*length/(6*pi)*(1-cos(6*pi*x/length)));
forwardkin=[integral(@(x) forward_x(a,x),0,length);integral(@(x) forward_y(a,x),0,length)];
end

function [Jacobian] = jacob(n,a)
% this function computes the jacobian matrix for mode parameters a
length=1;
epsilon=10^-5;
for count=1:n
    position=zeros(n,1);
    position(count)=1;
    forwardkin1 = forward(a+0.5*epsilon*position);
    forwardkin0 = forward(a-0.5*epsilon*position);
    Jacobian1(:,count)=1/epsilon*(forwardkin1-forwardkin0);
end
    Jacobian2=[1*length/pi*(1-cos(pi*1/length)),1*length/(2*pi)*(1-cos(2*pi*1/length)),1*length/(3*pi)*(1-cos(3*pi*1/length)),1*length/(4*pi)*(1-cos(4*pi*1/length)),1*length/(5*pi)*(1-cos(5*pi*1/length)),1*length/(6*pi)*(1-cos(6*pi*1/length))];
    Jacobian=[Jacobian1;Jacobian2];
end