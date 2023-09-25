clear
w=5; %number of modes
for i=1:200
    p=0.1*rand(w,1);%mode parameters. 
    a1=p(1);
    a2=p(2);
    a3=p(3);
    a4=p(4);
    a5=p(5);
    L=10;%height of chamber
    r0=[1.6 1.8];%radius of inter and outer sidewall.
    %choose sample points during the numerical integration
    thetao=0:pi/2:3/2*pi;
    zo=0:1:10;
    roo=repelem(r0,length(thetao));
    R=repmat(roo,1,length(zo));
    thetaoo=repmat(thetao,1,length(r0));
    THETA=repmat(thetaoo,1,length(zo));
    z=repelem(zo,1,length(r0)*length(thetao));
    x=R.*cos(THETA);
    y=R.*sin(THETA);
    %strain energy (numerical integration of strain energy density over the body. Strain energy density can be compute from StrainEnergy_Analytical.m
    phi(i)=sum((119.*((x.*y)./((x.^2 + y.^2).^(1./2).*(a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L) + x.^2 + y.^2).^(1./2)) - (x.*y.*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2))./(x.^2 + y.^2).^(3./2)).^2)./50 - (5757.*(((x.*y)./((x.^2 + y.^2).^(1./2).*(a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L) + x.^2 + y.^2).^(1./2)) - (x.*y.*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2))./(x.^2 + y.^2).^(3./2)).*((x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2)./(x.^2 + y.^2).^(1./2) + x.^2./((x.^2 + y.^2).^(1./2).*(a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L) + x.^2 + y.^2).^(1./2)) - (x.^2.*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2))./(x.^2 + y.^2).^(3./2)) + ((x.*y)./((x.^2 + y.^2).^(1./2).*(a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L) + x.^2 + y.^2).^(1./2)) - (x.*y.*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2))./(x.^2 + y.^2).^(3./2)).*((x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2)./(x.^2 + y.^2).^(1./2) + y.^2./((x.^2 + y.^2).^(1./2).*(a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L) + x.^2 + y.^2).^(1./2)) - (y.^2.*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2))./(x.^2 + y.^2).^(3./2)) + (x.*y.*((a1.*pi.*cos((pi.*z)./L))./L + (3.*a2.*pi.*cos((3.*pi.*z)./L))./L + (5.*a3.*pi.*cos((5.*pi.*z)./L))./L + (7.*a4.*pi.*cos((7.*pi.*z)./L))./L + (9.*a5.*pi.*cos((9.*pi.*z)./L))./L).^2)./(4.*(x.^2 + y.^2).*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)))).^2)./250 + (5757.*(2.*((x.*y)./((x.^2 + y.^2).^(1./2).*(a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L) + x.^2 + y.^2).^(1./2)) - (x.*y.*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2))./(x.^2 + y.^2).^(3./2)).^2 + ((x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2)./(x.^2 + y.^2).^(1./2) + x.^2./((x.^2 + y.^2).^(1./2).*(a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L) + x.^2 + y.^2).^(1./2)) - (x.^2.*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2))./(x.^2 + y.^2).^(3./2)).^2 + ((x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2)./(x.^2 + y.^2).^(1./2) + y.^2./((x.^2 + y.^2).^(1./2).*(a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L) + x.^2 + y.^2).^(1./2)) - (y.^2.*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2))./(x.^2 + y.^2).^(3./2)).^2 + (x.^2.*((a1.*pi.*cos((pi.*z)./L))./L + (3.*a2.*pi.*cos((3.*pi.*z)./L))./L + (5.*a3.*pi.*cos((5.*pi.*z)./L))./L + (7.*a4.*pi.*cos((7.*pi.*z)./L))./L + (9.*a5.*pi.*cos((9.*pi.*z)./L))./L).^2)./(4.*(x.^2 + y.^2).*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L))) + (y.^2.*((a1.*pi.*cos((pi.*z)./L))./L + (3.*a2.*pi.*cos((3.*pi.*z)./L))./L + (5.*a3.*pi.*cos((5.*pi.*z)./L))./L + (7.*a4.*pi.*cos((7.*pi.*z)./L))./L + (9.*a5.*pi.*cos((9.*pi.*z)./L))./L).^2)./(4.*(x.^2 + y.^2).*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L))) + 1).^2)./500 - (5757.*(((x.*y)./((x.^2 + y.^2).^(1./2).*(a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L) + x.^2 + y.^2).^(1./2)) - (x.*y.*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2))./(x.^2 + y.^2).^(3./2)).^2 + ((x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2)./(x.^2 + y.^2).^(1./2) + x.^2./((x.^2 + y.^2).^(1./2).*(a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L) + x.^2 + y.^2).^(1./2)) - (x.^2.*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2))./(x.^2 + y.^2).^(3./2)).^2 + (x.^2.*((a1.*pi.*cos((pi.*z)./L))./L + (3.*a2.*pi.*cos((3.*pi.*z)./L))./L + (5.*a3.*pi.*cos((5.*pi.*z)./L))./L + (7.*a4.*pi.*cos((7.*pi.*z)./L))./L + (9.*a5.*pi.*cos((9.*pi.*z)./L))./L).^2)./(4.*(x.^2 + y.^2).*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)))).^2)./500 - (5757.*(((x.*y)./((x.^2 + y.^2).^(1./2).*(a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L) + x.^2 + y.^2).^(1./2)) - (x.*y.*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2))./(x.^2 + y.^2).^(3./2)).^2 + ((x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2)./(x.^2 + y.^2).^(1./2) + y.^2./((x.^2 + y.^2).^(1./2).*(a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L) + x.^2 + y.^2).^(1./2)) - (y.^2.*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2))./(x.^2 + y.^2).^(3./2)).^2 + (y.^2.*((a1.*pi.*cos((pi.*z)./L))./L + (3.*a2.*pi.*cos((3.*pi.*z)./L))./L + (5.*a3.*pi.*cos((5.*pi.*z)./L))./L + (7.*a4.*pi.*cos((7.*pi.*z)./L))./L + (9.*a5.*pi.*cos((9.*pi.*z)./L))./L).^2)./(4.*(x.^2 + y.^2).*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)))).^2)./500 + (119.*((x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2)./(x.^2 + y.^2).^(1./2) + x.^2./((x.^2 + y.^2).^(1./2).*(a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L) + x.^2 + y.^2).^(1./2)) - (x.^2.*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2))./(x.^2 + y.^2).^(3./2)).^2)./100 + (119.*((x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2)./(x.^2 + y.^2).^(1./2) + y.^2./((x.^2 + y.^2).^(1./2).*(a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L) + x.^2 + y.^2).^(1./2)) - (y.^2.*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L)).^(1./2))./(x.^2 + y.^2).^(3./2)).^2)./100 - (10919.*x.^2.*((a1.*pi.*cos((pi.*z)./L))./L + (3.*a2.*pi.*cos((3.*pi.*z)./L))./L + (5.*a3.*pi.*cos((5.*pi.*z)./L))./L + (7.*a4.*pi.*cos((7.*pi.*z)./L))./L + (9.*a5.*pi.*cos((9.*pi.*z)./L))./L).^2)./(2000.*(x.^2 + y.^2).*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L))) - (10919.*y.^2.*((a1.*pi.*cos((pi.*z)./L))./L + (3.*a2.*pi.*cos((3.*pi.*z)./L))./L + (5.*a3.*pi.*cos((5.*pi.*z)./L))./L + (7.*a4.*pi.*cos((7.*pi.*z)./L))./L + (9.*a5.*pi.*cos((9.*pi.*z)./L))./L).^2)./(2000.*(x.^2 + y.^2).*(x.^2 + y.^2 + a1.*sin((z.*pi)./L) + a2.*sin((3.*z.*pi)./L) + a3.*sin((5.*z.*pi)./L) + a4.*sin((7.*z.*pi)./L) + a5.*sin((9.*z.*pi)./L))) - 41489./500,2);
    start=1;
    tail=w;
    for count=1:w
        for k=start:tail
            M=zeros(w);
            M(count,k-start+count)=1;
            M(k-start+count,count)=1;
            A(i,k)=p'*M*p;
        end
        start=tail+1;
        tail=start+w-count-1;
    end
end
W1=pinv(A)*phi';% W1 is a vector contains all the entries of the weight matrix needed.
 
%get weight matrix in pseudoinverse.
W=zeros(w);
start=1;
tail=w;
for count=1:w
    for k=start:tail
        M=zeros(w);
        M(count,k-start+count)=1;
        M(k-start+count,count)=1;
        W=M*W1(k)+W;
    end
    start=tail+1;
    tail=start+w-count-1;
end
W=W/norm(W);  
