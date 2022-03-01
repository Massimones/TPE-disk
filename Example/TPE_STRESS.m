function [x,tau11s,tau22s,tau33s,tau13s]=TPE_STRESS(H,alfa,dp,dT,a,db,ni,mu,lambda,c,limiteplot,zlm,k)

%%%The function generates stresses for a 1 disk TPE inclusion 
%%% Input parameters are in S.I.
%%% Massimo Nespoli 01/03/2022

%%%E.g.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H=10*10^9;        % Constant of Biot
% alfa=3*10^(-5);   % thermal expansion
% dp=1e6;           % pore pressure change of disk 1
% dT=100;           % Temperature change of disk 1
% dpB=0.5e6;        % pore pressure change of disk 2
% dTB=0;            % Temperature change of disk 1
% a=2000;           % disk radius
% db=100;           % disk height 
% ni=0.2;           % Poisson modulus
% mu=6*10^9;        % Shear modulus
% lambda=4*10^9;    % Lamè constant
% MedianPlane=3000; % TPE inclusion depth   of median plane  
% limiteplot=4000;  % Limit in plot (max(x))
% k=100;            % step for plot in x
% Zmin=MedianPlane-2000;   % min z for computation
% Zmax=MedianPlane+100;    % max z for computation
% Zstep=100;        % step for plot in z
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

high=db; 
e0=(1/(3*H))*dp+(1/3)*alfa*dT;
e1=e0*(1+ni)/(1-ni);
A=(e1*high)/(2*a);
%DEFINING RECEIVER
x=0.1:k:limiteplot;
y=0;
z=c-zlm;
r=sqrt(x.^2+y.^2+z.^2);
ct=z./r;
st=sqrt(x.^2+y.^2)./r;
sp=y./sqrt(x.^2+y.^2);
cp=x./sqrt(x.^2+y.^2);
errl(length(r))=0;
err(length(r))=0;
ettl(length(r))=0;
ett(length(r))=0;
ertl(length(r))=0;
ert(length(r))=0;
effl(length(r))=0;
eff(length(r))=0;
M=1000;
L=2*M;
cc(L)=0;
cr(L)=0; 
tau11s(length(r))=0;
tau22s(length(r))=0;
tau33s(length(r))=0;
tau13s(length(r))=0;

P(L)=0;
pp(L)=0;
ppc(L)=0;
ppp(L)=0;
 

for i=1:length(r)
for m=2:L
P(1,i)=1;
P(2,i)=ct(i);
P(m+1,i)=((2*m-1)*ct(i)*P(m,i)-(m-1)*P(m-1,i))/m; % polinomio di Legendre
pp_zero(i)=0;
ppc(m,i)=(m*ct(i)*P(m+1,i)-(m)*P(m,i))/(1-ct(i)^2);
pp(m,i)=ppc(m,i)*st(i); %First derivative of  Legendre p
ppp(m,i)=ct(i)*pp(m,i)/(-st(i))-m*(m+1)*P(m+1,i);  %Second derivative of  Legendre p
end
for m=1:M
l=2*m;
if m<80
            cc(l)=(((-1)^m)/4^m)*factorial(2*m)/(factorial(m))^2;
         else
             cc(l)=((-1)^m)/(sqrt(pi*m));
 end
cr(l)=cc(l);
cr(1)=1;
%INTERNAL DOMAIN
if abs(r(i))<a 
errl(l,i)=l*cc(l).*P(l+1,i).*(abs(r(i))./a).^(l-2);
ettl(l,i)=(cc(l)/(l-1)).*(ppp(l,i)+l.*P(l+1,i)).*(abs(r(i))./a).^(l-2);
effl(l,i)=(cc(l)/(l-1)).*((ct(i)./st(i)).*pp(l,i)+l.*P(l+1,i)).*(abs(r(i))./a).^(l-2);
ertl(l,i)=cc(l).*pp(l,i).*(abs(r(i))./a).^(l-2);
ert(i)=A*sum(ertl(:,i));
% if abs(z)==high/2
% ett(i)=A*(sum(ettl(:,i)))+2*A*a/high;
% else
ett(i)=A*(sum(ettl(:,i)));
% end
err(i)=A*sum(errl(:,i));
eff(i)=A*sum(effl(:,i));
% %CONFINE
% elseif abs(r(i))==a
% err(i)=err(i-1);
% ett(i)=ett(i-1);
% eff(i)=eff(i-1);
% ert(i)=ert(i-1);
%EXTERNAL DOMAIN
elseif abs(r(i))>a
errl(1,i)=(1)*cr(1).*P(1,i).*(a./abs(r(i))).^(3);
effl(1,i)=(cr(1)./(2)).*((-ct(i)./st(i)).*pp_zero(i)+(1).*P(1,i)).*(a./abs(r(i))).^(3); %il meno va nel coseno
errl(l,i)=(l+1)*cr(l).*P(l+1,i).*(a./abs(r(i))).^(l+3);
ertl(l,i)=cc(l).*pp(l,i).*(a./abs(r(i))).^(l+3);
ettl(l,i)=(cc(l)/(m+1)).*(ppp(l,i)-(l+1).*P(l+1,i)).*(a./abs(r(i))).^l;
effl(l,i)=(cc(l)./(l+2)).*((-ct(i)./st(i)).*pp(l,i)+(l+1).*P(l+1,i)).*(a./abs(r(i))).^(l+3); %il meno va nel coseno
ett(i)=(A/2)*(a./abs(r(i))).^3.*(1-sum(ettl(:,i)));
ert(i)=A*sum(ertl(:,i));
err(i)=-A*sum(errl(:,i));
eff(i)=A*sum(effl(:,i));
end
end
end



for ju=1:length(err)
   matrixR(1,1)=st(ju)*cp(ju);
   matrixR(1,2)=ct(ju)*cp(ju);
   matrixR(1,3)=-sp(ju);
   matrixR(2,1)=st(ju)*sp(ju);
   matrixR(2,2)=ct(ju)*sp(ju);
   matrixR(2,3)=cp(ju);
   matrixR(3,1)=ct(ju);
   matrixR(3,2)=-st(ju);
   matrixR(3,3)=0;
  
  e_polar=[err(ju),0,ert(ju);0,eff(ju),0;ert(ju),0,ett(ju)];
  e_cartesian=matrixR'*e_polar*matrixR;  
  e11smio(ju)=e_cartesian(1,1);
  e12smio(ju)=e_cartesian(1,3);
  e13smio(ju)=e_cartesian(1,2);
  e21smio(ju)=e_cartesian(3,1);
  e22smio(ju)=e_cartesian(3,3);
  e23smio(ju)=e_cartesian(3,2);
  e31smio(ju)=e_cartesian(2,1);
  e32smio(ju)=e_cartesian(2,3);
  e33smio(ju)=e_cartesian(2,2); 
end


e11s=err.*(st.^2).*(cp.^2)+2.*ert.*st.*ct.*(cp.^2)+ett.*(ct.^2).*(cp.^2)+eff.*(sp.^2);
e22s=err.*(st.^2).*(sp.^2)+2.*ert.*st.*ct.*(sp.^2)+ett.*(ct.^2).*(sp.^2)+eff.*(cp.^2);
e33s=err.*(ct.^2)-2*ert.*st.*ct+ett.*(st.^2);
   if abs(z)<=high/2
 e33s(x<a)=e33s(x<a)+2*A*a/high;
   end

e13s=e13smio;

for i=1:length(x)
  if abs(z)<=high/2 

 if x(i)<=a
 tau11s(i)=2*mu.*(-e1+e11s(i));
 tau22s(i)=2*mu.*(-e1+e22s(i));
 tau33s(i)=2*mu.*(-e1+e33s(i));
 tau13s(i)=2*mu.*(e13s(i));


elseif x(i)>a
tau11s(i)=2*mu*e11s(i);
 tau22s(i)=2*mu*e22s(i);
 tau33s(i)=2*mu*e33s(i);
 tau13s(i)=2*mu.*(e13s(i));


 end
  else
  %TAU OUT    
 tau11s(i)=2*mu*e11s(i);
 tau22s(i)=2*mu*e22s(i);
 tau33s(i)=2*mu*e33s(i);
 tau13s(i)=2*mu.*(e13s(i));

 end
end

end
