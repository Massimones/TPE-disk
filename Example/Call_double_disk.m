clear all

%%%The script generates stresses for a TPE inclusion 
%%% composed by two disks with different p and T
%%% Input parameters are in S.I.
%%% Massimo Nespoli 01/03/2022
%E.g.  2disks_dT=100K_dp=1MPa_dT=0K_dp=05MPa
%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=10*10^9;               % Constant of Biot
alfa=3*10^(-5);          % thermal expansion
dp=1e6;                  % pore pressure change of disk 1
dT=100;                  % Temperature change of disk 1
dpB=0.5e6;               % pore pressure change of disk 2
dTB=0;                   % Temperature change of disk 1
a=2000;                  % disk radius
db=100;                  % disk height 
ni=0.2;                  % Poisson modulus
mu=6*10^9;               % Shear modulus
lambda=4*10^9;           % Lamè constant
MedianPlane=3000;        % TPE inclusion, depth   of median plane  
limiteplot=4000;         % Limit in plot (max(x))
k=100;                   % step for plot in x
Zmin=MedianPlane-400;   % min z for computation
Zmax=MedianPlane+100;    % max z for computation
Zstep=10;               % step for plot in z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=MedianPlane+db/2; 
cB=MedianPlane-db/2;
Zetav=Zmin:Zstep:Zmax;

for i=1:length(Zetav)
    disp(i)
    zlm=Zetav(i);
[xA(i,:),tau11A(i,:),tau22A(i,:),tau33A(i,:),tau13A(i,:)]=TPE_STRESS(H,alfa,dp,dT,a,db,ni,mu,lambda,c,limiteplot,zlm,k);

[xB(i,:),tau11B(i,:),tau22B(i,:),tau33B(i,:),tau13B(i,:)]=TPE_STRESS(H,alfa,dpB,dTB,a,db,ni,mu,lambda,cB,limiteplot,zlm,k);
x(i,:)=xA(i,:);
tau11(i,:)=tau11A(i,:)+tau11B(i,:);
tau22(i,:)=tau22A(i,:)+tau22B(i,:);
tau33(i,:)=tau33A(i,:)+tau33B(i,:);
tau13(i,:)=tau13A(i,:)+tau13B(i,:);

z(i,1:length(x(1,:)))=zlm;
end

save CaseTEST
