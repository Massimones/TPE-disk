clear all

%%%The function plot stresses of a TPE inclusion
%%% Input parameters are in S.I.
%%% Massimo Nespoli 01/03/2022
%%%%input%%%%%%
load CaseTEST
xlim=2000;
zlim1=0;
zlim2=100;
%%%%%%%%%%%%%%%%


 xr = reshape(x', 1, []);
 zr = 3000-reshape(z', 1, []);
 tau11r = reshape(tau11', 1, []);
 tau22r=reshape(tau22', 1, []);
 tau33r=reshape(tau33', 1, []);
 tau13r=reshape(tau13', 1, []);
 
 
figure('Position', [10 10 1300 200])
subplot(1,4,1)
xlin2=linspace(min(xr),max(xr),500);
ylin2=linspace(min(zr),max(zr),500);
[X,Y]=meshgrid(xlin2,ylin2);
disp =griddata(xr,zr,tau11r,X,Y); 
L=image(xlin2,ylin2,disp,'Cdatamapping','scaled');
 colorbar
      load('MyColormaps.mat','bluered');
set(gcf,'Colormap',bluered);
hold on
pos=[0 -100 2000 200];
rectangle('Position',pos,'LineWidth',5)

      caxis([-1 1]*1e6);
     set(gca, 'YDir', 'normal');
     
     %axis equal 
     axis([0 2*a 0-200/2 2*200]);
     xlabel('x (m)');
     ylabel('z (m)');
     title('\tau_1_1 (Pa)');

subplot(1,4,2)
xlin2=linspace(min(xr),max(xr),500);
ylin2=linspace(min(zr),max(zr),500);
[X,Y]=meshgrid(xlin2,ylin2);
disp =griddata(xr,zr,tau22r,X,Y); 
L=image(xlin2,ylin2,disp,'Cdatamapping','scaled');
 colorbar
      load('MyColormaps.mat','bluered');
set(gcf,'Colormap',bluered);
hold on
pos=[0 -100 2000 200];
rectangle('Position',pos,'LineWidth',5)

       caxis([-1 1]*1e6);
     set(gca, 'YDir', 'normal');
     
     %axis equal 
     axis([0 2*a 0-200/2 2*200]);
     xlabel('x (m)');
     ylabel('z (m)');
     title('\tau_2_2 (Pa)');

     subplot(1,4,3)
xlin2=linspace(min(xr),max(xr),500);
ylin2=linspace(min(zr),max(zr),500);
[X,Y]=meshgrid(xlin2,ylin2);
disp =griddata(xr,zr,tau33r,X,Y); 
L=image(xlin2,ylin2,disp,'Cdatamapping','scaled');
 colorbar
      load('MyColormaps.mat','bluered');
set(gcf,'Colormap',bluered);
hold on
pos=[0 -100 2000 200];
rectangle('Position',pos,'LineWidth',5)

       caxis([-1 1]*1e6);
     set(gca, 'YDir', 'normal');
     
     %axis equal 
     axis([0 2*a 0-200/2 2*200]);
     xlabel('x (m)');
     ylabel('z (m)');
     title('\tau_3_3 (Pa)');
     
     subplot(1,4,4)
xlin2=linspace(min(xr),max(xr),500);
ylin2=linspace(min(zr),max(zr),500);
[X,Y]=meshgrid(xlin2,ylin2);
disp =griddata(xr,zr,tau13r,X,Y); 
L=image(xlin2,ylin2,disp,'Cdatamapping','scaled');
 colorbar
      load('MyColormaps.mat','bluered');
set(gcf,'Colormap',bluered);
hold on
pos=[0 -100 2000 200];
rectangle('Position',pos,'LineWidth',5)

       caxis([-1 1]*1e6);
     set(gca, 'YDir', 'normal');
     
     %axis equal 
     axis([0 2*a 0-200/2 2*200]);
     xlabel('x (m)');
     ylabel('z (m)');
     title('\tau_1_3 (Pa)');

 

