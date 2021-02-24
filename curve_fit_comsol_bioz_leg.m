clc
clear
close all

wf = [4e-3 6e-3 8e-3 10e-3 12e-3];
wm = [50e-3 52e-3 54e-3 56e-3 58e-3];
wf0 = 4e-3;
wm0 = 50e-3;
lb_4cm = [4 0 0 0.98 0 0 0.64 0 0 0.83];
ub_4cm = [5 200 1e-3 0.99 200 1e-3 0.65 200 1e-3 0.84];
lb_2cm = [4 0 0 0.98 0 3e-5 0.64 0 0 0.83];
ub_2cm = [5 200 1.5e-4 0.99 200 7e-5 0.65 70 1e-3 0.84];
param0 = [10 31 1/7e3 0.99 16 1/17e3 0.71 35 1/1.5e6 0.8];
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxFunctionEvaluations',6000,'StepTolerance',1e-6,'MaxIterations',1000)

z_comsol_leg_wf = csvread('z_comsol_leg_wf.csv',5);
z_comsol_leg_wm = csvread('z_comsol_leg_wm.csv',5);

freq = z_comsol_leg_wf(:,1);
omega = 2*pi*freq;

fun_abs = @(param,freq)abs(param(1) + param(2)./(1+param(2)*param(3)*(j*omega).^param(4)) + param(5)./(1+param(5)*param(6)*(j*omega).^param(7)) + param(8)./(1+param(8)*param(9)*(j*omega).^param(10)));

% a_comsol_cad_2cm = csvread('mag_comsol_leg_2cm.csv',1,1);
% p_comsol_cad_2cm = csvread('pha_comsol_leg_2cm.csv',1,1);
% a_comsol_cad_4cm = csvread('mag_comsol_leg_4cm.csv',1,1);
% p_comsol_cad_4cm = csvread('pha_comsol_leg_4cm.csv',1,1);

Z_2cm_wf1 = z_comsol_leg_wf(:,2);
Z_2cm_wf2 = z_comsol_leg_wf(:,5);
Z_2cm_wf3 = z_comsol_leg_wf(:,8);
Z_2cm_wf4 = z_comsol_leg_wf(:,11);
Z_2cm_wf5 = z_comsol_leg_wf(:,14);

Z_4cm_wf1 = z_comsol_leg_wf(:,17);
Z_4cm_wf2 = z_comsol_leg_wf(:,20);
Z_4cm_wf3 = z_comsol_leg_wf(:,23);
Z_4cm_wf4 = z_comsol_leg_wf(:,26);
Z_4cm_wf5 = z_comsol_leg_wf(:,29);

param_fit_abs_4cm_wf1 = lsqcurvefit(fun_abs,param0,omega,abs(Z_4cm_wf1),lb_4cm,ub_4cm,options);
param_fit_abs_4cm_wf2 = lsqcurvefit(fun_abs,param0,omega,abs(Z_4cm_wf2),lb_4cm,ub_4cm,options);
param_fit_abs_4cm_wf3 = lsqcurvefit(fun_abs,param0,omega,abs(Z_4cm_wf3),lb_4cm,ub_4cm,options);
param_fit_abs_4cm_wf4 = lsqcurvefit(fun_abs,param0,omega,abs(Z_4cm_wf4),lb_4cm,ub_4cm,options);
param_fit_abs_4cm_wf5 = lsqcurvefit(fun_abs,param0,omega,abs(Z_4cm_wf5),lb_4cm,ub_4cm,options);

param_fit_abs_2cm_wf1 = lsqcurvefit(fun_abs,param0,omega,abs(Z_2cm_wf1),lb_2cm,ub_2cm,options);
param_fit_abs_2cm_wf2 = lsqcurvefit(fun_abs,param0,omega,abs(Z_2cm_wf2),lb_2cm,ub_2cm,options);
param_fit_abs_2cm_wf3 = lsqcurvefit(fun_abs,param0,omega,abs(Z_2cm_wf3),lb_2cm,ub_2cm,options);
param_fit_abs_2cm_wf4 = lsqcurvefit(fun_abs,param0,omega,abs(Z_2cm_wf4),lb_2cm,ub_2cm,options);
param_fit_abs_2cm_wf5 = lsqcurvefit(fun_abs,param0,omega,abs(Z_2cm_wf5),lb_2cm,ub_2cm,options);

Z_fit_abs_4cm_wf1 = param_fit_abs_4cm_wf1(1) + param_fit_abs_4cm_wf1(2)./(1+param_fit_abs_4cm_wf1(2)*param_fit_abs_4cm_wf1(3)*(j*omega).^param_fit_abs_4cm_wf1(4)) + param_fit_abs_4cm_wf1(5)./(1+param_fit_abs_4cm_wf1(5)*param_fit_abs_4cm_wf1(6)*(j*omega).^param_fit_abs_4cm_wf1(7)) + param_fit_abs_4cm_wf1(8)./(1+param_fit_abs_4cm_wf1(8)*param_fit_abs_4cm_wf1(9)*(j*omega).^param_fit_abs_4cm_wf1(10));
Z_fit_abs_4cm_wf2 = param_fit_abs_4cm_wf2(1) + param_fit_abs_4cm_wf2(2)./(1+param_fit_abs_4cm_wf2(2)*param_fit_abs_4cm_wf2(3)*(j*omega).^param_fit_abs_4cm_wf2(4)) + param_fit_abs_4cm_wf2(5)./(1+param_fit_abs_4cm_wf2(5)*param_fit_abs_4cm_wf2(6)*(j*omega).^param_fit_abs_4cm_wf2(7)) + param_fit_abs_4cm_wf2(8)./(1+param_fit_abs_4cm_wf2(8)*param_fit_abs_4cm_wf2(9)*(j*omega).^param_fit_abs_4cm_wf2(10));
Z_fit_abs_4cm_wf3 = param_fit_abs_4cm_wf3(1) + param_fit_abs_4cm_wf3(2)./(1+param_fit_abs_4cm_wf3(2)*param_fit_abs_4cm_wf3(3)*(j*omega).^param_fit_abs_4cm_wf3(4)) + param_fit_abs_4cm_wf3(5)./(1+param_fit_abs_4cm_wf3(5)*param_fit_abs_4cm_wf3(6)*(j*omega).^param_fit_abs_4cm_wf3(7)) + param_fit_abs_4cm_wf3(8)./(1+param_fit_abs_4cm_wf3(8)*param_fit_abs_4cm_wf3(9)*(j*omega).^param_fit_abs_4cm_wf3(10));
Z_fit_abs_4cm_wf4 = param_fit_abs_4cm_wf4(1) + param_fit_abs_4cm_wf4(2)./(1+param_fit_abs_4cm_wf4(2)*param_fit_abs_4cm_wf4(3)*(j*omega).^param_fit_abs_4cm_wf4(4)) + param_fit_abs_4cm_wf4(5)./(1+param_fit_abs_4cm_wf4(5)*param_fit_abs_4cm_wf4(6)*(j*omega).^param_fit_abs_4cm_wf4(7)) + param_fit_abs_4cm_wf4(8)./(1+param_fit_abs_4cm_wf4(8)*param_fit_abs_4cm_wf4(9)*(j*omega).^param_fit_abs_4cm_wf4(10));
Z_fit_abs_4cm_wf5 = param_fit_abs_4cm_wf5(1) + param_fit_abs_4cm_wf5(2)./(1+param_fit_abs_4cm_wf5(2)*param_fit_abs_4cm_wf5(3)*(j*omega).^param_fit_abs_4cm_wf5(4)) + param_fit_abs_4cm_wf5(5)./(1+param_fit_abs_4cm_wf5(5)*param_fit_abs_4cm_wf5(6)*(j*omega).^param_fit_abs_4cm_wf5(7)) + param_fit_abs_4cm_wf5(8)./(1+param_fit_abs_4cm_wf5(8)*param_fit_abs_4cm_wf5(9)*(j*omega).^param_fit_abs_4cm_wf5(10));

Z_fit_abs_2cm_wf1 = param_fit_abs_2cm_wf1(1) + param_fit_abs_2cm_wf1(2)./(1+param_fit_abs_2cm_wf1(2)*param_fit_abs_2cm_wf1(3)*(j*omega).^param_fit_abs_2cm_wf1(4)) + param_fit_abs_2cm_wf1(5)./(1+param_fit_abs_2cm_wf1(5)*param_fit_abs_2cm_wf1(6)*(j*omega).^param_fit_abs_2cm_wf1(7)) + param_fit_abs_2cm_wf1(8)./(1+param_fit_abs_2cm_wf1(8)*param_fit_abs_2cm_wf1(9)*(j*omega).^param_fit_abs_2cm_wf1(10));
Z_fit_abs_2cm_wf2 = param_fit_abs_2cm_wf2(1) + param_fit_abs_2cm_wf2(2)./(1+param_fit_abs_2cm_wf2(2)*param_fit_abs_2cm_wf2(3)*(j*omega).^param_fit_abs_2cm_wf2(4)) + param_fit_abs_2cm_wf2(5)./(1+param_fit_abs_2cm_wf2(5)*param_fit_abs_2cm_wf2(6)*(j*omega).^param_fit_abs_2cm_wf2(7)) + param_fit_abs_2cm_wf2(8)./(1+param_fit_abs_2cm_wf2(8)*param_fit_abs_2cm_wf2(9)*(j*omega).^param_fit_abs_2cm_wf2(10));
Z_fit_abs_2cm_wf3 = param_fit_abs_2cm_wf3(1) + param_fit_abs_2cm_wf3(2)./(1+param_fit_abs_2cm_wf3(2)*param_fit_abs_2cm_wf3(3)*(j*omega).^param_fit_abs_2cm_wf3(4)) + param_fit_abs_2cm_wf3(5)./(1+param_fit_abs_2cm_wf3(5)*param_fit_abs_2cm_wf3(6)*(j*omega).^param_fit_abs_2cm_wf3(7)) + param_fit_abs_2cm_wf3(8)./(1+param_fit_abs_2cm_wf3(8)*param_fit_abs_2cm_wf3(9)*(j*omega).^param_fit_abs_2cm_wf3(10));
Z_fit_abs_2cm_wf4 = param_fit_abs_2cm_wf4(1) + param_fit_abs_2cm_wf4(2)./(1+param_fit_abs_2cm_wf4(2)*param_fit_abs_2cm_wf4(3)*(j*omega).^param_fit_abs_2cm_wf4(4)) + param_fit_abs_2cm_wf4(5)./(1+param_fit_abs_2cm_wf4(5)*param_fit_abs_2cm_wf4(6)*(j*omega).^param_fit_abs_2cm_wf4(7)) + param_fit_abs_2cm_wf4(8)./(1+param_fit_abs_2cm_wf4(8)*param_fit_abs_2cm_wf4(9)*(j*omega).^param_fit_abs_2cm_wf4(10));
Z_fit_abs_2cm_wf5 = param_fit_abs_2cm_wf5(1) + param_fit_abs_2cm_wf5(2)./(1+param_fit_abs_2cm_wf5(2)*param_fit_abs_2cm_wf5(3)*(j*omega).^param_fit_abs_2cm_wf5(4)) + param_fit_abs_2cm_wf5(5)./(1+param_fit_abs_2cm_wf5(5)*param_fit_abs_2cm_wf5(6)*(j*omega).^param_fit_abs_2cm_wf5(7)) + param_fit_abs_2cm_wf5(8)./(1+param_fit_abs_2cm_wf5(8)*param_fit_abs_2cm_wf5(9)*(j*omega).^param_fit_abs_2cm_wf5(10));

figure('WindowState','maximized')
hAx=axes;
hAx.XScale='log';
% hAx.ColorOrder = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1];
hAx.ColorOrder = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
grid minor
title('Lower Leg Bioimpedance Spectrum - Magnitude - Changes in Fat Thickness')
ylabel('Magnitude [Ohm]')
xlabel('Frequency [Hz]')
set(gca,'FontSize',20); hold all
plot(freq,abs(Z_4cm_wf1),'o--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_4cm_wf2),'o--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_4cm_wf3),'o--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_4cm_wf4),'o--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_4cm_wf5),'o--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_2cm_wf1),'s--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_2cm_wf2),'s--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_2cm_wf3),'s--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_2cm_wf4),'s--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_2cm_wf5),'s--','LineWidth',2,'MarkerSize',16);
   
plot(freq,abs(Z_fit_abs_4cm_wf1),'.-','LineWidth',3,'MarkerSize',45);
plot(freq,abs(Z_fit_abs_4cm_wf2),'.-','LineWidth',3,'MarkerSize',45);
plot(freq,abs(Z_fit_abs_4cm_wf3),'.-','LineWidth',3,'MarkerSize',45);
plot(freq,abs(Z_fit_abs_4cm_wf4),'.-','LineWidth',3,'MarkerSize',45);
plot(freq,abs(Z_fit_abs_4cm_wf5),'.-','LineWidth',3,'MarkerSize',45);
plot(freq,abs(Z_fit_abs_2cm_wf1),'x-','LineWidth',3,'MarkerSize',16);
plot(freq,abs(Z_fit_abs_2cm_wf2),'x-','LineWidth',3,'MarkerSize',16);
plot(freq,abs(Z_fit_abs_2cm_wf3),'x-','LineWidth',3,'MarkerSize',16);
plot(freq,abs(Z_fit_abs_2cm_wf4),'x-','LineWidth',3,'MarkerSize',16);
plot(freq,abs(Z_fit_abs_2cm_wf5),'x-','LineWidth',3,'MarkerSize',16);
legend('FEA - 4 cm IED - t_f = 4mm',...
       'FEA - 4 cm IED - t_f = 6mm',...
       'FEA - 4 cm IED - t_f = 8mm',...
       'FEA - 4 cm IED - t_f = 10mm',...
       'FEA - 4 cm IED - t_f = 12mm',...
       'FEA - 2.5 cm IED - t_f = 4mm',...
       'FEA - 2.5 cm IED - t_f = 6mm',...
       'FEA - 2.5 cm IED - t_f = 8mm',...
       'FEA - 2.5 cm IED - t_f = 10mm',...
       'FEA - 2.5 cm IED - t_f = 12mm',...
       'Circuit Model - 4 cm IED - t_f = 4mm',...
       'Circuit Model - 4 cm IED - t_f = 6mm',...
       'Circuit Model - 4 cm IED - t_f = 8mm',...
       'Circuit Model - 4 cm IED - t_f = 10mm',...
       'Circuit Model - 4 cm IED - t_f = 12mm',...
       'Circuit Model - 2.5 cm IED - t_f = 4mm',...
       'Circuit Model - 2.5 cm IED - t_f = 6mm',...
       'Circuit Model - 2.5 cm IED - t_f = 8mm',...
       'Circuit Model - 2.5 cm IED - t_f = 10mm',...
       'Circuit Model - 2.5 cm IED - t_f = 12mm','Location','northeastoutside')
% legend('FEA - 4 cm IED - tf = 4mm',...
%        'FEA - 4 cm IED - tf = 6mm',...
%        'FEA - 4 cm IED - tf = 8mm',...
%        'FEA - 4 cm IED - tf = 10mm',...
%        'FEA - 4 cm IED - tf = 12mm',...
%        'Circuit Model - 4 cm IED - tf = 4mm',...
%        'Circuit Model - 4 cm IED - tf = 6mm',...
%        'Circuit Model - 4 cm IED - tf = 8mm',...
%        'Circuit Model - 4 cm IED - tf = 10mm',...
%        'Circuit Model - 4 cm IED - tf = 12mm','Location','northeast')
% print -depsc m_leg_tf

figure('WindowState','maximized')
hAx=axes;
hAx.XScale='log';
% hAx.ColorOrder = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1];
hAx.ColorOrder = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
grid minor
title('Lower Leg Bioimpedance Spectrum - Phase - Changes in Fat Thickness')
ylabel('Phase [Degrees]')
xlabel('Frequency [Hz]')
set(gca,'FontSize',20); hold all
plot(freq,-180*angle(Z_4cm_wf1)/pi,'o--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_4cm_wf2)/pi,'o--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_4cm_wf3)/pi,'o--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_4cm_wf4)/pi,'o--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_4cm_wf5)/pi,'o--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_2cm_wf1)/pi,'s--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_2cm_wf2)/pi,'s--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_2cm_wf3)/pi,'s--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_2cm_wf4)/pi,'s--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_2cm_wf5)/pi,'s--','LineWidth',2,'MarkerSize',16);

plot(freq,-180*angle(Z_fit_abs_4cm_wf1)/pi,'.-','LineWidth',3,'MarkerSize',45);
plot(freq,-180*angle(Z_fit_abs_4cm_wf2)/pi,'.-','LineWidth',3,'MarkerSize',45);
plot(freq,-180*angle(Z_fit_abs_4cm_wf3)/pi,'.-','LineWidth',3,'MarkerSize',45);
plot(freq,-180*angle(Z_fit_abs_4cm_wf4)/pi,'.-','LineWidth',3,'MarkerSize',45);
plot(freq,-180*angle(Z_fit_abs_4cm_wf5)/pi,'.-','LineWidth',3,'MarkerSize',45);
plot(freq,-180*angle(Z_fit_abs_2cm_wf1)/pi,'x-','LineWidth',3,'MarkerSize',16);
plot(freq,-180*angle(Z_fit_abs_2cm_wf2)/pi,'x-','LineWidth',3,'MarkerSize',16);
plot(freq,-180*angle(Z_fit_abs_2cm_wf3)/pi,'x-','LineWidth',3,'MarkerSize',16);
plot(freq,-180*angle(Z_fit_abs_2cm_wf4)/pi,'x-','LineWidth',3,'MarkerSize',16);
plot(freq,-180*angle(Z_fit_abs_2cm_wf5)/pi,'x-','LineWidth',3,'MarkerSize',16);
legend('FEA - 4 cm IED - t_f = 4mm',...
       'FEA - 4 cm IED - t_f = 6mm',...
       'FEA - 4 cm IED - t_f = 8mm',...
       'FEA - 4 cm IED - t_f = 10mm',...
       'FEA - 4 cm IED - t_f = 12mm',...
       'FEA - 2.5 cm IED - t_f = 4mm',...
       'FEA - 2.5 cm IED - t_f = 6mm',...
       'FEA - 2.5 cm IED - t_f = 8mm',...
       'FEA - 2.5 cm IED - t_f = 10mm',...
       'FEA - 2.5 cm IED - t_f = 12mm',...
       'Circuit Model - 4 cm IED - t_f = 4mm',...
       'Circuit Model - 4 cm IED - t_f = 6mm',...
       'Circuit Model - 4 cm IED - t_f = 8mm',...
       'Circuit Model - 4 cm IED - t_f = 10mm',...
       'Circuit Model - 4 cm IED - t_f = 12mm',...
       'Circuit Model - 2.5 cm IED - t_f = 4mm',...
       'Circuit Model - 2.5 cm IED - t_f = 6mm',...
       'Circuit Model - 2.5 cm IED - t_f = 8mm',...
       'Circuit Model - 2.5 cm IED - t_f = 10mm',...
       'Circuit Model - 2.5 cm IED - t_f = 12mm','Location','northeastoutside')
% legend('FEA - 4 cm IED - tf = 4mm',...
%        'FEA - 4 cm IED - tf = 6mm',...
%        'FEA - 4 cm IED - tf = 8mm',...
%        'FEA - 4 cm IED - tf = 10mm',...
%        'FEA - 4 cm IED - tf = 12mm',...
%        'Circuit Model - 4 cm IED - tf = 4mm',...
%        'Circuit Model - 4 cm IED - tf = 6mm',...
%        'Circuit Model - 4 cm IED - tf = 8mm',...
%        'Circuit Model - 4 cm IED - tf = 10mm',...
%        'Circuit Model - 4 cm IED - tf = 12mm','Location','northwest')
% print -depsc p_leg_tf

wf_r1 = [param_fit_abs_2cm_wf1(2) param_fit_abs_2cm_wf2(2) param_fit_abs_2cm_wf3(2) param_fit_abs_2cm_wf4(2) param_fit_abs_2cm_wf5(2)];
wf_r2 = [param_fit_abs_2cm_wf1(5) param_fit_abs_2cm_wf2(5) param_fit_abs_2cm_wf3(5) param_fit_abs_2cm_wf4(5) param_fit_abs_2cm_wf5(5)];
wf_r3 = [param_fit_abs_2cm_wf1(8) param_fit_abs_2cm_wf2(8) param_fit_abs_2cm_wf3(8) param_fit_abs_2cm_wf4(8) param_fit_abs_2cm_wf5(8)];
wf_z1 = [param_fit_abs_2cm_wf1(3) param_fit_abs_2cm_wf2(3) param_fit_abs_2cm_wf3(3) param_fit_abs_2cm_wf4(3) param_fit_abs_2cm_wf5(3)];
wf_z2 = [param_fit_abs_2cm_wf1(6) param_fit_abs_2cm_wf2(6) param_fit_abs_2cm_wf3(6) param_fit_abs_2cm_wf4(6) param_fit_abs_2cm_wf5(6)];
wf_z3 = [param_fit_abs_2cm_wf1(9) param_fit_abs_2cm_wf2(9) param_fit_abs_2cm_wf3(9) param_fit_abs_2cm_wf4(9) param_fit_abs_2cm_wf5(9)];

fun_wf_r1 = @(a,wf)a(2)*(wf-wf0).^2 + a(1)*(wf-wf0) + wf_r1(1);
fun_wf_r2 = @(a,wf)a(2)*(wf-wf0).^2 + a(1)*(wf-wf0) + wf_r2(1);
fun_wf_r3 = @(a,wf)a(2)*(wf-wf0).^2 + a(1)*(wf-wf0) + wf_r3(1);
fun_wf_z1 = @(a,wf)a(2)*(wf-wf0).^2 + a(1)*(wf-wf0) + wf_z1(1);
fun_wf_z2 = @(a,wf)a(2)*(wf-wf0).^2 + a(1)*(wf-wf0) + wf_z2(1);
fun_wf_z3 = @(a,wf)a(2)*(wf-wf0).^2 + a(1)*(wf-wf0) + wf_z3(1);

a_wf_r1 = lsqcurvefit(fun_wf_r1,[0 0],wf,wf_r1,[],[],options);
a_wf_r2 = lsqcurvefit(fun_wf_r2,[0 0],wf,wf_r2,[],[],options);
a_wf_r3 = lsqcurvefit(fun_wf_r3,[0 0],wf,wf_r3,[],[],options);
a_wf_z1 = lsqcurvefit(fun_wf_z1,[0 0],wf,wf_z1,[],[],options);
a_wf_z2 = lsqcurvefit(fun_wf_z2,[0 0],wf,wf_z2,[],[],options);
a_wf_z3 = lsqcurvefit(fun_wf_z3,[0 0],wf,wf_z3,[],[],options);

wf_r1_fit = a_wf_r1(2)*(wf-wf0).^2 + a_wf_r1(1)*(wf-wf0) + wf_r1(1);
wf_r2_fit = a_wf_r2(2)*(wf-wf0).^2 + a_wf_r2(1)*(wf-wf0) + wf_r2(1);
wf_r3_fit = a_wf_r3(2)*(wf-wf0).^2 + a_wf_r3(1)*(wf-wf0) + wf_r3(1);
wf_z1_fit = a_wf_z1(2)*(wf-wf0).^2 + a_wf_z1(1)*(wf-wf0) + wf_z1(1);
wf_z2_fit = a_wf_z2(2)*(wf-wf0).^2 + a_wf_z2(1)*(wf-wf0) + wf_z2(1);
wf_z3_fit = a_wf_z3(2)*(wf-wf0).^2 + a_wf_z3(1)*(wf-wf0) + wf_z3(1);

figure('WindowState','maximized')
sgtitle('Lower Leg Circuit Model Parameters - 4 cm IED','FontSize',20)
subplot(2,3,1)
grid minor
ylabel('R_{1} [Ohm]')
xlabel('t_f [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
plot(wf/1e-3,wf_r1,'o','LineWidth',3,'MarkerSize',12);
plot(wf/1e-3,wf_r1_fit,'LineWidth',3,'MarkerSize',12);
subplot(2,3,2)
grid minor
ylabel('R_{2} [Ohm]')
xlabel('t_f [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
p1 = plot(wf/1e-3,wf_r2,'o','LineWidth',3,'MarkerSize',12);
p2 = plot(wf/1e-3,wf_r2_fit,'LineWidth',3,'MarkerSize',12);
subplot(2,3,3)
grid minor
ylabel('R_{3} [Ohm]')
xlabel('t_f [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
plot(wf/1e-3,wf_r3,'o','LineWidth',3,'MarkerSize',12);
plot(wf/1e-3,wf_r3_fit,'LineWidth',3,'MarkerSize',12);
subplot(2,3,4)
grid minor
ylabel('C_{1} [F]')
xlabel('t_f [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
plot(wf/1e-3,wf_z1,'o','LineWidth',3,'MarkerSize',12);
plot(wf/1e-3,wf_z1_fit,'LineWidth',3,'MarkerSize',12);
subplot(2,3,5)
grid minor
ylabel('C_{2} [F]')
xlabel('t_f [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
plot(wf/1e-3,wf_z2,'o','LineWidth',3,'MarkerSize',12);
plot(wf/1e-3,wf_z2_fit,'LineWidth',3,'MarkerSize',12);
subplot(2,3,6)
grid minor
ylabel('C_{3} [F]')
xlabel('t_f [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
plot(wf/1e-3,wf_z3,'o','LineWidth',3,'MarkerSize',12);
plot(wf/1e-3,wf_z3_fit,'LineWidth',3,'MarkerSize',12);
hL = legend([p1,p2],{'Circuit Model Parameter','Scaling Model'});
newLocation = [0.19 0.91 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newLocation,'Units', newUnits);
% print -depsc param_leg_4cm_tf

wf_r1_2 = [param_fit_abs_2cm_wf1(2) param_fit_abs_2cm_wf2(2) param_fit_abs_2cm_wf3(2) param_fit_abs_2cm_wf4(2) param_fit_abs_2cm_wf5(2)];
wf_r2_2 = [param_fit_abs_2cm_wf1(5) param_fit_abs_2cm_wf2(5) param_fit_abs_2cm_wf3(5) param_fit_abs_2cm_wf4(5) param_fit_abs_2cm_wf5(5)];
wf_r3_2 = [param_fit_abs_2cm_wf1(8) param_fit_abs_2cm_wf2(8) param_fit_abs_2cm_wf3(8) param_fit_abs_2cm_wf4(8) param_fit_abs_2cm_wf5(8)];
wf_z1_2 = [param_fit_abs_2cm_wf1(3) param_fit_abs_2cm_wf2(3) param_fit_abs_2cm_wf3(3) param_fit_abs_2cm_wf4(3) param_fit_abs_2cm_wf5(3)];
wf_z2_2 = [param_fit_abs_2cm_wf1(6) param_fit_abs_2cm_wf2(6) param_fit_abs_2cm_wf3(6) param_fit_abs_2cm_wf4(6) param_fit_abs_2cm_wf5(6)];
wf_z3_2 = [param_fit_abs_2cm_wf1(9) param_fit_abs_2cm_wf2(9) param_fit_abs_2cm_wf3(9) param_fit_abs_2cm_wf4(9) param_fit_abs_2cm_wf5(9)];

fun_wf_r1_2 = @(a,wf)a(2)*(wf-wf0).^2 + a(1)*(wf-wf0) + wf_r1_2(1);
fun_wf_r2_2 = @(a,wf)a(2)*(wf-wf0).^2 + a(1)*(wf-wf0) + wf_r2_2(1);
fun_wf_r3_2 = @(a,wf)a(2)*(wf-wf0).^2 + a(1)*(wf-wf0) + wf_r3_2(1);
fun_wf_z1_2 = @(a,wf)a(2)*(wf-wf0).^2 + a(1)*(wf-wf0) + wf_z1_2(1);
fun_wf_z2_2 = @(a,wf)a(2)*(wf-wf0).^2 + a(1)*(wf-wf0) + wf_z2_2(1);
fun_wf_z3_2 = @(a,wf)a(2)*(wf-wf0).^2 + a(1)*(wf-wf0) + wf_z3_2(1);

a_wf_r1_2 = lsqcurvefit(fun_wf_r1_2,[0 0],wf,wf_r1_2,[],[],options);
a_wf_r2_2 = lsqcurvefit(fun_wf_r2_2,[0 0],wf,wf_r2_2,[],[],options);
a_wf_r3_2 = lsqcurvefit(fun_wf_r3_2,[0 0],wf,wf_r3_2,[],[],options);
a_wf_z1_2 = lsqcurvefit(fun_wf_z1_2,[0 0],wf,wf_z1_2,[],[],options);
a_wf_z2_2 = lsqcurvefit(fun_wf_z2_2,[0 0],wf,wf_z2_2,[],[],options);
a_wf_z3_2 = lsqcurvefit(fun_wf_z3_2,[0 0],wf,wf_z3_2,[],[],options);

wf_r1_fit_2 = a_wf_r1_2(2)*(wf-wf0).^2 + a_wf_r1_2(1)*(wf-wf0) + wf_r1_2(1);
wf_r2_fit_2 = a_wf_r2_2(2)*(wf-wf0).^2 + a_wf_r2_2(1)*(wf-wf0) + wf_r2_2(1);
wf_r3_fit_2 = a_wf_r3_2(2)*(wf-wf0).^2 + a_wf_r3_2(1)*(wf-wf0) + wf_r3_2(1);
wf_z1_fit_2 = a_wf_z1_2(2)*(wf-wf0).^2 + a_wf_z1_2(1)*(wf-wf0) + wf_z1_2(1);
wf_z2_fit_2 = a_wf_z2_2(2)*(wf-wf0).^2 + a_wf_z2_2(1)*(wf-wf0) + wf_z2_2(1);
wf_z3_fit_2 = a_wf_z3_2(2)*(wf-wf0).^2 + a_wf_z3_2(1)*(wf-wf0) + wf_z3_2(1);

figure('WindowState','maximized')
sgtitle('Lower Leg Circuit Model Parameters - 2.5 cm IED','FontSize',20)
subplot(2,3,1)
grid minor
ylabel('R_{1} [Ohm]')
xlabel('t_f [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
plot(wf/1e-3,wf_r1_2,'o','LineWidth',3,'MarkerSize',12);
plot(wf/1e-3,wf_r1_fit_2,'LineWidth',3,'MarkerSize',12);
subplot(2,3,2)
grid minor
ylabel('R_{2} [Ohm]')
xlabel('t_f [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
p1 = plot(wf/1e-3,wf_r2_2,'o','LineWidth',3,'MarkerSize',12);
p2 = plot(wf/1e-3,wf_r2_fit_2,'LineWidth',3,'MarkerSize',12);
subplot(2,3,3)
grid minor
ylabel('R_{3} [Ohm]')
xlabel('t_f [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
plot(wf/1e-3,wf_r3_2,'o','LineWidth',3,'MarkerSize',12);
plot(wf/1e-3,wf_r3_fit_2,'LineWidth',3,'MarkerSize',12);
subplot(2,3,4)
grid minor
ylabel('C_{1} [F]')
xlabel('t_f [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
plot(wf/1e-3,wf_z1_2,'o','LineWidth',3,'MarkerSize',12);
plot(wf/1e-3,wf_z1_fit_2,'LineWidth',3,'MarkerSize',12);
subplot(2,3,5)
grid minor
ylabel('C_{2} [F]')
xlabel('t_f [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
plot(wf/1e-3,wf_z2_2,'o','LineWidth',3,'MarkerSize',12);
plot(wf/1e-3,wf_z2_fit_2,'LineWidth',3,'MarkerSize',12);
subplot(2,3,6)
grid minor
ylabel('C_{3} [F]')
xlabel('t_f [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
plot(wf/1e-3,wf_z3_2,'o','LineWidth',3,'MarkerSize',12);
plot(wf/1e-3,wf_z3_fit_2,'LineWidth',3,'MarkerSize',12);
hL = legend([p1,p2],{'Circuit Model Parameter','Scaling Model'});
newLocation = [0.19 0.91 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newLocation,'Units', newUnits);
% print -depsc param_leg_2cm_tf

%%%%%%%%%%%%%%%%%%%%%      Muscle Variations    %%%%%%%%%%%%%%%%%%%%%%%%%%%
lb_4cm = [4 0 0 0.98 0 0 0.64 0 0 0.83];
ub_4cm = [5 200 1e-3 0.99 200 1e-3 0.65 200 1e-3 0.84];
lb_2cm = [4 0 0 0.98 0 1e-4 0.64 0 0 0.83];
ub_2cm = [5 200 1.5e-3 0.99 200 1.5e-4 0.65 70 1e-3 0.84];

Z_2cm_wm1 = z_comsol_leg_wm(:,2);
Z_2cm_wm2 = z_comsol_leg_wm(:,5);
Z_2cm_wm3 = z_comsol_leg_wm(:,8);
Z_2cm_wm4 = z_comsol_leg_wm(:,11);
Z_2cm_wm5 = z_comsol_leg_wm(:,14);

Z_4cm_wm1 = z_comsol_leg_wm(:,17);
Z_4cm_wm2 = z_comsol_leg_wm(:,20);
Z_4cm_wm3 = z_comsol_leg_wm(:,23);
Z_4cm_wm4 = z_comsol_leg_wm(:,26);
Z_4cm_wm5 = z_comsol_leg_wm(:,29);

param_fit_abs_4cm_wm1 = lsqcurvefit(fun_abs,param0,omega,abs(Z_4cm_wm1),lb_4cm,ub_4cm,options);
param_fit_abs_4cm_wm2 = lsqcurvefit(fun_abs,param0,omega,abs(Z_4cm_wm2),lb_4cm,ub_4cm,options);
param_fit_abs_4cm_wm3 = lsqcurvefit(fun_abs,param0,omega,abs(Z_4cm_wm3),lb_4cm,ub_4cm,options);
param_fit_abs_4cm_wm4 = lsqcurvefit(fun_abs,param0,omega,abs(Z_4cm_wm4),lb_4cm,ub_4cm,options);
param_fit_abs_4cm_wm5 = lsqcurvefit(fun_abs,param0,omega,abs(Z_4cm_wm5),lb_4cm,ub_4cm,options);

param_fit_abs_2cm_wm1 = lsqcurvefit(fun_abs,param0,omega,abs(Z_2cm_wm1),lb_2cm,ub_2cm,options);
param_fit_abs_2cm_wm2 = lsqcurvefit(fun_abs,param0,omega,abs(Z_2cm_wm2),lb_2cm,ub_2cm,options);
param_fit_abs_2cm_wm3 = lsqcurvefit(fun_abs,param0,omega,abs(Z_2cm_wm3),lb_2cm,ub_2cm,options);
param_fit_abs_2cm_wm4 = lsqcurvefit(fun_abs,param0,omega,abs(Z_2cm_wm4),lb_2cm,ub_2cm,options);
param_fit_abs_2cm_wm5 = lsqcurvefit(fun_abs,param0,omega,abs(Z_2cm_wm5),lb_2cm,ub_2cm,options);

Z_fit_abs_4cm_wm1 = param_fit_abs_4cm_wm1(1) + param_fit_abs_4cm_wm1(2)./(1+param_fit_abs_4cm_wm1(2)*param_fit_abs_4cm_wm1(3)*(j*omega).^param_fit_abs_4cm_wm1(4)) + param_fit_abs_4cm_wm1(5)./(1+param_fit_abs_4cm_wm1(5)*param_fit_abs_4cm_wm1(6)*(j*omega).^param_fit_abs_4cm_wm1(7)) + param_fit_abs_4cm_wm1(8)./(1+param_fit_abs_4cm_wm1(8)*param_fit_abs_4cm_wm1(9)*(j*omega).^param_fit_abs_4cm_wm1(10));
Z_fit_abs_4cm_wm2 = param_fit_abs_4cm_wm2(1) + param_fit_abs_4cm_wm2(2)./(1+param_fit_abs_4cm_wm2(2)*param_fit_abs_4cm_wm2(3)*(j*omega).^param_fit_abs_4cm_wm2(4)) + param_fit_abs_4cm_wm2(5)./(1+param_fit_abs_4cm_wm2(5)*param_fit_abs_4cm_wm2(6)*(j*omega).^param_fit_abs_4cm_wm2(7)) + param_fit_abs_4cm_wm2(8)./(1+param_fit_abs_4cm_wm2(8)*param_fit_abs_4cm_wm2(9)*(j*omega).^param_fit_abs_4cm_wm2(10));
Z_fit_abs_4cm_wm3 = param_fit_abs_4cm_wm3(1) + param_fit_abs_4cm_wm3(2)./(1+param_fit_abs_4cm_wm3(2)*param_fit_abs_4cm_wm3(3)*(j*omega).^param_fit_abs_4cm_wm3(4)) + param_fit_abs_4cm_wm3(5)./(1+param_fit_abs_4cm_wm3(5)*param_fit_abs_4cm_wm3(6)*(j*omega).^param_fit_abs_4cm_wm3(7)) + param_fit_abs_4cm_wm3(8)./(1+param_fit_abs_4cm_wm3(8)*param_fit_abs_4cm_wm3(9)*(j*omega).^param_fit_abs_4cm_wm3(10));
Z_fit_abs_4cm_wm4 = param_fit_abs_4cm_wm4(1) + param_fit_abs_4cm_wm4(2)./(1+param_fit_abs_4cm_wm4(2)*param_fit_abs_4cm_wm4(3)*(j*omega).^param_fit_abs_4cm_wm4(4)) + param_fit_abs_4cm_wm4(5)./(1+param_fit_abs_4cm_wm4(5)*param_fit_abs_4cm_wm4(6)*(j*omega).^param_fit_abs_4cm_wm4(7)) + param_fit_abs_4cm_wm4(8)./(1+param_fit_abs_4cm_wm4(8)*param_fit_abs_4cm_wm4(9)*(j*omega).^param_fit_abs_4cm_wm4(10));
Z_fit_abs_4cm_wm5 = param_fit_abs_4cm_wm5(1) + param_fit_abs_4cm_wm5(2)./(1+param_fit_abs_4cm_wm5(2)*param_fit_abs_4cm_wm5(3)*(j*omega).^param_fit_abs_4cm_wm5(4)) + param_fit_abs_4cm_wm5(5)./(1+param_fit_abs_4cm_wm5(5)*param_fit_abs_4cm_wm5(6)*(j*omega).^param_fit_abs_4cm_wm5(7)) + param_fit_abs_4cm_wm5(8)./(1+param_fit_abs_4cm_wm5(8)*param_fit_abs_4cm_wm5(9)*(j*omega).^param_fit_abs_4cm_wm5(10));

Z_fit_abs_2cm_wm1 = param_fit_abs_2cm_wm1(1) + param_fit_abs_2cm_wm1(2)./(1+param_fit_abs_2cm_wm1(2)*param_fit_abs_2cm_wm1(3)*(j*omega).^param_fit_abs_2cm_wm1(4)) + param_fit_abs_2cm_wm1(5)./(1+param_fit_abs_2cm_wm1(5)*param_fit_abs_2cm_wm1(6)*(j*omega).^param_fit_abs_2cm_wm1(7)) + param_fit_abs_2cm_wm1(8)./(1+param_fit_abs_2cm_wm1(8)*param_fit_abs_2cm_wm1(9)*(j*omega).^param_fit_abs_2cm_wm1(10));
Z_fit_abs_2cm_wm2 = param_fit_abs_2cm_wm2(1) + param_fit_abs_2cm_wm2(2)./(1+param_fit_abs_2cm_wm2(2)*param_fit_abs_2cm_wm2(3)*(j*omega).^param_fit_abs_2cm_wm2(4)) + param_fit_abs_2cm_wm2(5)./(1+param_fit_abs_2cm_wm2(5)*param_fit_abs_2cm_wm2(6)*(j*omega).^param_fit_abs_2cm_wm2(7)) + param_fit_abs_2cm_wm2(8)./(1+param_fit_abs_2cm_wm2(8)*param_fit_abs_2cm_wm2(9)*(j*omega).^param_fit_abs_2cm_wm2(10));
Z_fit_abs_2cm_wm3 = param_fit_abs_2cm_wm3(1) + param_fit_abs_2cm_wm3(2)./(1+param_fit_abs_2cm_wm3(2)*param_fit_abs_2cm_wm3(3)*(j*omega).^param_fit_abs_2cm_wm3(4)) + param_fit_abs_2cm_wm3(5)./(1+param_fit_abs_2cm_wm3(5)*param_fit_abs_2cm_wm3(6)*(j*omega).^param_fit_abs_2cm_wm3(7)) + param_fit_abs_2cm_wm3(8)./(1+param_fit_abs_2cm_wm3(8)*param_fit_abs_2cm_wm3(9)*(j*omega).^param_fit_abs_2cm_wm3(10));
Z_fit_abs_2cm_wm4 = param_fit_abs_2cm_wm4(1) + param_fit_abs_2cm_wm4(2)./(1+param_fit_abs_2cm_wm4(2)*param_fit_abs_2cm_wm4(3)*(j*omega).^param_fit_abs_2cm_wm4(4)) + param_fit_abs_2cm_wm4(5)./(1+param_fit_abs_2cm_wm4(5)*param_fit_abs_2cm_wm4(6)*(j*omega).^param_fit_abs_2cm_wm4(7)) + param_fit_abs_2cm_wm4(8)./(1+param_fit_abs_2cm_wm4(8)*param_fit_abs_2cm_wm4(9)*(j*omega).^param_fit_abs_2cm_wm4(10));
Z_fit_abs_2cm_wm5 = param_fit_abs_2cm_wm5(1) + param_fit_abs_2cm_wm5(2)./(1+param_fit_abs_2cm_wm5(2)*param_fit_abs_2cm_wm5(3)*(j*omega).^param_fit_abs_2cm_wm5(4)) + param_fit_abs_2cm_wm5(5)./(1+param_fit_abs_2cm_wm5(5)*param_fit_abs_2cm_wm5(6)*(j*omega).^param_fit_abs_2cm_wm5(7)) + param_fit_abs_2cm_wm5(8)./(1+param_fit_abs_2cm_wm5(8)*param_fit_abs_2cm_wm5(9)*(j*omega).^param_fit_abs_2cm_wm5(10));

figure('WindowState','maximized')
hAx=axes;
hAx.XScale='log';
% hAx.ColorOrder = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1];
hAx.ColorOrder = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
grid minor
title('Lower Leg Bioimpedance Spectrum - Magnitude - Changes in Muscle Thickness')
ylabel('Magnitude [Ohm]')
xlabel('Frequency [Hz]')
set(gca,'FontSize',20); hold all
plot(freq,abs(Z_4cm_wm1),'o--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_4cm_wm2),'o--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_4cm_wm3),'o--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_4cm_wm4),'o--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_4cm_wm5),'o--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_2cm_wm1),'s--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_2cm_wm2),'s--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_2cm_wm3),'s--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_2cm_wm4),'s--','LineWidth',2,'MarkerSize',16);
plot(freq,abs(Z_2cm_wm5),'s--','LineWidth',2,'MarkerSize',16);

plot(freq,abs(Z_fit_abs_4cm_wm1),'.-','LineWidth',3,'MarkerSize',45);
plot(freq,abs(Z_fit_abs_4cm_wm2),'.-','LineWidth',3,'MarkerSize',45);
plot(freq,abs(Z_fit_abs_4cm_wm3),'.-','LineWidth',3,'MarkerSize',45);
plot(freq,abs(Z_fit_abs_4cm_wm4),'.-','LineWidth',3,'MarkerSize',45);
plot(freq,abs(Z_fit_abs_4cm_wm5),'.-','LineWidth',3,'MarkerSize',45);
plot(freq,abs(Z_fit_abs_2cm_wm1),'x-','LineWidth',3,'MarkerSize',16);
plot(freq,abs(Z_fit_abs_2cm_wm2),'x-','LineWidth',3,'MarkerSize',16);
plot(freq,abs(Z_fit_abs_2cm_wm3),'x-','LineWidth',3,'MarkerSize',16);
plot(freq,abs(Z_fit_abs_2cm_wm4),'x-','LineWidth',3,'MarkerSize',16);
plot(freq,abs(Z_fit_abs_2cm_wm5),'x-','LineWidth',3,'MarkerSize',16);
legend('FEA - 4 cm IED - t_m = 4mm',...
       'FEA - 4 cm IED - t_m = 6mm',...
       'FEA - 4 cm IED - t_m = 8mm',...
       'FEA - 4 cm IED - t_m = 10mm',...
       'FEA - 4 cm IED - t_m = 12mm',...
       'FEA - 2.5 cm IED - t_m = 4mm',...
       'FEA - 2.5 cm IED - t_m = 6mm',...
       'FEA - 2.5 cm IED - t_m = 8mm',...
       'FEA - 2.5 cm IED - t_m = 10mm',...
       'FEA - 2.5 cm IED - t_m = 12mm',...
       'Circuit Model - 4 cm IED - t_m = 4mm',...
       'Circuit Model - 4 cm IED - t_m = 6mm',...
       'Circuit Model - 4 cm IED - t_m = 8mm',...
       'Circuit Model - 4 cm IED - t_m = 10mm',...
       'Circuit Model - 4 cm IED - t_m = 12mm',...
       'Circuit Model - 2.5 cm IED - t_m = 4mm',...
       'Circuit Model - 2.5 cm IED - t_m = 6mm',...
       'Circuit Model - 2.5 cm IED - t_m = 8mm',...
       'Circuit Model - 2.5 cm IED - t_m = 10mm',...
       'Circuit Model - 2.5 cm IED - t_m = 12mm','Location','northeastoutside')
% legend('FEA - 4 cm IED - tm = 50mm',...
%        'FEA - 4 cm IED - tm = 52mm',...
%        'FEA - 4 cm IED - tm = 54mm',...
%        'FEA - 4 cm IED - tm = 56mm',...
%        'FEA - 4 cm IED - tm = 58mm',...
%        'Circuit Model - 4 cm IED - tm = 50mm',...
%        'Circuit Model - 4 cm IED - tm = 52mm',...
%        'Circuit Model - 4 cm IED - tm = 54mm',...
%        'Circuit Model - 4 cm IED - tm = 56mm',...
%        'Circuit Model - 4 cm IED - tm = 58mm','Location','northeast')
% print -depsc m_leg_tm
   
figure('WindowState','maximized')
hAx=axes;
hAx.XScale='log';
% hAx.ColorOrder = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1];
hAx.ColorOrder = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
grid minor
title('Lower Leg Bioimpedance Spectrum - Phase - Changes in Muscle Thickness')
ylabel('Phase [Degrees]')
xlabel('Frequency [Hz]')
set(gca,'FontSize',20); hold all
plot(freq,-180*angle(Z_4cm_wm1)/pi,'o--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_4cm_wm2)/pi,'o--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_4cm_wm3)/pi,'o--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_4cm_wm4)/pi,'o--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_4cm_wm5)/pi,'o--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_2cm_wm1)/pi,'s--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_2cm_wm2)/pi,'s--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_2cm_wm3)/pi,'s--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_2cm_wm4)/pi,'s--','LineWidth',2,'MarkerSize',16);
plot(freq,-180*angle(Z_2cm_wm5)/pi,'s--','LineWidth',2,'MarkerSize',16);

plot(freq,-180*angle(Z_fit_abs_4cm_wm1)/pi,'.-','LineWidth',3,'MarkerSize',45);
plot(freq,-180*angle(Z_fit_abs_4cm_wm2)/pi,'.-','LineWidth',3,'MarkerSize',45);
plot(freq,-180*angle(Z_fit_abs_4cm_wm3)/pi,'.-','LineWidth',3,'MarkerSize',45);
plot(freq,-180*angle(Z_fit_abs_4cm_wm4)/pi,'.-','LineWidth',3,'MarkerSize',45);
plot(freq,-180*angle(Z_fit_abs_4cm_wm5)/pi,'.-','LineWidth',3,'MarkerSize',45);
plot(freq,-180*angle(Z_fit_abs_2cm_wm1)/pi,'x-','LineWidth',3,'MarkerSize',16);
plot(freq,-180*angle(Z_fit_abs_2cm_wm2)/pi,'x-','LineWidth',3,'MarkerSize',16);
plot(freq,-180*angle(Z_fit_abs_2cm_wm3)/pi,'x-','LineWidth',3,'MarkerSize',16);
plot(freq,-180*angle(Z_fit_abs_2cm_wm4)/pi,'x-','LineWidth',3,'MarkerSize',16);
plot(freq,-180*angle(Z_fit_abs_2cm_wm5)/pi,'x-','LineWidth',3,'MarkerSize',16);
legend('FEA - 4 cm IED - t_m = 4mm',...
       'FEA - 4 cm IED - t_m = 6mm',...
       'FEA - 4 cm IED - t_m = 8mm',...
       'FEA - 4 cm IED - t_m = 10mm',...
       'FEA - 4 cm IED - t_m = 12mm',...
       'FEA - 2.5 cm IED - t_m = 4mm',...
       'FEA - 2.5 cm IED - t_m = 6mm',...
       'FEA - 2.5 cm IED - t_m = 8mm',...
       'FEA - 2.5 cm IED - t_m = 10mm',...
       'FEA - 2.5 cm IED - t_m = 12mm',...
       'Circuit Model - 4 cm IED - t_m = 4mm',...
       'Circuit Model - 4 cm IED - t_m = 6mm',...
       'Circuit Model - 4 cm IED - t_m = 8mm',...
       'Circuit Model - 4 cm IED - t_m = 10mm',...
       'Circuit Model - 4 cm IED - t_m = 12mm',...
       'Circuit Model - 2.5 cm IED - t_m = 4mm',...
       'Circuit Model - 2.5 cm IED - t_m = 6mm',...
       'Circuit Model - 2.5 cm IED - t_m = 8mm',...
       'Circuit Model - 2.5 cm IED - t_m = 10mm',...
       'Circuit Model - 2.5 cm IED - t_m = 12mm','Location','northeastoutside')
% legend('FEA - 4 cm IED - tm = 50mm',...
%        'FEA - 4 cm IED - tm = 52mm',...
%        'FEA - 4 cm IED - tm = 54mm',...
%        'FEA - 4 cm IED - tm = 56mm',...
%        'FEA - 4 cm IED - tm = 58mm',...
%        'Circuit Model - 4 cm IED - tm = 50mm',...
%        'Circuit Model - 4 cm IED - tm = 52mm',...
%        'Circuit Model - 4 cm IED - tm = 54mm',...
%        'Circuit Model - 4 cm IED - tm = 56mm',...
%        'Circuit Model - 4 cm IED - tm = 58mm','Location','northwest')
% print -depsc p_leg_tm

wm_r1 = [param_fit_abs_2cm_wm1(2) param_fit_abs_2cm_wm2(2) param_fit_abs_2cm_wm3(2) param_fit_abs_2cm_wm4(2) param_fit_abs_2cm_wm5(2)];
wm_r2 = [param_fit_abs_2cm_wm1(5) param_fit_abs_2cm_wm2(5) param_fit_abs_2cm_wm3(5) param_fit_abs_2cm_wm4(5) param_fit_abs_2cm_wm5(5)];
wm_r3 = [param_fit_abs_2cm_wm1(8) param_fit_abs_2cm_wm2(8) param_fit_abs_2cm_wm3(8) param_fit_abs_2cm_wm4(8) param_fit_abs_2cm_wm5(8)];
wm_z1 = [param_fit_abs_2cm_wm1(3) param_fit_abs_2cm_wm2(3) param_fit_abs_2cm_wm3(3) param_fit_abs_2cm_wm4(3) param_fit_abs_2cm_wm5(3)];
wm_z2 = [param_fit_abs_2cm_wm1(6) param_fit_abs_2cm_wm2(6) param_fit_abs_2cm_wm3(6) param_fit_abs_2cm_wm4(6) param_fit_abs_2cm_wm5(6)];
wm_z3 = [param_fit_abs_2cm_wm1(9) param_fit_abs_2cm_wm2(9) param_fit_abs_2cm_wm3(9) param_fit_abs_2cm_wm4(9) param_fit_abs_2cm_wm5(9)];

fun_wm_r1 = @(a,wm)a(2)*(wm-wm0).^2 + a(1)*(wm-wm0) + wm_r1(1);
fun_wm_r2 = @(a,wm)a(2)*(wm-wm0).^2 + a(1)*(wm-wm0) + wm_r2(1);
fun_wm_r3 = @(a,wm)a(2)*(wm-wm0).^2 + a(1)*(wm-wm0) + wm_r3(1);
fun_wm_z1 = @(a,wm)a(2)*(wm-wm0).^2 + a(1)*(wm-wm0) + wm_z1(1);
fun_wm_z2 = @(a,wm)a(2)*(wm-wm0).^2 + a(1)*(wm-wm0) + wm_z2(1);
fun_wm_z3 = @(a,wm)a(2)*(wm-wm0).^2 + a(1)*(wm-wm0) + wm_z3(1);

a_wm_r1 = lsqcurvefit(fun_wm_r1,[0 0],wm,wm_r1,[],[],options);
a_wm_r2 = lsqcurvefit(fun_wm_r2,[0 0],wm,wm_r2,[],[],options);
a_wm_r3 = lsqcurvefit(fun_wm_r3,[0 0],wm,wm_r3,[],[],options);
a_wm_z1 = lsqcurvefit(fun_wm_z1,[0 0],wm,wm_z1,[],[],options);
a_wm_z2 = lsqcurvefit(fun_wm_z2,[0 0],wm,wm_z2,[],[],options);
a_wm_z3 = lsqcurvefit(fun_wm_z3,[0 0],wm,wm_z3,[],[],options);

wm_r1_fit = a_wm_r1(2)*(wm-wm0).^2 + a_wm_r1(1)*(wm-wm0) + wm_r1(1);
wm_r2_fit = a_wm_r2(2)*(wm-wm0).^2 + a_wm_r2(1)*(wm-wm0) + wm_r2(1);
wm_r3_fit = a_wm_r3(2)*(wm-wm0).^2 + a_wm_r3(1)*(wm-wm0) + wm_r3(1);
wm_z1_fit = a_wm_z1(2)*(wm-wm0).^2 + a_wm_z1(1)*(wm-wm0) + wm_z1(1);
wm_z2_fit = a_wm_z2(2)*(wm-wm0).^2 + a_wm_z2(1)*(wm-wm0) + wm_z2(1);
wm_z3_fit = a_wm_z3(2)*(wm-wm0).^2 + a_wm_z3(1)*(wm-wm0) + wm_z3(1);

figure('WindowState','maximized')
sgtitle('Lower Leg Model Parameters -4 cm IED','FontSize',20)
subplot(2,3,1)
grid minor
ylabel('R_{1} [Ohm]')
xlabel('t_m [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
plot(wm/1e-3,wm_r1,'o','LineWidth',3,'MarkerSize',12);
plot(wm/1e-3,wm_r1_fit,'LineWidth',3,'MarkerSize',12);
xlim([50 58])
subplot(2,3,2)
grid minor
ylabel('R_{2} [Ohm]')
xlabel('t_m [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
p1 = plot(wm/1e-3,wm_r2,'o','LineWidth',3,'MarkerSize',12);
p2 = plot(wm/1e-3,wm_r2_fit,'LineWidth',3,'MarkerSize',12);
xlim([50 58])
subplot(2,3,3)
grid minor
ylabel('R_{3} [Ohm]')
xlabel('t_m [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
plot(wm/1e-3,wm_r3,'o','LineWidth',3,'MarkerSize',12);
plot(wm/1e-3,wm_r3_fit,'LineWidth',3,'MarkerSize',12);
xlim([50 58])
subplot(2,3,4)
grid minor
ylabel('C_{1} [F]')
xlabel('t_m [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
plot(wm/1e-3,wm_z1,'o','LineWidth',3,'MarkerSize',12);
plot(wm/1e-3,wm_z1_fit,'LineWidth',3,'MarkerSize',12);
xlim([50 58])
subplot(2,3,5)
grid minor
ylabel('C_{2} [F]')
xlabel('t_m [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
ytickformat('%.1f')
plot(wm/1e-3,wm_z2,'o','LineWidth',3,'MarkerSize',12);
plot(wm/1e-3,wm_z2_fit,'LineWidth',3,'MarkerSize',12);
xlim([50 58])
subplot(2,3,6)
grid minor
ylabel('C_{3} [F]')
xlabel('t_m [mm]')
hold on; set(gca,'FontSize',28); pbaspect([1.5 1 1])
plot(wm/1e-3,wm_z3,'o','LineWidth',3,'MarkerSize',12);
plot(wm/1e-3,wm_z3_fit,'LineWidth',3,'MarkerSize',12);
xlim([50 58])
hL = legend([p1,p2],{'Circuit Model Parameter','Scaling Model'});
newLocation = [0.19 0.91 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newLocation,'Units', newUnits);
% print -depsc param_leg_4cm_tm

wm_r1_2 = [param_fit_abs_2cm_wm1(2) param_fit_abs_2cm_wm2(2) param_fit_abs_2cm_wm3(2) param_fit_abs_2cm_wm4(2) param_fit_abs_2cm_wm5(2)];
wm_r2_2 = [param_fit_abs_2cm_wm1(5) param_fit_abs_2cm_wm2(5) param_fit_abs_2cm_wm3(5) param_fit_abs_2cm_wm4(5) param_fit_abs_2cm_wm5(5)];
wm_r3_2 = [param_fit_abs_2cm_wm1(8) param_fit_abs_2cm_wm2(8) param_fit_abs_2cm_wm3(8) param_fit_abs_2cm_wm4(8) param_fit_abs_2cm_wm5(8)];
wm_z1_2 = [param_fit_abs_2cm_wm1(3) param_fit_abs_2cm_wm2(3) param_fit_abs_2cm_wm3(3) param_fit_abs_2cm_wm4(3) param_fit_abs_2cm_wm5(3)];
wm_z2_2 = [param_fit_abs_2cm_wm1(6) param_fit_abs_2cm_wm2(6) param_fit_abs_2cm_wm3(6) param_fit_abs_2cm_wm4(6) param_fit_abs_2cm_wm5(6)];
wm_z3_2 = [param_fit_abs_2cm_wm1(9) param_fit_abs_2cm_wm2(9) param_fit_abs_2cm_wm3(9) param_fit_abs_2cm_wm4(9) param_fit_abs_2cm_wm5(9)];

fun_wm_r1_2 = @(a,wm)a(2)*(wm-wm0).^2 + a(1)*(wm-wm0) + wm_r1_2(1);
fun_wm_r2_2 = @(a,wm)a(2)*(wm-wm0).^2 + a(1)*(wm-wm0) + wm_r2_2(1);
fun_wm_r3_2 = @(a,wm)a(2)*(wm-wm0).^2 + a(1)*(wm-wm0) + wm_r3_2(1);
fun_wm_z1_2 = @(a,wm)a(2)*(wm-wm0).^2 + a(1)*(wm-wm0) + wm_z1_2(1);
fun_wm_z2_2 = @(a,wm)a(2)*(wm-wm0).^2 + a(1)*(wm-wm0) + wm_z2_2(1);
fun_wm_z3_2 = @(a,wm)a(2)*(wm-wm0).^2 + a(1)*(wm-wm0) + wm_z3_2(1);

a_wm_r1_2 = lsqcurvefit(fun_wm_r1_2,[0 0],wm,wm_r1_2,[],[],options);
a_wm_r2_2 = lsqcurvefit(fun_wm_r2_2,[0 0],wm,wm_r2_2,[],[],options);
a_wm_r3_2 = lsqcurvefit(fun_wm_r3_2,[0 0],wm,wm_r3_2,[],[],options);
a_wm_z1_2 = lsqcurvefit(fun_wm_z1_2,[0 0],wm,wm_z1_2,[],[],options);
a_wm_z2_2 = lsqcurvefit(fun_wm_z2_2,[0 0],wm,wm_z2_2,[],[],options);
a_wm_z3_2 = lsqcurvefit(fun_wm_z3_2,[0 0],wm,wm_z3_2,[],[],options);

wm_r1_fit_2 = a_wm_r1_2(2)*(wm-wm0).^2 + a_wm_r1_2(1)*(wm-wm0) + wm_r1_2(1);
wm_r2_fit_2 = a_wm_r2_2(2)*(wm-wm0).^2 + a_wm_r2_2(1)*(wm-wm0) + wm_r2_2(1);
wm_r3_fit_2 = a_wm_r3_2(2)*(wm-wm0).^2 + a_wm_r3_2(1)*(wm-wm0) + wm_r3_2(1);
wm_z1_fit_2 = a_wm_z1_2(2)*(wm-wm0).^2 + a_wm_z1_2(1)*(wm-wm0) + wm_z1_2(1);
wm_z2_fit_2 = a_wm_z2_2(2)*(wm-wm0).^2 + a_wm_z2_2(1)*(wm-wm0) + wm_z2_2(1);
wm_z3_fit_2 = a_wm_z3_2(2)*(wm-wm0).^2 + a_wm_z3_2(1)*(wm-wm0) + wm_z3_2(1);

figure('WindowState','maximized')
sgtitle('Lower Leg Circuit Model Parameters - 2.5 cm IED','FontSize',20)
subplot(2,3,1)
grid minor
ylabel('R_{1} [Ohm]')
xlabel('t_m [mm]')
hold on; set(gca,'FontSize',28);
plot(wm/1e-3,wm_r1_2,'o','LineWidth',3,'MarkerSize',12);
plot(wm/1e-3,wm_r1_fit_2,'LineWidth',3,'MarkerSize',12);
subplot(2,3,2)
grid minor
ylabel('R_{2} [Ohm]')
xlabel('t_m [mm]')
hold on; set(gca,'FontSize',28);
p1 = plot(wm/1e-3,wm_r2_2,'o','LineWidth',3,'MarkerSize',12);
p2 = plot(wm/1e-3,wm_r2_fit_2,'LineWidth',3,'MarkerSize',12);
subplot(2,3,3)
grid minor
ylabel('R_{3} [Ohm]')
xlabel('t_m [mm]')
hold on; set(gca,'FontSize',28);
plot(wm/1e-3,wm_r3_2,'o','LineWidth',3,'MarkerSize',12);
plot(wm/1e-3,wm_r3_fit_2,'LineWidth',3,'MarkerSize',12);
subplot(2,3,4)
grid minor
ylabel('C_{1} [F]')
xlabel('t_m [mm]')
hold on; set(gca,'FontSize',28);
plot(wm/1e-3,wm_z1_2,'o','LineWidth',3,'MarkerSize',12);
plot(wm/1e-3,wm_z1_fit_2,'LineWidth',3,'MarkerSize',12);
subplot(2,3,5)
grid minor
ylabel('C_{2} [F]')
xlabel('t_m [mm]')
ytickformat('%.1f')
hold on; set(gca,'FontSize',28);
plot(wm/1e-3,wm_z2_2,'o','LineWidth',3,'MarkerSize',12);
plot(wm/1e-3,wm_z2_fit_2,'LineWidth',3,'MarkerSize',12);
subplot(2,3,6)
grid minor
ylabel('C_{3} [F]')
xlabel('t_m [mm]')
hold on; set(gca,'FontSize',28);
plot(wm/1e-3,wm_z3_2,'o','LineWidth',3,'MarkerSize',12);
plot(wm/1e-3,wm_z3_fit_2,'LineWidth',3,'MarkerSize',12);
hL = legend([p1,p2],{'Circuit Model Parameter','Scaling Model'});
newLocation = [0.19 0.91 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newLocation,'Units', newUnits);
% print -depsc param_leg_2cm_tm