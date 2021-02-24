close all
clear all
clc

for a = 1:2
    for b = 1:2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subject = a; % 
limb = b;    % 1 = Arm ; 2 = Leg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f=[976.563	1953.13	3906.25	7812.5	15625	31250	62500	125000	250000	500000	1000000];

if (subject == 1) && (limb == 1)
    
    left_2cm_1 = parseXML('EIM_left_bicep_2cm_11_2021-02-22_T17-29-48.xml');
    left_2cm_2 = parseXML('EIM_left_bicep_2cm_12_2021-02-22_T17-30-29.xml');
    left_2cm_3 = parseXML('EIM_left_bicep_2cm_13_2021-02-22_T17-30-59.xml');
    left_2cm_4 = parseXML('EIM_left_bicep_2cm_14_2021-02-22_T17-31-34.xml');
    left_2cm_5 = parseXML('EIM_left_bicep_2cm_15_2021-02-22_T17-32-02.xml');
    left_4cm_1 = parseXML('EIM_left_bicep_4cm_11_2021-02-22_T17-38-00.xml');
    left_4cm_2 = parseXML('EIM_left_bicep_4cm_12_2021-02-22_T17-38-32.xml');
    left_4cm_3 = parseXML('EIM_left_bicep_4cm_13_2021-02-22_T17-39-20.xml');
    left_4cm_4 = parseXML('EIM_left_bicep_4cm_14_2021-02-22_T17-39-49.xml');
    left_4cm_5 = parseXML('EIM_left_bicep_4cm_15_2021-02-22_T17-40-21.xml');

    right_2cm_1 = parseXML('EIM_right_bicep_2cm_11_2021-02-22_T17-51-04.xml');
    right_2cm_2 = parseXML('EIM_right_bicep_2cm_12_2021-02-22_T17-51-34.xml');
    right_2cm_3 = parseXML('EIM_right_bicep_2cm_13_2021-02-22_T17-53-05.xml');
    right_2cm_4 = parseXML('EIM_right_bicep_2cm_14_2021-02-22_T17-54-56.xml');
    right_2cm_5 = parseXML('EIM_right_bicep_2cm_15_2021-02-22_T17-55-42.xml');
    right_4cm_1 = parseXML('EIM_right_bicep_4cm_11_2021-02-22_T18-00-49.xml');
    right_4cm_2 = parseXML('EIM_right_bicep_4cm_12_2021-02-22_T18-01-36.xml');
    right_4cm_3 = parseXML('EIM_right_bicep_4cm_13_2021-02-22_T18-02-14.xml');
    right_4cm_4 = parseXML('EIM_right_bicep_4cm_14_2021-02-22_T18-04-35.xml');
    right_4cm_5 = parseXML('EIM_right_bicep_4cm_15_2021-02-22_T18-05-11.xml');

    lb_4cm = [4 0 0 0.98 0 0 0.64 0 0 0.83];
    ub_4cm = [5 200 1e-3 0.99 200 1e-3 0.65 200 1e-3 0.84];
    lb_2cm = [4 0 0 0.98 0 3e-5 0.64 0 0 0.83];
    ub_2cm = [5 200 1.5e-4 0.99 200 7e-5 0.65 70 1e-3 0.84];
    param0 = [10 31 1/7e3 0.99 16 1/17e3 0.71 35 1/1.5e6 0.8];
    options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxFunctionEvaluations',6000,'StepTolerance',1e-6,'MaxIterations',1000)

    z_comsol_arm_wf = csvread('z_comsol_arm_wf.csv',5);
    z_comsol_arm_wm = csvread('z_comsol_arm_wm.csv',5);
    
    freq = z_comsol_arm_wf(:,1);
    omega = 2*pi*freq;

    fun_abs = @(param,freq)abs(param(1) + param(2)./(1+param(2)*param(3)*(j*omega).^param(4)) + param(5)./(1+param(5)*param(6)*(j*omega).^param(7)) + param(8)./(1+param(8)*param(9)*(j*omega).^param(10)));

    Z_2cm_wf1 = z_comsol_arm_wf(:,2);
    Z_2cm_wf2 = z_comsol_arm_wf(:,5);
    Z_2cm_wf3 = z_comsol_arm_wf(:,8);
    Z_2cm_wf4 = z_comsol_arm_wf(:,11);
    Z_2cm_wf5 = z_comsol_arm_wf(:,14);

    Z_4cm_wf1 = z_comsol_arm_wf(:,17);
    Z_4cm_wf2 = z_comsol_arm_wf(:,20);
    Z_4cm_wf3 = z_comsol_arm_wf(:,23);
    Z_4cm_wf4 = z_comsol_arm_wf(:,26);
    Z_4cm_wf5 = z_comsol_arm_wf(:,29);

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

    lb_4cm = [4 0 0 0.98 0 0 0.64 0 0 0.83];
    ub_4cm = [5 200 1e-3 0.99 200 1e-3 0.65 200 1e-3 0.84];
    lb_2cm = [4 0 0 0.98 0 1e-4 0.64 0 0 0.83];
    ub_2cm = [5 200 1.5e-3 0.99 200 1.5e-4 0.65 70 1e-3 0.84];

    Z_2cm_wm1 = z_comsol_arm_wm(:,2);
    Z_2cm_wm2 = z_comsol_arm_wm(:,5);
    Z_2cm_wm3 = z_comsol_arm_wm(:,8);
    Z_2cm_wm4 = z_comsol_arm_wm(:,11);
    Z_2cm_wm5 = z_comsol_arm_wm(:,14);

    Z_4cm_wm1 = z_comsol_arm_wm(:,17);
    Z_4cm_wm2 = z_comsol_arm_wm(:,20);
    Z_4cm_wm3 = z_comsol_arm_wm(:,23);
    Z_4cm_wm4 = z_comsol_arm_wm(:,26);
    Z_4cm_wm5 = z_comsol_arm_wm(:,29);

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
    
    a_sim_2cm = abs(Z_2cm_wf1(11:21));
    p_sim_2cm = -(180/pi)*phase(Z_2cm_wf1(11:21));
    a_sim_4cm = abs(Z_4cm_wf1(11:21));
    p_sim_4cm = -(180/pi)*phase(Z_4cm_wf1(11:21));
    a_cadence_2cm = abs(Z_fit_abs_2cm_wf1(11:21));
    p_cadence_2cm = -(180/pi)*phase(Z_fit_abs_2cm_wf1(11:21));
    a_cadence_4cm = abs(Z_fit_abs_4cm_wf1(11:21));
    p_cadence_4cm = -(180/pi)*phase(Z_fit_abs_4cm_wf1(11:21));
    
elseif (subject == 1) && (limb == 2)
    
    left_2cm_1 = parseXML('EIM_left_leg_2cm_11_2021-02-22_T18-15-16.xml');
    left_2cm_2 = parseXML('EIM_left_leg_2cm_12_2021-02-22_T18-16-15.xml');
    left_2cm_3 = parseXML('EIM_left_leg_2cm_13_2021-02-22_T18-16-41.xml');
    left_2cm_4 = parseXML('EIM_left_leg_2cm_14_2021-02-22_T18-17-08.xml');
    left_2cm_5 = parseXML('EIM_left_leg_2cm_15_2021-02-22_T18-17-37.xml');
    left_4cm_1 = parseXML('EIM_left_leg_4cm_11_2021-02-22_T18-21-38.xml');
    left_4cm_2 = parseXML('EIM_left_leg_4cm_12_2021-02-22_T18-22-09.xml');
    left_4cm_3 = parseXML('EIM_left_leg_4cm_13_2021-02-22_T18-22-49.xml');
    left_4cm_4 = parseXML('EIM_left_leg_4cm_14_2021-02-22_T18-23-18.xml');
    left_4cm_5 = parseXML('EIM_left_leg_4cm_15_2021-02-22_T18-24-01.xml');

    right_2cm_1 = parseXML('EIM_right_leg_2cm_11_2021-02-22_T18-34-59.xml');
    right_2cm_2 = parseXML('EIM_right_leg_2cm_12_2021-02-22_T18-35-32.xml');
    right_2cm_3 = parseXML('EIM_right_leg_2cm_13_2021-02-22_T18-36-00.xml');
    right_2cm_4 = parseXML('EIM_right_leg_2cm_14_2021-02-22_T18-37-45.xml');
    right_2cm_5 = parseXML('EIM_right_leg_2cm_15_2021-02-22_T18-38-24.xml');
    right_4cm_1 = parseXML('EIM_right_leg_4cm_11_2021-02-22_T18-43-26.xml');
    right_4cm_2 = parseXML('EIM_right_leg_4cm_12_2021-02-22_T18-43-56.xml');
    right_4cm_3 = parseXML('EIM_right_leg_4cm_13_2021-02-22_T18-44-23.xml');
    right_4cm_4 = parseXML('EIM_right_leg_4cm_14_2021-02-22_T18-44-53.xml');
    right_4cm_5 = parseXML('EIM_right_leg_4cm_15_2021-02-22_T18-45-25.xml');

    lb_4cm = [4 0 0 0.98 0 0 0.64 0 0 0.83];
    ub_4cm = [5 200 1e-3 0.99 200 1e-3 0.65 200 1e-3 0.84];
    lb_2cm = [4 0 0 0.98 0 3e-5 0.64 0 0 0.83];
    ub_2cm = [5 200 1.5e-4 0.99 200 7e-5 0.65 70 1e-3 0.84];
    param0 = [10 31 1/7e3 0.99 16 1/17e3 0.71 35 1/1.5e6 0.8];
    options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxFunctionEvaluations',6000,'StepTolerance',1e-6,'MaxIterations',1000)

    z_comsol_leg_wf = csvread('z_comsol_leg_wf.csv',5);
    z_comsol_leg_wm = csvread('z_comsol_leg_wm.csv',5);

    freq = z_comsol_arm_wf(:,1);
    omega = 2*pi*freq;

    fun_abs = @(param,freq)abs(param(1) + param(2)./(1+param(2)*param(3)*(j*omega).^param(4)) + param(5)./(1+param(5)*param(6)*(j*omega).^param(7)) + param(8)./(1+param(8)*param(9)*(j*omega).^param(10)));

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

    a_sim_2cm = abs(Z_2cm_wf2(11:21));
    p_sim_2cm = -(180/pi)*phase(Z_2cm_wf2(11:21));
    a_sim_4cm = abs(Z_4cm_wf2(11:21));
    p_sim_4cm = -(180/pi)*phase(Z_4cm_wf2(11:21));
    a_cadence_2cm = abs(Z_fit_abs_2cm_wf2(11:21));
    p_cadence_2cm = -(180/pi)*phase(Z_fit_abs_2cm_wf2(11:21));
    a_cadence_4cm = abs(Z_fit_abs_4cm_wf2(11:21));
    p_cadence_4cm = -(180/pi)*phase(Z_fit_abs_4cm_wf2(11:21));
    
elseif (subject == 2) && (limb == 1)

    left_2cm_1 = parseXML('EIM_left_bicep_2cm_6_2021-02-22_T10-35-45.xml');
    left_2cm_2 = parseXML('EIM_left_bicep_2cm_7_2021-02-22_T10-36-16.xml');
    left_2cm_3 = parseXML('EIM_left_bicep_2cm_8_2021-02-22_T10-36-42.xml');
    left_2cm_4 = parseXML('EIM_left_bicep_2cm_9_2021-02-22_T10-37-16.xml');
    left_2cm_5 = parseXML('EIM_left_bicep_2cm_10_2021-02-22_T10-37-43.xml');
    left_4cm_1 = parseXML('EIM_left_bicep_4cm_6_2021-02-22_T10-41-16.xml');
    left_4cm_2 = parseXML('EIM_left_bicep_4cm_7_2021-02-22_T10-41-43.xml');
    left_4cm_3 = parseXML('EIM_left_bicep_4cm_8_2021-02-22_T10-42-10.xml');
    left_4cm_4 = parseXML('EIM_left_bicep_4cm_9_2021-02-22_T10-42-36.xml');
    left_4cm_5 = parseXML('EIM_left_bicep_4cm_10_2021-02-22_T10-43-02.xml');

    right_2cm_1 = parseXML('EIM_right_bicep_2cm_6_2021-02-22_T10-51-11.xml');
    right_2cm_2 = parseXML('EIM_right_bicep_2cm_7_2021-02-22_T10-51-37.xml');
    right_2cm_3 = parseXML('EIM_right_bicep_2cm_8_2021-02-22_T10-52-04.xml');
    right_2cm_4 = parseXML('EIM_right_bicep_2cm_9_2021-02-22_T10-52-32.xml');
    right_2cm_5 = parseXML('EIM_right_bicep_2cm_10_2021-02-22_T10-52-58.xml');
    right_4cm_1 = parseXML('EIM_right_bicep_4cm_6_2021-02-22_T10-56-59.xml');
    right_4cm_2 = parseXML('EIM_right_bicep_4cm_7_2021-02-22_T10-57-25.xml');
    right_4cm_3 = parseXML('EIM_right_bicep_4cm_8_2021-02-22_T10-57-53.xml');
    right_4cm_4 = parseXML('EIM_right_bicep_4cm_9_2021-02-22_T10-58-19.xml');
    right_4cm_5 = parseXML('EIM_right_bicep_4cm_10_2021-02-22_T10-58-47.xml');
    
    lb_4cm = [4 0 0 0.98 0 0 0.64 0 0 0.83];
    ub_4cm = [5 200 1e-3 0.99 200 1e-3 0.65 200 1e-3 0.84];
    lb_2cm = [4 0 0 0.98 0 3e-5 0.64 0 0 0.83];
    ub_2cm = [5 200 1.5e-4 0.99 200 7e-5 0.65 70 1e-3 0.84];
    param0 = [10 31 1/7e3 0.99 16 1/17e3 0.71 35 1/1.5e6 0.8];
    options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxFunctionEvaluations',6000,'StepTolerance',1e-6,'MaxIterations',1000)

    z_comsol_arm_wf = csvread('z_comsol_arm_wf.csv',5);
    z_comsol_arm_wm = csvread('z_comsol_arm_wm.csv',5);
    
    freq = z_comsol_arm_wf(:,1);
    omega = 2*pi*freq;

    fun_abs = @(param,freq)abs(param(1) + param(2)./(1+param(2)*param(3)*(j*omega).^param(4)) + param(5)./(1+param(5)*param(6)*(j*omega).^param(7)) + param(8)./(1+param(8)*param(9)*(j*omega).^param(10)));

    Z_2cm_wf1 = z_comsol_arm_wf(:,2);
    Z_2cm_wf2 = z_comsol_arm_wf(:,5);
    Z_2cm_wf3 = z_comsol_arm_wf(:,8);
    Z_2cm_wf4 = z_comsol_arm_wf(:,11);
    Z_2cm_wf5 = z_comsol_arm_wf(:,14);

    Z_4cm_wf1 = z_comsol_arm_wf(:,17);
    Z_4cm_wf2 = z_comsol_arm_wf(:,20);
    Z_4cm_wf3 = z_comsol_arm_wf(:,23);
    Z_4cm_wf4 = z_comsol_arm_wf(:,26);
    Z_4cm_wf5 = z_comsol_arm_wf(:,29);

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

    lb_4cm = [4 0 0 0.98 0 0 0.64 0 0 0.83];
    ub_4cm = [5 200 1e-3 0.99 200 1e-3 0.65 200 1e-3 0.84];
    lb_2cm = [4 0 0 0.98 0 1e-4 0.64 0 0 0.83];
    ub_2cm = [5 200 1.5e-3 0.99 200 1.5e-4 0.65 70 1e-3 0.84];

    Z_2cm_wm1 = z_comsol_arm_wm(:,2);
    Z_2cm_wm2 = z_comsol_arm_wm(:,5);
    Z_2cm_wm3 = z_comsol_arm_wm(:,8);
    Z_2cm_wm4 = z_comsol_arm_wm(:,11);
    Z_2cm_wm5 = z_comsol_arm_wm(:,14);

    Z_4cm_wm1 = z_comsol_arm_wm(:,17);
    Z_4cm_wm2 = z_comsol_arm_wm(:,20);
    Z_4cm_wm3 = z_comsol_arm_wm(:,23);
    Z_4cm_wm4 = z_comsol_arm_wm(:,26);
    Z_4cm_wm5 = z_comsol_arm_wm(:,29);

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
    
    a_sim_2cm = abs(Z_2cm_wf5(11:21));
    p_sim_2cm = -(180/pi)*phase(Z_2cm_wf5(11:21));
    a_sim_4cm = abs(Z_4cm_wf5(11:21));
    p_sim_4cm = -(180/pi)*phase(Z_4cm_wf5(11:21));
    a_cadence_2cm = abs(Z_fit_abs_2cm_wf5(11:21));
    p_cadence_2cm = -(180/pi)*phase(Z_fit_abs_2cm_wf5(11:21));
    a_cadence_4cm = abs(Z_fit_abs_4cm_wf5(11:21));
    p_cadence_4cm = -(180/pi)*phase(Z_fit_abs_4cm_wf5(11:21));
    
elseif (subject == 2) && (limb == 2)
    
    left_2cm_1 = parseXML('EIM_left_leg_2cm_6_2021-02-22_T11-11-42.xml');
    left_2cm_2 = parseXML('EIM_left_leg_2cm_7_2021-02-22_T11-12-59.xml');
    left_2cm_3 = parseXML('EIM_left_leg_2cm_8_2021-02-22_T11-13-47.xml');
    left_2cm_4 = parseXML('EIM_left_leg_2cm_9_2021-02-22_T11-14-37.xml');
    left_2cm_5 = parseXML('EIM_left_leg_2cm_10_2021-02-22_T11-15-42.xml');
    left_4cm_1 = parseXML('EIM_left_leg_4cm_6_2021-02-22_T11-25-01.xml');
    left_4cm_2 = parseXML('EIM_left_leg_4cm_7_2021-02-22_T11-25-28.xml');
    left_4cm_3 = parseXML('EIM_left_leg_4cm_8_2021-02-22_T11-25-55.xml');
    left_4cm_4 = parseXML('EIM_left_leg_4cm_9_2021-02-22_T11-26-24.xml');
    left_4cm_5 = parseXML('EIM_left_leg_4cm_10_2021-02-22_T11-26-53.xml');

    right_2cm_1 = parseXML('EIM_right_leg_2cm_6_2021-02-22_T11-30-47.xml');
    right_2cm_2 = parseXML('EIM_right_leg_2cm_7_2021-02-22_T11-31-16.xml');
    right_2cm_3 = parseXML('EIM_right_leg_2cm_8_2021-02-22_T11-31-43.xml');
    right_2cm_4 = parseXML('EIM_right_leg_2cm_9_2021-02-22_T11-32-09.xml');
    right_2cm_5 = parseXML('EIM_right_leg_2cm_10_2021-02-22_T11-32-35.xml');
    right_4cm_1 = parseXML('EIM_right_leg_4cm_6_2021-02-22_T11-35-58.xml');
    right_4cm_2 = parseXML('EIM_right_leg_4cm_7_2021-02-22_T11-36-25.xml');
    right_4cm_3 = parseXML('EIM_right_leg_4cm_8_2021-02-22_T11-36-51.xml');
    right_4cm_4 = parseXML('EIM_right_leg_4cm_9_2021-02-22_T11-37-18.xml');
    right_4cm_5 = parseXML('EIM_right_leg_4cm_10_2021-02-22_T11-37-46.xml');
    
    lb_4cm = [4 0 0 0.98 0 0 0.64 0 0 0.83];
    ub_4cm = [5 200 1e-3 0.99 200 1e-3 0.65 200 1e-3 0.84];
    lb_2cm = [4 0 0 0.98 0 3e-5 0.64 0 0 0.83];
    ub_2cm = [5 200 1.5e-4 0.99 200 7e-5 0.65 70 1e-3 0.84];
    param0 = [10 31 1/7e3 0.99 16 1/17e3 0.71 35 1/1.5e6 0.8];
    options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxFunctionEvaluations',6000,'StepTolerance',1e-6,'MaxIterations',1000)

    z_comsol_leg_wf = csvread('z_comsol_leg_wf.csv',5);
    z_comsol_leg_wm = csvread('z_comsol_leg_wm.csv',5);

    freq = z_comsol_arm_wf(:,1);
    omega = 2*pi*freq;

    fun_abs = @(param,freq)abs(param(1) + param(2)./(1+param(2)*param(3)*(j*omega).^param(4)) + param(5)./(1+param(5)*param(6)*(j*omega).^param(7)) + param(8)./(1+param(8)*param(9)*(j*omega).^param(10)));

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

    a_sim_2cm = abs(Z_2cm_wf2(11:21));
    p_sim_2cm = -(180/pi)*phase(Z_2cm_wf2(11:21));
    a_sim_4cm = abs(Z_4cm_wf2(11:21));
    p_sim_4cm = -(180/pi)*phase(Z_4cm_wf2(11:21));
    a_cadence_2cm = abs(Z_fit_abs_2cm_wf2(11:21));
    p_cadence_2cm = -(180/pi)*phase(Z_fit_abs_2cm_wf2(11:21));
    a_cadence_4cm = abs(Z_fit_abs_4cm_wf2(11:21));
    p_cadence_4cm = -(180/pi)*phase(Z_fit_abs_4cm_wf2(11:21));
    
end
 
a_left_2cm = [str2double(left_2cm_1.Children(2).Children(10).Children(2).Children.Data) str2double(left_2cm_1.Children(2).Children(10).Children(4).Children.Data) str2double(left_2cm_1.Children(2).Children(10).Children(6).Children.Data) str2double(left_2cm_1.Children(2).Children(10).Children(8).Children.Data) str2double(left_2cm_1.Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm_1.Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm_1.Children(2).Children(10).Children(14).Children.Data) str2double(left_2cm_1.Children(2).Children(10).Children(16).Children.Data) str2double(left_2cm_1.Children(2).Children(10).Children(18).Children.Data) str2double(left_2cm_1.Children(2).Children(10).Children(20).Children.Data) str2double(left_2cm_1.Children(2).Children(10).Children(22).Children.Data)
                       str2double(left_2cm_2.Children(2).Children(10).Children(2).Children.Data) str2double(left_2cm_2.Children(2).Children(10).Children(4).Children.Data) str2double(left_2cm_2.Children(2).Children(10).Children(6).Children.Data) str2double(left_2cm_2.Children(2).Children(10).Children(8).Children.Data) str2double(left_2cm_2.Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm_2.Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm_2.Children(2).Children(10).Children(14).Children.Data) str2double(left_2cm_2.Children(2).Children(10).Children(16).Children.Data) str2double(left_2cm_2.Children(2).Children(10).Children(18).Children.Data) str2double(left_2cm_2.Children(2).Children(10).Children(20).Children.Data) str2double(left_2cm_2.Children(2).Children(10).Children(22).Children.Data)
                       str2double(left_2cm_3.Children(2).Children(10).Children(2).Children.Data) str2double(left_2cm_3.Children(2).Children(10).Children(4).Children.Data) str2double(left_2cm_3.Children(2).Children(10).Children(6).Children.Data) str2double(left_2cm_3.Children(2).Children(10).Children(8).Children.Data) str2double(left_2cm_3.Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm_3.Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm_3.Children(2).Children(10).Children(14).Children.Data) str2double(left_2cm_3.Children(2).Children(10).Children(16).Children.Data) str2double(left_2cm_3.Children(2).Children(10).Children(18).Children.Data) str2double(left_2cm_3.Children(2).Children(10).Children(20).Children.Data) str2double(left_2cm_3.Children(2).Children(10).Children(22).Children.Data)
                       str2double(left_2cm_4.Children(2).Children(10).Children(2).Children.Data) str2double(left_2cm_4.Children(2).Children(10).Children(4).Children.Data) str2double(left_2cm_4.Children(2).Children(10).Children(6).Children.Data) str2double(left_2cm_4.Children(2).Children(10).Children(8).Children.Data) str2double(left_2cm_4.Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm_4.Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm_4.Children(2).Children(10).Children(14).Children.Data) str2double(left_2cm_4.Children(2).Children(10).Children(16).Children.Data) str2double(left_2cm_4.Children(2).Children(10).Children(18).Children.Data) str2double(left_2cm_4.Children(2).Children(10).Children(20).Children.Data) str2double(left_2cm_4.Children(2).Children(10).Children(22).Children.Data)
                       str2double(left_2cm_5.Children(2).Children(10).Children(2).Children.Data) str2double(left_2cm_5.Children(2).Children(10).Children(4).Children.Data) str2double(left_2cm_5.Children(2).Children(10).Children(6).Children.Data) str2double(left_2cm_5.Children(2).Children(10).Children(8).Children.Data) str2double(left_2cm_5.Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm_5.Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm_5.Children(2).Children(10).Children(14).Children.Data) str2double(left_2cm_5.Children(2).Children(10).Children(16).Children.Data) str2double(left_2cm_5.Children(2).Children(10).Children(18).Children.Data) str2double(left_2cm_5.Children(2).Children(10).Children(20).Children.Data) str2double(left_2cm_5.Children(2).Children(10).Children(22).Children.Data)];

a_right_2cm = [str2double(right_2cm_1.Children(2).Children(10).Children(2).Children.Data) str2double(right_2cm_1.Children(2).Children(10).Children(4).Children.Data) str2double(right_2cm_1.Children(2).Children(10).Children(6).Children.Data) str2double(right_2cm_1.Children(2).Children(10).Children(8).Children.Data) str2double(right_2cm_1.Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm_1.Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm_1.Children(2).Children(10).Children(14).Children.Data) str2double(right_2cm_1.Children(2).Children(10).Children(16).Children.Data) str2double(right_2cm_1.Children(2).Children(10).Children(18).Children.Data) str2double(right_2cm_1.Children(2).Children(10).Children(20).Children.Data) str2double(right_2cm_1.Children(2).Children(10).Children(22).Children.Data)
                       str2double(right_2cm_2.Children(2).Children(10).Children(2).Children.Data) str2double(right_2cm_2.Children(2).Children(10).Children(4).Children.Data) str2double(right_2cm_2.Children(2).Children(10).Children(6).Children.Data) str2double(right_2cm_2.Children(2).Children(10).Children(8).Children.Data) str2double(right_2cm_2.Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm_2.Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm_2.Children(2).Children(10).Children(14).Children.Data) str2double(right_2cm_2.Children(2).Children(10).Children(16).Children.Data) str2double(right_2cm_2.Children(2).Children(10).Children(18).Children.Data) str2double(right_2cm_2.Children(2).Children(10).Children(20).Children.Data) str2double(right_2cm_2.Children(2).Children(10).Children(22).Children.Data)
                       str2double(right_2cm_3.Children(2).Children(10).Children(2).Children.Data) str2double(right_2cm_3.Children(2).Children(10).Children(4).Children.Data) str2double(right_2cm_3.Children(2).Children(10).Children(6).Children.Data) str2double(right_2cm_3.Children(2).Children(10).Children(8).Children.Data) str2double(right_2cm_3.Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm_3.Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm_3.Children(2).Children(10).Children(14).Children.Data) str2double(right_2cm_3.Children(2).Children(10).Children(16).Children.Data) str2double(right_2cm_3.Children(2).Children(10).Children(18).Children.Data) str2double(right_2cm_3.Children(2).Children(10).Children(20).Children.Data) str2double(right_2cm_3.Children(2).Children(10).Children(22).Children.Data)
                       str2double(right_2cm_4.Children(2).Children(10).Children(2).Children.Data) str2double(right_2cm_4.Children(2).Children(10).Children(4).Children.Data) str2double(right_2cm_4.Children(2).Children(10).Children(6).Children.Data) str2double(right_2cm_4.Children(2).Children(10).Children(8).Children.Data) str2double(right_2cm_4.Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm_4.Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm_4.Children(2).Children(10).Children(14).Children.Data) str2double(right_2cm_4.Children(2).Children(10).Children(16).Children.Data) str2double(right_2cm_4.Children(2).Children(10).Children(18).Children.Data) str2double(right_2cm_4.Children(2).Children(10).Children(20).Children.Data) str2double(right_2cm_4.Children(2).Children(10).Children(22).Children.Data)
                       str2double(right_2cm_5.Children(2).Children(10).Children(2).Children.Data) str2double(right_2cm_5.Children(2).Children(10).Children(4).Children.Data) str2double(right_2cm_5.Children(2).Children(10).Children(6).Children.Data) str2double(right_2cm_5.Children(2).Children(10).Children(8).Children.Data) str2double(right_2cm_5.Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm_5.Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm_5.Children(2).Children(10).Children(14).Children.Data) str2double(right_2cm_5.Children(2).Children(10).Children(16).Children.Data) str2double(right_2cm_5.Children(2).Children(10).Children(18).Children.Data) str2double(right_2cm_5.Children(2).Children(10).Children(20).Children.Data) str2double(right_2cm_5.Children(2).Children(10).Children(22).Children.Data)];
                   
a_left_4cm = [str2double(left_4cm_1.Children(2).Children(10).Children(2).Children.Data) str2double(left_4cm_1.Children(2).Children(10).Children(4).Children.Data) str2double(left_4cm_1.Children(2).Children(10).Children(6).Children.Data) str2double(left_4cm_1.Children(2).Children(10).Children(8).Children.Data) str2double(left_4cm_1.Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm_1.Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm_1.Children(2).Children(10).Children(14).Children.Data) str2double(left_4cm_1.Children(2).Children(10).Children(16).Children.Data) str2double(left_4cm_1.Children(2).Children(10).Children(18).Children.Data) str2double(left_4cm_1.Children(2).Children(10).Children(20).Children.Data) str2double(left_4cm_1.Children(2).Children(10).Children(22).Children.Data)
                       str2double(left_4cm_2.Children(2).Children(10).Children(2).Children.Data) str2double(left_4cm_2.Children(2).Children(10).Children(4).Children.Data) str2double(left_4cm_2.Children(2).Children(10).Children(6).Children.Data) str2double(left_4cm_2.Children(2).Children(10).Children(8).Children.Data) str2double(left_4cm_2.Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm_2.Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm_2.Children(2).Children(10).Children(14).Children.Data) str2double(left_4cm_2.Children(2).Children(10).Children(16).Children.Data) str2double(left_4cm_2.Children(2).Children(10).Children(18).Children.Data) str2double(left_4cm_2.Children(2).Children(10).Children(20).Children.Data) str2double(left_4cm_2.Children(2).Children(10).Children(22).Children.Data)
                       str2double(left_4cm_3.Children(2).Children(10).Children(2).Children.Data) str2double(left_4cm_3.Children(2).Children(10).Children(4).Children.Data) str2double(left_4cm_3.Children(2).Children(10).Children(6).Children.Data) str2double(left_4cm_3.Children(2).Children(10).Children(8).Children.Data) str2double(left_4cm_3.Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm_3.Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm_3.Children(2).Children(10).Children(14).Children.Data) str2double(left_4cm_3.Children(2).Children(10).Children(16).Children.Data) str2double(left_4cm_3.Children(2).Children(10).Children(18).Children.Data) str2double(left_4cm_3.Children(2).Children(10).Children(20).Children.Data) str2double(left_4cm_3.Children(2).Children(10).Children(22).Children.Data)
                       str2double(left_4cm_4.Children(2).Children(10).Children(2).Children.Data) str2double(left_4cm_4.Children(2).Children(10).Children(4).Children.Data) str2double(left_4cm_4.Children(2).Children(10).Children(6).Children.Data) str2double(left_4cm_4.Children(2).Children(10).Children(8).Children.Data) str2double(left_4cm_4.Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm_4.Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm_4.Children(2).Children(10).Children(14).Children.Data) str2double(left_4cm_4.Children(2).Children(10).Children(16).Children.Data) str2double(left_4cm_4.Children(2).Children(10).Children(18).Children.Data) str2double(left_4cm_4.Children(2).Children(10).Children(20).Children.Data) str2double(left_4cm_4.Children(2).Children(10).Children(22).Children.Data)
                       str2double(left_4cm_5.Children(2).Children(10).Children(2).Children.Data) str2double(left_4cm_5.Children(2).Children(10).Children(4).Children.Data) str2double(left_4cm_5.Children(2).Children(10).Children(6).Children.Data) str2double(left_4cm_5.Children(2).Children(10).Children(8).Children.Data) str2double(left_4cm_5.Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm_5.Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm_5.Children(2).Children(10).Children(14).Children.Data) str2double(left_4cm_5.Children(2).Children(10).Children(16).Children.Data) str2double(left_4cm_5.Children(2).Children(10).Children(18).Children.Data) str2double(left_4cm_5.Children(2).Children(10).Children(20).Children.Data) str2double(left_4cm_5.Children(2).Children(10).Children(22).Children.Data)];
                   
a_right_4cm = [str2double(right_4cm_1.Children(2).Children(10).Children(2).Children.Data) str2double(right_4cm_1.Children(2).Children(10).Children(4).Children.Data) str2double(right_4cm_1.Children(2).Children(10).Children(6).Children.Data) str2double(right_4cm_1.Children(2).Children(10).Children(8).Children.Data) str2double(right_4cm_1.Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm_1.Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm_1.Children(2).Children(10).Children(14).Children.Data) str2double(right_4cm_1.Children(2).Children(10).Children(16).Children.Data) str2double(right_4cm_1.Children(2).Children(10).Children(18).Children.Data) str2double(right_4cm_1.Children(2).Children(10).Children(20).Children.Data) str2double(right_4cm_1.Children(2).Children(10).Children(22).Children.Data)
                       str2double(right_4cm_2.Children(2).Children(10).Children(2).Children.Data) str2double(right_4cm_2.Children(2).Children(10).Children(4).Children.Data) str2double(right_4cm_2.Children(2).Children(10).Children(6).Children.Data) str2double(right_4cm_2.Children(2).Children(10).Children(8).Children.Data) str2double(right_4cm_2.Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm_2.Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm_2.Children(2).Children(10).Children(14).Children.Data) str2double(right_4cm_2.Children(2).Children(10).Children(16).Children.Data) str2double(right_4cm_2.Children(2).Children(10).Children(18).Children.Data) str2double(right_4cm_2.Children(2).Children(10).Children(20).Children.Data) str2double(right_4cm_2.Children(2).Children(10).Children(22).Children.Data)
                       str2double(right_4cm_3.Children(2).Children(10).Children(2).Children.Data) str2double(right_4cm_3.Children(2).Children(10).Children(4).Children.Data) str2double(right_4cm_3.Children(2).Children(10).Children(6).Children.Data) str2double(right_4cm_3.Children(2).Children(10).Children(8).Children.Data) str2double(right_4cm_3.Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm_3.Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm_3.Children(2).Children(10).Children(14).Children.Data) str2double(right_4cm_3.Children(2).Children(10).Children(16).Children.Data) str2double(right_4cm_3.Children(2).Children(10).Children(18).Children.Data) str2double(right_4cm_3.Children(2).Children(10).Children(20).Children.Data) str2double(right_4cm_3.Children(2).Children(10).Children(22).Children.Data)
                       str2double(right_4cm_4.Children(2).Children(10).Children(2).Children.Data) str2double(right_4cm_4.Children(2).Children(10).Children(4).Children.Data) str2double(right_4cm_4.Children(2).Children(10).Children(6).Children.Data) str2double(right_4cm_4.Children(2).Children(10).Children(8).Children.Data) str2double(right_4cm_4.Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm_4.Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm_4.Children(2).Children(10).Children(14).Children.Data) str2double(right_4cm_4.Children(2).Children(10).Children(16).Children.Data) str2double(right_4cm_4.Children(2).Children(10).Children(18).Children.Data) str2double(right_4cm_4.Children(2).Children(10).Children(20).Children.Data) str2double(right_4cm_4.Children(2).Children(10).Children(22).Children.Data)
                       str2double(right_4cm_5.Children(2).Children(10).Children(2).Children.Data) str2double(right_4cm_5.Children(2).Children(10).Children(4).Children.Data) str2double(right_4cm_5.Children(2).Children(10).Children(6).Children.Data) str2double(right_4cm_5.Children(2).Children(10).Children(8).Children.Data) str2double(right_4cm_5.Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm_5.Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm_5.Children(2).Children(10).Children(14).Children.Data) str2double(right_4cm_5.Children(2).Children(10).Children(16).Children.Data) str2double(right_4cm_5.Children(2).Children(10).Children(18).Children.Data) str2double(right_4cm_5.Children(2).Children(10).Children(20).Children.Data) str2double(right_4cm_5.Children(2).Children(10).Children(22).Children.Data)];
                   
p_left_2cm = [str2double(left_2cm_1.Children(2).Children(12).Children(2).Children.Data) str2double(left_2cm_1.Children(2).Children(12).Children(4).Children.Data) str2double(left_2cm_1.Children(2).Children(12).Children(6).Children.Data) str2double(left_2cm_1.Children(2).Children(12).Children(8).Children.Data) str2double(left_2cm_1.Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm_1.Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm_1.Children(2).Children(12).Children(14).Children.Data) str2double(left_2cm_1.Children(2).Children(12).Children(16).Children.Data) str2double(left_2cm_1.Children(2).Children(12).Children(18).Children.Data) str2double(left_2cm_1.Children(2).Children(12).Children(20).Children.Data) str2double(left_2cm_1.Children(2).Children(12).Children(22).Children.Data)
                       str2double(left_2cm_2.Children(2).Children(12).Children(2).Children.Data) str2double(left_2cm_2.Children(2).Children(12).Children(4).Children.Data) str2double(left_2cm_2.Children(2).Children(12).Children(6).Children.Data) str2double(left_2cm_2.Children(2).Children(12).Children(8).Children.Data) str2double(left_2cm_2.Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm_2.Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm_2.Children(2).Children(12).Children(14).Children.Data) str2double(left_2cm_2.Children(2).Children(12).Children(16).Children.Data) str2double(left_2cm_2.Children(2).Children(12).Children(18).Children.Data) str2double(left_2cm_2.Children(2).Children(12).Children(20).Children.Data) str2double(left_2cm_2.Children(2).Children(12).Children(22).Children.Data)
                       str2double(left_2cm_3.Children(2).Children(12).Children(2).Children.Data) str2double(left_2cm_3.Children(2).Children(12).Children(4).Children.Data) str2double(left_2cm_3.Children(2).Children(12).Children(6).Children.Data) str2double(left_2cm_3.Children(2).Children(12).Children(8).Children.Data) str2double(left_2cm_3.Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm_3.Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm_3.Children(2).Children(12).Children(14).Children.Data) str2double(left_2cm_3.Children(2).Children(12).Children(16).Children.Data) str2double(left_2cm_3.Children(2).Children(12).Children(18).Children.Data) str2double(left_2cm_3.Children(2).Children(12).Children(20).Children.Data) str2double(left_2cm_3.Children(2).Children(12).Children(22).Children.Data)
                       str2double(left_2cm_4.Children(2).Children(12).Children(2).Children.Data) str2double(left_2cm_4.Children(2).Children(12).Children(4).Children.Data) str2double(left_2cm_4.Children(2).Children(12).Children(6).Children.Data) str2double(left_2cm_4.Children(2).Children(12).Children(8).Children.Data) str2double(left_2cm_4.Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm_4.Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm_4.Children(2).Children(12).Children(14).Children.Data) str2double(left_2cm_4.Children(2).Children(12).Children(16).Children.Data) str2double(left_2cm_4.Children(2).Children(12).Children(18).Children.Data) str2double(left_2cm_4.Children(2).Children(12).Children(20).Children.Data) str2double(left_2cm_4.Children(2).Children(12).Children(22).Children.Data)
                       str2double(left_2cm_5.Children(2).Children(12).Children(2).Children.Data) str2double(left_2cm_5.Children(2).Children(12).Children(4).Children.Data) str2double(left_2cm_5.Children(2).Children(12).Children(6).Children.Data) str2double(left_2cm_5.Children(2).Children(12).Children(8).Children.Data) str2double(left_2cm_5.Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm_5.Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm_5.Children(2).Children(12).Children(14).Children.Data) str2double(left_2cm_5.Children(2).Children(12).Children(16).Children.Data) str2double(left_2cm_5.Children(2).Children(12).Children(18).Children.Data) str2double(left_2cm_5.Children(2).Children(12).Children(20).Children.Data) str2double(left_2cm_5.Children(2).Children(12).Children(22).Children.Data)];

p_right_2cm = [str2double(right_2cm_1.Children(2).Children(12).Children(2).Children.Data) str2double(right_2cm_1.Children(2).Children(12).Children(4).Children.Data) str2double(right_2cm_1.Children(2).Children(12).Children(6).Children.Data) str2double(right_2cm_1.Children(2).Children(12).Children(8).Children.Data) str2double(right_2cm_1.Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm_1.Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm_1.Children(2).Children(12).Children(14).Children.Data) str2double(right_2cm_1.Children(2).Children(12).Children(16).Children.Data) str2double(right_2cm_1.Children(2).Children(12).Children(18).Children.Data) str2double(right_2cm_1.Children(2).Children(12).Children(20).Children.Data) str2double(right_2cm_1.Children(2).Children(12).Children(22).Children.Data)
                       str2double(right_2cm_2.Children(2).Children(12).Children(2).Children.Data) str2double(right_2cm_2.Children(2).Children(12).Children(4).Children.Data) str2double(right_2cm_2.Children(2).Children(12).Children(6).Children.Data) str2double(right_2cm_2.Children(2).Children(12).Children(8).Children.Data) str2double(right_2cm_2.Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm_2.Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm_2.Children(2).Children(12).Children(14).Children.Data) str2double(right_2cm_2.Children(2).Children(12).Children(16).Children.Data) str2double(right_2cm_2.Children(2).Children(12).Children(18).Children.Data) str2double(right_2cm_2.Children(2).Children(12).Children(20).Children.Data) str2double(right_2cm_2.Children(2).Children(12).Children(22).Children.Data)
                       str2double(right_2cm_3.Children(2).Children(12).Children(2).Children.Data) str2double(right_2cm_3.Children(2).Children(12).Children(4).Children.Data) str2double(right_2cm_3.Children(2).Children(12).Children(6).Children.Data) str2double(right_2cm_3.Children(2).Children(12).Children(8).Children.Data) str2double(right_2cm_3.Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm_3.Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm_3.Children(2).Children(12).Children(14).Children.Data) str2double(right_2cm_3.Children(2).Children(12).Children(16).Children.Data) str2double(right_2cm_3.Children(2).Children(12).Children(18).Children.Data) str2double(right_2cm_3.Children(2).Children(12).Children(20).Children.Data) str2double(right_2cm_3.Children(2).Children(12).Children(22).Children.Data)
                       str2double(right_2cm_4.Children(2).Children(12).Children(2).Children.Data) str2double(right_2cm_4.Children(2).Children(12).Children(4).Children.Data) str2double(right_2cm_4.Children(2).Children(12).Children(6).Children.Data) str2double(right_2cm_4.Children(2).Children(12).Children(8).Children.Data) str2double(right_2cm_4.Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm_4.Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm_4.Children(2).Children(12).Children(14).Children.Data) str2double(right_2cm_4.Children(2).Children(12).Children(16).Children.Data) str2double(right_2cm_4.Children(2).Children(12).Children(18).Children.Data) str2double(right_2cm_4.Children(2).Children(12).Children(20).Children.Data) str2double(right_2cm_4.Children(2).Children(12).Children(22).Children.Data)
                       str2double(right_2cm_5.Children(2).Children(12).Children(2).Children.Data) str2double(right_2cm_5.Children(2).Children(12).Children(4).Children.Data) str2double(right_2cm_5.Children(2).Children(12).Children(6).Children.Data) str2double(right_2cm_5.Children(2).Children(12).Children(8).Children.Data) str2double(right_2cm_5.Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm_5.Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm_5.Children(2).Children(12).Children(14).Children.Data) str2double(right_2cm_5.Children(2).Children(12).Children(16).Children.Data) str2double(right_2cm_5.Children(2).Children(12).Children(18).Children.Data) str2double(right_2cm_5.Children(2).Children(12).Children(20).Children.Data) str2double(right_2cm_5.Children(2).Children(12).Children(22).Children.Data)];
                   
p_left_4cm = [str2double(left_4cm_1.Children(2).Children(12).Children(2).Children.Data) str2double(left_4cm_1.Children(2).Children(12).Children(4).Children.Data) str2double(left_4cm_1.Children(2).Children(12).Children(6).Children.Data) str2double(left_4cm_1.Children(2).Children(12).Children(8).Children.Data) str2double(left_4cm_1.Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm_1.Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm_1.Children(2).Children(12).Children(14).Children.Data) str2double(left_4cm_1.Children(2).Children(12).Children(16).Children.Data) str2double(left_4cm_1.Children(2).Children(12).Children(18).Children.Data) str2double(left_4cm_1.Children(2).Children(12).Children(20).Children.Data) str2double(left_4cm_1.Children(2).Children(12).Children(22).Children.Data)
                       str2double(left_4cm_2.Children(2).Children(12).Children(2).Children.Data) str2double(left_4cm_2.Children(2).Children(12).Children(4).Children.Data) str2double(left_4cm_2.Children(2).Children(12).Children(6).Children.Data) str2double(left_4cm_2.Children(2).Children(12).Children(8).Children.Data) str2double(left_4cm_2.Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm_2.Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm_2.Children(2).Children(12).Children(14).Children.Data) str2double(left_4cm_2.Children(2).Children(12).Children(16).Children.Data) str2double(left_4cm_2.Children(2).Children(12).Children(18).Children.Data) str2double(left_4cm_2.Children(2).Children(12).Children(20).Children.Data) str2double(left_4cm_2.Children(2).Children(12).Children(22).Children.Data)
                       str2double(left_4cm_3.Children(2).Children(12).Children(2).Children.Data) str2double(left_4cm_3.Children(2).Children(12).Children(4).Children.Data) str2double(left_4cm_3.Children(2).Children(12).Children(6).Children.Data) str2double(left_4cm_3.Children(2).Children(12).Children(8).Children.Data) str2double(left_4cm_3.Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm_3.Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm_3.Children(2).Children(12).Children(14).Children.Data) str2double(left_4cm_3.Children(2).Children(12).Children(16).Children.Data) str2double(left_4cm_3.Children(2).Children(12).Children(18).Children.Data) str2double(left_4cm_3.Children(2).Children(12).Children(20).Children.Data) str2double(left_4cm_3.Children(2).Children(12).Children(22).Children.Data)
                       str2double(left_4cm_4.Children(2).Children(12).Children(2).Children.Data) str2double(left_4cm_4.Children(2).Children(12).Children(4).Children.Data) str2double(left_4cm_4.Children(2).Children(12).Children(6).Children.Data) str2double(left_4cm_4.Children(2).Children(12).Children(8).Children.Data) str2double(left_4cm_4.Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm_4.Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm_4.Children(2).Children(12).Children(14).Children.Data) str2double(left_4cm_4.Children(2).Children(12).Children(16).Children.Data) str2double(left_4cm_4.Children(2).Children(12).Children(18).Children.Data) str2double(left_4cm_4.Children(2).Children(12).Children(20).Children.Data) str2double(left_4cm_4.Children(2).Children(12).Children(22).Children.Data)
                       str2double(left_4cm_5.Children(2).Children(12).Children(2).Children.Data) str2double(left_4cm_5.Children(2).Children(12).Children(4).Children.Data) str2double(left_4cm_5.Children(2).Children(12).Children(6).Children.Data) str2double(left_4cm_5.Children(2).Children(12).Children(8).Children.Data) str2double(left_4cm_5.Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm_5.Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm_5.Children(2).Children(12).Children(14).Children.Data) str2double(left_4cm_5.Children(2).Children(12).Children(16).Children.Data) str2double(left_4cm_5.Children(2).Children(12).Children(18).Children.Data) str2double(left_4cm_5.Children(2).Children(12).Children(20).Children.Data) str2double(left_4cm_5.Children(2).Children(12).Children(22).Children.Data)];
                   
p_right_4cm = [str2double(right_4cm_1.Children(2).Children(12).Children(2).Children.Data) str2double(right_4cm_1.Children(2).Children(12).Children(4).Children.Data) str2double(right_4cm_1.Children(2).Children(12).Children(6).Children.Data) str2double(right_4cm_1.Children(2).Children(12).Children(8).Children.Data) str2double(right_4cm_1.Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm_1.Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm_1.Children(2).Children(12).Children(14).Children.Data) str2double(right_4cm_1.Children(2).Children(12).Children(16).Children.Data) str2double(right_4cm_1.Children(2).Children(12).Children(18).Children.Data) str2double(right_4cm_1.Children(2).Children(12).Children(20).Children.Data) str2double(right_4cm_1.Children(2).Children(12).Children(22).Children.Data)
                       str2double(right_4cm_2.Children(2).Children(12).Children(2).Children.Data) str2double(right_4cm_2.Children(2).Children(12).Children(4).Children.Data) str2double(right_4cm_2.Children(2).Children(12).Children(6).Children.Data) str2double(right_4cm_2.Children(2).Children(12).Children(8).Children.Data) str2double(right_4cm_2.Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm_2.Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm_2.Children(2).Children(12).Children(14).Children.Data) str2double(right_4cm_2.Children(2).Children(12).Children(16).Children.Data) str2double(right_4cm_2.Children(2).Children(12).Children(18).Children.Data) str2double(right_4cm_2.Children(2).Children(12).Children(20).Children.Data) str2double(right_4cm_2.Children(2).Children(12).Children(22).Children.Data)
                       str2double(right_4cm_3.Children(2).Children(12).Children(2).Children.Data) str2double(right_4cm_3.Children(2).Children(12).Children(4).Children.Data) str2double(right_4cm_3.Children(2).Children(12).Children(6).Children.Data) str2double(right_4cm_3.Children(2).Children(12).Children(8).Children.Data) str2double(right_4cm_3.Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm_3.Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm_3.Children(2).Children(12).Children(14).Children.Data) str2double(right_4cm_3.Children(2).Children(12).Children(16).Children.Data) str2double(right_4cm_3.Children(2).Children(12).Children(18).Children.Data) str2double(right_4cm_3.Children(2).Children(12).Children(20).Children.Data) str2double(right_4cm_3.Children(2).Children(12).Children(22).Children.Data)
                       str2double(right_4cm_4.Children(2).Children(12).Children(2).Children.Data) str2double(right_4cm_4.Children(2).Children(12).Children(4).Children.Data) str2double(right_4cm_4.Children(2).Children(12).Children(6).Children.Data) str2double(right_4cm_4.Children(2).Children(12).Children(8).Children.Data) str2double(right_4cm_4.Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm_4.Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm_4.Children(2).Children(12).Children(14).Children.Data) str2double(right_4cm_4.Children(2).Children(12).Children(16).Children.Data) str2double(right_4cm_4.Children(2).Children(12).Children(18).Children.Data) str2double(right_4cm_4.Children(2).Children(12).Children(20).Children.Data) str2double(right_4cm_4.Children(2).Children(12).Children(22).Children.Data)
                       str2double(right_4cm_5.Children(2).Children(12).Children(2).Children.Data) str2double(right_4cm_5.Children(2).Children(12).Children(4).Children.Data) str2double(right_4cm_5.Children(2).Children(12).Children(6).Children.Data) str2double(right_4cm_5.Children(2).Children(12).Children(8).Children.Data) str2double(right_4cm_5.Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm_5.Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm_5.Children(2).Children(12).Children(14).Children.Data) str2double(right_4cm_5.Children(2).Children(12).Children(16).Children.Data) str2double(right_4cm_5.Children(2).Children(12).Children(18).Children.Data) str2double(right_4cm_5.Children(2).Children(12).Children(20).Children.Data) str2double(right_4cm_5.Children(2).Children(12).Children(22).Children.Data)];
                   
ar_left_2cm = [str2double(left_2cm_1.Children(2).Children(14).Children(2).Children.Data) str2double(left_2cm_1.Children(2).Children(14).Children(4).Children.Data) str2double(left_2cm_1.Children(2).Children(14).Children(6).Children.Data) str2double(left_2cm_1.Children(2).Children(14).Children(8).Children.Data) str2double(left_2cm_1.Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm_1.Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm_1.Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm_1.Children(2).Children(14).Children(16).Children.Data) str2double(left_2cm_1.Children(2).Children(14).Children(18).Children.Data) str2double(left_2cm_1.Children(2).Children(14).Children(20).Children.Data) str2double(left_2cm_1.Children(2).Children(14).Children(22).Children.Data)
                       str2double(left_2cm_2.Children(2).Children(14).Children(2).Children.Data) str2double(left_2cm_2.Children(2).Children(14).Children(4).Children.Data) str2double(left_2cm_2.Children(2).Children(14).Children(6).Children.Data) str2double(left_2cm_2.Children(2).Children(14).Children(8).Children.Data) str2double(left_2cm_2.Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm_2.Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm_2.Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm_2.Children(2).Children(14).Children(16).Children.Data) str2double(left_2cm_2.Children(2).Children(14).Children(18).Children.Data) str2double(left_2cm_2.Children(2).Children(14).Children(20).Children.Data) str2double(left_2cm_2.Children(2).Children(14).Children(22).Children.Data)
                       str2double(left_2cm_3.Children(2).Children(14).Children(2).Children.Data) str2double(left_2cm_3.Children(2).Children(14).Children(4).Children.Data) str2double(left_2cm_3.Children(2).Children(14).Children(6).Children.Data) str2double(left_2cm_3.Children(2).Children(14).Children(8).Children.Data) str2double(left_2cm_3.Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm_3.Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm_3.Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm_3.Children(2).Children(14).Children(16).Children.Data) str2double(left_2cm_3.Children(2).Children(14).Children(18).Children.Data) str2double(left_2cm_3.Children(2).Children(14).Children(20).Children.Data) str2double(left_2cm_3.Children(2).Children(14).Children(22).Children.Data)
                       str2double(left_2cm_4.Children(2).Children(14).Children(2).Children.Data) str2double(left_2cm_4.Children(2).Children(14).Children(4).Children.Data) str2double(left_2cm_4.Children(2).Children(14).Children(6).Children.Data) str2double(left_2cm_4.Children(2).Children(14).Children(8).Children.Data) str2double(left_2cm_4.Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm_4.Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm_4.Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm_4.Children(2).Children(14).Children(16).Children.Data) str2double(left_2cm_4.Children(2).Children(14).Children(18).Children.Data) str2double(left_2cm_4.Children(2).Children(14).Children(20).Children.Data) str2double(left_2cm_4.Children(2).Children(14).Children(22).Children.Data)
                       str2double(left_2cm_5.Children(2).Children(14).Children(2).Children.Data) str2double(left_2cm_5.Children(2).Children(14).Children(4).Children.Data) str2double(left_2cm_5.Children(2).Children(14).Children(6).Children.Data) str2double(left_2cm_5.Children(2).Children(14).Children(8).Children.Data) str2double(left_2cm_5.Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm_5.Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm_5.Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm_5.Children(2).Children(14).Children(16).Children.Data) str2double(left_2cm_5.Children(2).Children(14).Children(18).Children.Data) str2double(left_2cm_5.Children(2).Children(14).Children(20).Children.Data) str2double(left_2cm_5.Children(2).Children(14).Children(22).Children.Data)];

ar_right_2cm = [str2double(right_2cm_1.Children(2).Children(14).Children(2).Children.Data) str2double(right_2cm_1.Children(2).Children(14).Children(4).Children.Data) str2double(right_2cm_1.Children(2).Children(14).Children(6).Children.Data) str2double(right_2cm_1.Children(2).Children(14).Children(8).Children.Data) str2double(right_2cm_1.Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm_1.Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm_1.Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm_1.Children(2).Children(14).Children(16).Children.Data) str2double(right_2cm_1.Children(2).Children(14).Children(18).Children.Data) str2double(right_2cm_1.Children(2).Children(14).Children(20).Children.Data) str2double(right_2cm_1.Children(2).Children(14).Children(22).Children.Data)
                       str2double(right_2cm_2.Children(2).Children(14).Children(2).Children.Data) str2double(right_2cm_2.Children(2).Children(14).Children(4).Children.Data) str2double(right_2cm_2.Children(2).Children(14).Children(6).Children.Data) str2double(right_2cm_2.Children(2).Children(14).Children(8).Children.Data) str2double(right_2cm_2.Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm_2.Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm_2.Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm_2.Children(2).Children(14).Children(16).Children.Data) str2double(right_2cm_2.Children(2).Children(14).Children(18).Children.Data) str2double(right_2cm_2.Children(2).Children(14).Children(20).Children.Data) str2double(right_2cm_2.Children(2).Children(14).Children(22).Children.Data)
                       str2double(right_2cm_3.Children(2).Children(14).Children(2).Children.Data) str2double(right_2cm_3.Children(2).Children(14).Children(4).Children.Data) str2double(right_2cm_3.Children(2).Children(14).Children(6).Children.Data) str2double(right_2cm_3.Children(2).Children(14).Children(8).Children.Data) str2double(right_2cm_3.Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm_3.Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm_3.Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm_3.Children(2).Children(14).Children(16).Children.Data) str2double(right_2cm_3.Children(2).Children(14).Children(18).Children.Data) str2double(right_2cm_3.Children(2).Children(14).Children(20).Children.Data) str2double(right_2cm_3.Children(2).Children(14).Children(22).Children.Data)
                       str2double(right_2cm_4.Children(2).Children(14).Children(2).Children.Data) str2double(right_2cm_4.Children(2).Children(14).Children(4).Children.Data) str2double(right_2cm_4.Children(2).Children(14).Children(6).Children.Data) str2double(right_2cm_4.Children(2).Children(14).Children(8).Children.Data) str2double(right_2cm_4.Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm_4.Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm_4.Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm_4.Children(2).Children(14).Children(16).Children.Data) str2double(right_2cm_4.Children(2).Children(14).Children(18).Children.Data) str2double(right_2cm_4.Children(2).Children(14).Children(20).Children.Data) str2double(right_2cm_4.Children(2).Children(14).Children(22).Children.Data)
                       str2double(right_2cm_5.Children(2).Children(14).Children(2).Children.Data) str2double(right_2cm_5.Children(2).Children(14).Children(4).Children.Data) str2double(right_2cm_5.Children(2).Children(14).Children(6).Children.Data) str2double(right_2cm_5.Children(2).Children(14).Children(8).Children.Data) str2double(right_2cm_5.Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm_5.Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm_5.Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm_5.Children(2).Children(14).Children(16).Children.Data) str2double(right_2cm_5.Children(2).Children(14).Children(18).Children.Data) str2double(right_2cm_5.Children(2).Children(14).Children(20).Children.Data) str2double(right_2cm_5.Children(2).Children(14).Children(22).Children.Data)];
                   
ar_left_4cm = [str2double(left_4cm_1.Children(2).Children(14).Children(2).Children.Data) str2double(left_4cm_1.Children(2).Children(14).Children(4).Children.Data) str2double(left_4cm_1.Children(2).Children(14).Children(6).Children.Data) str2double(left_4cm_1.Children(2).Children(14).Children(8).Children.Data) str2double(left_4cm_1.Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm_1.Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm_1.Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm_1.Children(2).Children(14).Children(16).Children.Data) str2double(left_4cm_1.Children(2).Children(14).Children(18).Children.Data) str2double(left_4cm_1.Children(2).Children(14).Children(20).Children.Data) str2double(left_4cm_1.Children(2).Children(14).Children(22).Children.Data)
                       str2double(left_4cm_2.Children(2).Children(14).Children(2).Children.Data) str2double(left_4cm_2.Children(2).Children(14).Children(4).Children.Data) str2double(left_4cm_2.Children(2).Children(14).Children(6).Children.Data) str2double(left_4cm_2.Children(2).Children(14).Children(8).Children.Data) str2double(left_4cm_2.Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm_2.Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm_2.Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm_2.Children(2).Children(14).Children(16).Children.Data) str2double(left_4cm_2.Children(2).Children(14).Children(18).Children.Data) str2double(left_4cm_2.Children(2).Children(14).Children(20).Children.Data) str2double(left_4cm_2.Children(2).Children(14).Children(22).Children.Data)
                       str2double(left_4cm_3.Children(2).Children(14).Children(2).Children.Data) str2double(left_4cm_3.Children(2).Children(14).Children(4).Children.Data) str2double(left_4cm_3.Children(2).Children(14).Children(6).Children.Data) str2double(left_4cm_3.Children(2).Children(14).Children(8).Children.Data) str2double(left_4cm_3.Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm_3.Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm_3.Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm_3.Children(2).Children(14).Children(16).Children.Data) str2double(left_4cm_3.Children(2).Children(14).Children(18).Children.Data) str2double(left_4cm_3.Children(2).Children(14).Children(20).Children.Data) str2double(left_4cm_3.Children(2).Children(14).Children(22).Children.Data)
                       str2double(left_4cm_4.Children(2).Children(14).Children(2).Children.Data) str2double(left_4cm_4.Children(2).Children(14).Children(4).Children.Data) str2double(left_4cm_4.Children(2).Children(14).Children(6).Children.Data) str2double(left_4cm_4.Children(2).Children(14).Children(8).Children.Data) str2double(left_4cm_4.Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm_4.Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm_4.Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm_4.Children(2).Children(14).Children(16).Children.Data) str2double(left_4cm_4.Children(2).Children(14).Children(18).Children.Data) str2double(left_4cm_4.Children(2).Children(14).Children(20).Children.Data) str2double(left_4cm_4.Children(2).Children(14).Children(22).Children.Data)
                       str2double(left_4cm_5.Children(2).Children(14).Children(2).Children.Data) str2double(left_4cm_5.Children(2).Children(14).Children(4).Children.Data) str2double(left_4cm_5.Children(2).Children(14).Children(6).Children.Data) str2double(left_4cm_5.Children(2).Children(14).Children(8).Children.Data) str2double(left_4cm_5.Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm_5.Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm_5.Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm_5.Children(2).Children(14).Children(16).Children.Data) str2double(left_4cm_5.Children(2).Children(14).Children(18).Children.Data) str2double(left_4cm_5.Children(2).Children(14).Children(20).Children.Data) str2double(left_4cm_5.Children(2).Children(14).Children(22).Children.Data)];
                   
ar_right_4cm = [str2double(right_4cm_1.Children(2).Children(14).Children(2).Children.Data) str2double(right_4cm_1.Children(2).Children(14).Children(4).Children.Data) str2double(right_4cm_1.Children(2).Children(14).Children(6).Children.Data) str2double(right_4cm_1.Children(2).Children(14).Children(8).Children.Data) str2double(right_4cm_1.Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm_1.Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm_1.Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm_1.Children(2).Children(14).Children(16).Children.Data) str2double(right_4cm_1.Children(2).Children(14).Children(18).Children.Data) str2double(right_4cm_1.Children(2).Children(14).Children(20).Children.Data) str2double(right_4cm_1.Children(2).Children(14).Children(22).Children.Data)
                       str2double(right_4cm_2.Children(2).Children(14).Children(2).Children.Data) str2double(right_4cm_2.Children(2).Children(14).Children(4).Children.Data) str2double(right_4cm_2.Children(2).Children(14).Children(6).Children.Data) str2double(right_4cm_2.Children(2).Children(14).Children(8).Children.Data) str2double(right_4cm_2.Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm_2.Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm_2.Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm_2.Children(2).Children(14).Children(16).Children.Data) str2double(right_4cm_2.Children(2).Children(14).Children(18).Children.Data) str2double(right_4cm_2.Children(2).Children(14).Children(20).Children.Data) str2double(right_4cm_2.Children(2).Children(14).Children(22).Children.Data)
                       str2double(right_4cm_3.Children(2).Children(14).Children(2).Children.Data) str2double(right_4cm_3.Children(2).Children(14).Children(4).Children.Data) str2double(right_4cm_3.Children(2).Children(14).Children(6).Children.Data) str2double(right_4cm_3.Children(2).Children(14).Children(8).Children.Data) str2double(right_4cm_3.Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm_3.Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm_3.Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm_3.Children(2).Children(14).Children(16).Children.Data) str2double(right_4cm_3.Children(2).Children(14).Children(18).Children.Data) str2double(right_4cm_3.Children(2).Children(14).Children(20).Children.Data) str2double(right_4cm_3.Children(2).Children(14).Children(22).Children.Data)
                       str2double(right_4cm_4.Children(2).Children(14).Children(2).Children.Data) str2double(right_4cm_4.Children(2).Children(14).Children(4).Children.Data) str2double(right_4cm_4.Children(2).Children(14).Children(6).Children.Data) str2double(right_4cm_4.Children(2).Children(14).Children(8).Children.Data) str2double(right_4cm_4.Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm_4.Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm_4.Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm_4.Children(2).Children(14).Children(16).Children.Data) str2double(right_4cm_4.Children(2).Children(14).Children(18).Children.Data) str2double(right_4cm_4.Children(2).Children(14).Children(20).Children.Data) str2double(right_4cm_4.Children(2).Children(14).Children(22).Children.Data)
                       str2double(right_4cm_5.Children(2).Children(14).Children(2).Children.Data) str2double(right_4cm_5.Children(2).Children(14).Children(4).Children.Data) str2double(right_4cm_5.Children(2).Children(14).Children(6).Children.Data) str2double(right_4cm_5.Children(2).Children(14).Children(8).Children.Data) str2double(right_4cm_5.Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm_5.Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm_5.Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm_5.Children(2).Children(14).Children(16).Children.Data) str2double(right_4cm_5.Children(2).Children(14).Children(18).Children.Data) str2double(right_4cm_5.Children(2).Children(14).Children(20).Children.Data) str2double(right_4cm_5.Children(2).Children(14).Children(22).Children.Data)];
                   
pr_left_2cm = [str2double(left_2cm_1.Children(2).Children(16).Children(2).Children.Data) str2double(left_2cm_1.Children(2).Children(16).Children(4).Children.Data) str2double(left_2cm_1.Children(2).Children(16).Children(6).Children.Data) str2double(left_2cm_1.Children(2).Children(16).Children(8).Children.Data) str2double(left_2cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_1.Children(2).Children(16).Children(18).Children.Data) str2double(left_2cm_1.Children(2).Children(16).Children(20).Children.Data) str2double(left_2cm_1.Children(2).Children(16).Children(22).Children.Data)
                       str2double(left_2cm_2.Children(2).Children(16).Children(2).Children.Data) str2double(left_2cm_2.Children(2).Children(16).Children(4).Children.Data) str2double(left_2cm_2.Children(2).Children(16).Children(6).Children.Data) str2double(left_2cm_2.Children(2).Children(16).Children(8).Children.Data) str2double(left_2cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_2.Children(2).Children(16).Children(18).Children.Data) str2double(left_2cm_2.Children(2).Children(16).Children(20).Children.Data) str2double(left_2cm_2.Children(2).Children(16).Children(22).Children.Data)
                       str2double(left_2cm_3.Children(2).Children(16).Children(2).Children.Data) str2double(left_2cm_3.Children(2).Children(16).Children(4).Children.Data) str2double(left_2cm_3.Children(2).Children(16).Children(6).Children.Data) str2double(left_2cm_3.Children(2).Children(16).Children(8).Children.Data) str2double(left_2cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_3.Children(2).Children(16).Children(18).Children.Data) str2double(left_2cm_3.Children(2).Children(16).Children(20).Children.Data) str2double(left_2cm_3.Children(2).Children(16).Children(22).Children.Data)
                       str2double(left_2cm_4.Children(2).Children(16).Children(2).Children.Data) str2double(left_2cm_4.Children(2).Children(16).Children(4).Children.Data) str2double(left_2cm_4.Children(2).Children(16).Children(6).Children.Data) str2double(left_2cm_4.Children(2).Children(16).Children(8).Children.Data) str2double(left_2cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_4.Children(2).Children(16).Children(18).Children.Data) str2double(left_2cm_4.Children(2).Children(16).Children(20).Children.Data) str2double(left_2cm_4.Children(2).Children(16).Children(22).Children.Data)
                       str2double(left_2cm_5.Children(2).Children(16).Children(2).Children.Data) str2double(left_2cm_5.Children(2).Children(16).Children(4).Children.Data) str2double(left_2cm_5.Children(2).Children(16).Children(6).Children.Data) str2double(left_2cm_5.Children(2).Children(16).Children(8).Children.Data) str2double(left_2cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm_5.Children(2).Children(16).Children(18).Children.Data) str2double(left_2cm_5.Children(2).Children(16).Children(20).Children.Data) str2double(left_2cm_5.Children(2).Children(16).Children(22).Children.Data)];

pr_right_2cm = [str2double(right_2cm_1.Children(2).Children(16).Children(2).Children.Data) str2double(right_2cm_1.Children(2).Children(16).Children(4).Children.Data) str2double(right_2cm_1.Children(2).Children(16).Children(6).Children.Data) str2double(right_2cm_1.Children(2).Children(16).Children(8).Children.Data) str2double(right_2cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_1.Children(2).Children(16).Children(18).Children.Data) str2double(right_2cm_1.Children(2).Children(16).Children(20).Children.Data) str2double(right_2cm_1.Children(2).Children(16).Children(22).Children.Data)
                       str2double(right_2cm_2.Children(2).Children(16).Children(2).Children.Data) str2double(right_2cm_2.Children(2).Children(16).Children(4).Children.Data) str2double(right_2cm_2.Children(2).Children(16).Children(6).Children.Data) str2double(right_2cm_2.Children(2).Children(16).Children(8).Children.Data) str2double(right_2cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_2.Children(2).Children(16).Children(18).Children.Data) str2double(right_2cm_2.Children(2).Children(16).Children(20).Children.Data) str2double(right_2cm_2.Children(2).Children(16).Children(22).Children.Data)
                       str2double(right_2cm_3.Children(2).Children(16).Children(2).Children.Data) str2double(right_2cm_3.Children(2).Children(16).Children(4).Children.Data) str2double(right_2cm_3.Children(2).Children(16).Children(6).Children.Data) str2double(right_2cm_3.Children(2).Children(16).Children(8).Children.Data) str2double(right_2cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_3.Children(2).Children(16).Children(18).Children.Data) str2double(right_2cm_3.Children(2).Children(16).Children(20).Children.Data) str2double(right_2cm_3.Children(2).Children(16).Children(22).Children.Data)
                       str2double(right_2cm_4.Children(2).Children(16).Children(2).Children.Data) str2double(right_2cm_4.Children(2).Children(16).Children(4).Children.Data) str2double(right_2cm_4.Children(2).Children(16).Children(6).Children.Data) str2double(right_2cm_4.Children(2).Children(16).Children(8).Children.Data) str2double(right_2cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_4.Children(2).Children(16).Children(18).Children.Data) str2double(right_2cm_4.Children(2).Children(16).Children(20).Children.Data) str2double(right_2cm_4.Children(2).Children(16).Children(22).Children.Data)
                       str2double(right_2cm_5.Children(2).Children(16).Children(2).Children.Data) str2double(right_2cm_5.Children(2).Children(16).Children(4).Children.Data) str2double(right_2cm_5.Children(2).Children(16).Children(6).Children.Data) str2double(right_2cm_5.Children(2).Children(16).Children(8).Children.Data) str2double(right_2cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm_5.Children(2).Children(16).Children(18).Children.Data) str2double(right_2cm_5.Children(2).Children(16).Children(20).Children.Data) str2double(right_2cm_5.Children(2).Children(16).Children(22).Children.Data)];
                   
pr_left_4cm = [str2double(left_4cm_1.Children(2).Children(16).Children(2).Children.Data) str2double(left_4cm_1.Children(2).Children(16).Children(4).Children.Data) str2double(left_4cm_1.Children(2).Children(16).Children(6).Children.Data) str2double(left_4cm_1.Children(2).Children(16).Children(8).Children.Data) str2double(left_4cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_1.Children(2).Children(16).Children(18).Children.Data) str2double(left_4cm_1.Children(2).Children(16).Children(20).Children.Data) str2double(left_4cm_1.Children(2).Children(16).Children(22).Children.Data)
                       str2double(left_4cm_2.Children(2).Children(16).Children(2).Children.Data) str2double(left_4cm_2.Children(2).Children(16).Children(4).Children.Data) str2double(left_4cm_2.Children(2).Children(16).Children(6).Children.Data) str2double(left_4cm_2.Children(2).Children(16).Children(8).Children.Data) str2double(left_4cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_2.Children(2).Children(16).Children(18).Children.Data) str2double(left_4cm_2.Children(2).Children(16).Children(20).Children.Data) str2double(left_4cm_2.Children(2).Children(16).Children(22).Children.Data)
                       str2double(left_4cm_3.Children(2).Children(16).Children(2).Children.Data) str2double(left_4cm_3.Children(2).Children(16).Children(4).Children.Data) str2double(left_4cm_3.Children(2).Children(16).Children(6).Children.Data) str2double(left_4cm_3.Children(2).Children(16).Children(8).Children.Data) str2double(left_4cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_3.Children(2).Children(16).Children(18).Children.Data) str2double(left_4cm_3.Children(2).Children(16).Children(20).Children.Data) str2double(left_4cm_3.Children(2).Children(16).Children(22).Children.Data)
                       str2double(left_4cm_4.Children(2).Children(16).Children(2).Children.Data) str2double(left_4cm_4.Children(2).Children(16).Children(4).Children.Data) str2double(left_4cm_4.Children(2).Children(16).Children(6).Children.Data) str2double(left_4cm_4.Children(2).Children(16).Children(8).Children.Data) str2double(left_4cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_4.Children(2).Children(16).Children(18).Children.Data) str2double(left_4cm_4.Children(2).Children(16).Children(20).Children.Data) str2double(left_4cm_4.Children(2).Children(16).Children(22).Children.Data)
                       str2double(left_4cm_5.Children(2).Children(16).Children(2).Children.Data) str2double(left_4cm_5.Children(2).Children(16).Children(4).Children.Data) str2double(left_4cm_5.Children(2).Children(16).Children(6).Children.Data) str2double(left_4cm_5.Children(2).Children(16).Children(8).Children.Data) str2double(left_4cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm_5.Children(2).Children(16).Children(18).Children.Data) str2double(left_4cm_5.Children(2).Children(16).Children(20).Children.Data) str2double(left_4cm_5.Children(2).Children(16).Children(22).Children.Data)];
                   
pr_right_4cm = [str2double(right_4cm_1.Children(2).Children(16).Children(2).Children.Data) str2double(right_4cm_1.Children(2).Children(16).Children(4).Children.Data) str2double(right_4cm_1.Children(2).Children(16).Children(6).Children.Data) str2double(right_4cm_1.Children(2).Children(16).Children(8).Children.Data) str2double(right_4cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_1.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_1.Children(2).Children(16).Children(18).Children.Data) str2double(right_4cm_1.Children(2).Children(16).Children(20).Children.Data) str2double(right_4cm_1.Children(2).Children(16).Children(22).Children.Data)
                       str2double(right_4cm_2.Children(2).Children(16).Children(2).Children.Data) str2double(right_4cm_2.Children(2).Children(16).Children(4).Children.Data) str2double(right_4cm_2.Children(2).Children(16).Children(6).Children.Data) str2double(right_4cm_2.Children(2).Children(16).Children(8).Children.Data) str2double(right_4cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_2.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_2.Children(2).Children(16).Children(18).Children.Data) str2double(right_4cm_2.Children(2).Children(16).Children(20).Children.Data) str2double(right_4cm_2.Children(2).Children(16).Children(22).Children.Data)
                       str2double(right_4cm_3.Children(2).Children(16).Children(2).Children.Data) str2double(right_4cm_3.Children(2).Children(16).Children(4).Children.Data) str2double(right_4cm_3.Children(2).Children(16).Children(6).Children.Data) str2double(right_4cm_3.Children(2).Children(16).Children(8).Children.Data) str2double(right_4cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_3.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_3.Children(2).Children(16).Children(18).Children.Data) str2double(right_4cm_3.Children(2).Children(16).Children(20).Children.Data) str2double(right_4cm_3.Children(2).Children(16).Children(22).Children.Data)
                       str2double(right_4cm_4.Children(2).Children(16).Children(2).Children.Data) str2double(right_4cm_4.Children(2).Children(16).Children(4).Children.Data) str2double(right_4cm_4.Children(2).Children(16).Children(6).Children.Data) str2double(right_4cm_4.Children(2).Children(16).Children(8).Children.Data) str2double(right_4cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_4.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_4.Children(2).Children(16).Children(18).Children.Data) str2double(right_4cm_4.Children(2).Children(16).Children(20).Children.Data) str2double(right_4cm_4.Children(2).Children(16).Children(22).Children.Data)
                       str2double(right_4cm_5.Children(2).Children(16).Children(2).Children.Data) str2double(right_4cm_5.Children(2).Children(16).Children(4).Children.Data) str2double(right_4cm_5.Children(2).Children(16).Children(6).Children.Data) str2double(right_4cm_5.Children(2).Children(16).Children(8).Children.Data) str2double(right_4cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_5.Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm_5.Children(2).Children(16).Children(18).Children.Data) str2double(right_4cm_5.Children(2).Children(16).Children(20).Children.Data) str2double(right_4cm_5.Children(2).Children(16).Children(22).Children.Data)];
                   
% f_cadence = logspace(3,6,31);
% a_cadence_spline = [92.09, 74.16, 58.78, 47.73, 40.43, 34.74, 31.12, 29.12, 28.62, 28.90, 30.55, 33.76, 38.53, 42.90, 40.89, 35.20, 37.46, 33.79, 33.84, 30.91, 28.66, 26.04, 23.49, 21.05, 18.29, 16.33, 14.21, 12.29, 10.86, 9.47, 8.60];
% p_cadence_spline = [71.63, 72.40, 72.98, 69.48, 60.29, 49.89, 39.76, 28.90, 19.02, 10.23, 2.74, -1.92, -2.17, 4.39, 14.53, 14.85, 13.87, 21.55, 22.60, 24.00, 27.12, 31.55, 31.98, 33.74, 34.40, 34.12, 32.88, 29.45, 24.06, 16.37, 5.36];
% a_cadence_bbspice = [232.73, 184.91, 145.82, 117.09, 92.33, 74.70, 62.14, 52.01, 44.03, 39.76, 36.66, 35.30, 35.50, 35.72, 36.28, 36.50, 36.29, 35.58, 34.28, 32.41, 29.93, 26.89, 23.56, 20.30, 17.50, 15.30, 13.73, 12.68, 12.00, 11.58, 11.31];
% p_cadence_bbspice = [-50.36, -41.59, -30.09, -20.34, -11.06, -1.70, 3.98, 9.01, 10.48, 9.69, 8.25, 6.11, 4.70, 4.75, 6.39, 8.57, 11.71, 15.43, 19.60, 24.01, 28.60, 32.60, 35.43, 36.58, 35.73, 32.96, 28.98, 24.56, 20.29, 16.49, 13.31];
                   
% %%% 2.5 cm IED %%%
%                    
% figure
% hAx=axes;
% hAx.XScale='log';
% % xlim([f(4) f(8)]);
% hold all
% for k = 1:5
%     errorbar(f,a_left_2cm(k,:),ar_left_2cm(k,:),'o-','LineWidth',1.5)
%     errorbar(f,a_right_2cm(k,:),ar_right_2cm(k,:),'*--','LineWidth',1.5)
% end
% grid minor
% title('Magnitude - 2.5cm IED')
% xlabel('Frequency [Hz]')
% ylabel('Magnitude [Ohm]')
% legend('left 1','left 2','left 3','left 4','left 5','right 1','right 2','right 3','right 4','right 5')
% 
% figure
% hAx=axes;
% hAx.XScale='log';
% % xlim([f(4) f(8)]);
% hold all
% for k = 1:5
%     errorbar(f,p_left_2cm(k,:),pr_left_2cm(k,:),'o-','LineWidth',1.5)
%     errorbar(f,p_right_2cm(k,:),pr_right_2cm(k,:),'*--','LineWidth',1.5)
% end
% grid minor
% title('Phase - 2.5cm IED')
% xlabel('Frequency [Hz]')
% ylabel('Phase [Degrees]')
% legend('left 1','left 2','left 3','left 4','left 5','right 1','right 2','right 3','right 4','right 5')
% 
% %%% 4cm IED %%%
% 
% figure
% hAx=axes;
% hAx.XScale='log';
% % xlim([f(4) f(8)]);
% hold all
% for k = 1:5
%     errorbar(f,a_left_4cm(k,:),ar_left_4cm(k,:),'o-','LineWidth',1.5)
%     errorbar(f,a_right_4cm(k,:),ar_right_4cm(k,:),'*--','LineWidth',1.5)
% end
% grid minor
% title('Magnitude - 4cm IED')
% xlabel('Frequency [Hz]')
% ylabel('Magnitude [Ohm]')
% legend('left 1','left 2','left 3','left 4','left 5','right 1','right 2','right 3','right 4','right 5')
% 
% figure
% hAx=axes;
% hAx.XScale='log';
% % xlim([f(4) f(8)]);
% hold all
% for k = 1:5
%     errorbar(f,p_left_4cm(k,:),pr_left_4cm(k,:),'o-','LineWidth',1.5)
%     errorbar(f,p_right_4cm(k,:),pr_right_4cm(k,:),'*--','LineWidth',1.5)
% end
% grid minor
% title('Phase - 4cm IED')
% xlabel('Frequency [Hz]')
% ylabel('Phase [Degrees]')
% legend('left 1','left 2','left 3','left 4','left 5','right 1','right 2','right 3','right 4','right 5')

%%% Mean - 2.5 cm IED %%%

if limb == 1 && subject == 1
    figure(1)
    hAx=axes;
    hAx.XScale='log';
    % xlim([f(4) f(8)]);
    % ylim([-10 70]);
    hold on
    set(gca,'FontSize',28);
    set(gcf, 'WindowState', 'maximized');
    errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#D95319')
    errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#EDB120')
    
    plot(f,a_sim_2cm,'*--r','LineWidth',2,'MarkerSize',16)
    plot(f,a_sim_4cm,'*--b','LineWidth',2,'MarkerSize',16)
    plot(f,a_cadence_2cm,'s:r','LineWidth',2,'MarkerSize',16)
    plot(f,a_cadence_4cm,'s:b','LineWidth',2,'MarkerSize',16)
elseif limb == 1 && subject == 2
    figure(1)
    hold on
    errorbar(f,mean([a_left_2cm ; a_right_4cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
    errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
    
%     plot(f,a_sim_2cm,'*--g','LineWidth',2,'MarkerSize',16)
%     plot(f,a_sim_4cm,'*--y','LineWidth',2,'MarkerSize',16)
%     plot(f,a_cadence_2cm,'s:g','LineWidth',2,'MarkerSize',16)
%     plot(f,a_cadence_4cm,'s:y','LineWidth',2,'MarkerSize',16)
    grid minor
    xlabel('Frequency [Hz]')
    ylabel('Magnitude [Ohm]')
    legend({'Subject 1 - Measured Mean - 2.5cm IED',...
            'Subject 1 - Measured Mean - 4cm IED',...
            'Subject 2 - Measured Mean - 2.5cm IED',...
            'Subject 2 - Measured Mean - 4cm IED',...
            'FEA - 2.5cm IED',...
            'FEA - 4cm IED',...
            'Circuit Model - 2.5cm IED',...
            'Circuit Model - 4cm IED',},...
            'Location','southwest','NumColumns',1)
%      print -depsc a_arm_1

elseif limb == 2 && subject == 1
    figure(3)
    hAx=axes;
    hAx.XScale='log';
    % xlim([f(4) f(8)]);
    ylim([-15 60]);
    hold on
    set(gca,'FontSize',28);
    set(gcf, 'WindowState', 'maximized');
    errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#D95319')
    errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#EDB120')
    
    plot(f,a_sim_2cm,'*--r','LineWidth',2,'MarkerSize',16)
    plot(f,a_sim_4cm,'*--b','LineWidth',2,'MarkerSize',16)
    plot(f,a_cadence_2cm,'s:r','LineWidth',2,'MarkerSize',16)
    plot(f,a_cadence_4cm,'s:b','LineWidth',2,'MarkerSize',16)
elseif limb == 2 && subject == 2
    figure(3)
    hold on
    errorbar(f,mean([a_left_2cm ; a_right_4cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
    errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
    
%     plot(f,a_sim_2cm,'*--g','LineWidth',2,'MarkerSize',16)
%     plot(f,a_sim_4cm,'*--y','LineWidth',2,'MarkerSize',16)
%     plot(f,a_cadence_2cm,'s:g','LineWidth',2,'MarkerSize',16)
%     plot(f,a_cadence_4cm,'s:y','LineWidth',2,'MarkerSize',16)
    grid minor
    xlabel('Frequency [Hz]')
    ylabel('Magnitude [Ohm]')
    legend({'Subject 1 - Measured Mean - 2.5cm IED',...
            'Subject 1 - Measured Mean - 4cm IED',...
            'Subject 2 - Measured Mean - 2.5cm IED',...
            'Subject 2 - Measured Mean - 4cm IED',...
            'FEA - 2.5cm IED',...
            'FEA - 4cm IED',...
            'Circuit Model - 2.5cm IED',...
            'Circuit Model - 4cm IED',},...
            'Location','southwest','NumColumns',1)
%      print -depsc a_leg_1
end


if limb == 1 && subject == 1
    figure(2)
    hAx=axes;
    hAx.XScale='log';
    % xlim([f(4) f(8)]);
    % ylim([-10 60]);
    set(gca,'FontSize',28);
    set(gcf, 'WindowState', 'maximized');
    hold on
    errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([p_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#D95319')
    errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([p_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#EDB120')
    
    plot(f,p_sim_2cm,'*--r','LineWidth',2,'MarkerSize',16)
    plot(f,p_sim_4cm,'*--b','LineWidth',2,'MarkerSize',16)
    plot(f,p_cadence_2cm,'s:r','LineWidth',2,'MarkerSize',16)
    plot(f,p_cadence_4cm,'s:b','LineWidth',2,'MarkerSize',16)
elseif limb == 1 && subject == 2
    figure(2)
    hold on
    errorbar(f,mean([p_left_2cm ; p_right_4cm]),std([p_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
    errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([p_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
    
%     plot(f,p_sim_2cm,'*--g','LineWidth',2,'MarkerSize',16)
%     plot(f,p_sim_4cm,'*--y','LineWidth',2,'MarkerSize',16)
%     plot(f,p_cadence_2cm,'s:g','LineWidth',2,'MarkerSize',16)
%     plot(f,p_cadence_4cm,'s:y','LineWidth',2,'MarkerSize',16)
    
    grid minor
    xlabel('Frequency [Hz]')
    ylabel('Phase [Degrees]')
    legend({'Subject 1 - Measured Mean - 2.5cm IED',...
            'Subject 1 - Measured Mean - 4cm IED',...
            'Subject 2 - Measured Mean - 2.5cm IED',...
            'Subject 2 - Measured Mean - 4cm IED',...
            'FEA - 2.5cm IED',...
            'FEA - 4cm IED',...
            'Circuit Model - 2.5cm IED',...
            'Circuit Model - 4cm IED',},...
            'Location','northwest','NumColumns',1)
%      print -depsc p_arm_1
 
elseif limb == 2 && subject == 1
    figure(4)
    hAx=axes;
    hAx.XScale='log';
    % xlim([f(4) f(8)]);
    % ylim([-10 60]);
    set(gca,'FontSize',28);
    set(gcf, 'WindowState', 'maximized');
    hold on
    errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([p_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#D95319')
    errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([p_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#EDB120')
    
    plot(f,p_sim_2cm,'*--r','LineWidth',2,'MarkerSize',16)
    plot(f,p_sim_4cm,'*--b','LineWidth',2,'MarkerSize',16)
    plot(f,p_cadence_2cm,'s:r','LineWidth',2,'MarkerSize',16)
    plot(f,p_cadence_4cm,'s:b','LineWidth',2,'MarkerSize',16)
elseif limb == 2 && subject == 2
    figure(4)
    hold on
    errorbar(f,mean([p_left_2cm ; p_right_4cm]),std([p_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
    errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([p_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
    
%     plot(f,p_sim_2cm,'*--g','LineWidth',2,'MarkerSize',16)
%     plot(f,p_sim_4cm,'*--y','LineWidth',2,'MarkerSize',16)
%     plot(f,p_cadence_2cm,'s:g','LineWidth',2,'MarkerSize',16)
%     plot(f,p_cadence_4cm,'s:y','LineWidth',2,'MarkerSize',16)
    
    grid minor
    xlabel('Frequency [Hz]')
    ylabel('Phase [Degrees]')
    legend({'Subject 1 - Measured Mean - 2.5cm IED',...
            'Subject 1 - Measured Mean - 4cm IED',...
            'Subject 2 - Measured Mean - 2.5cm IED',...
            'Subject 2 - Measured Mean - 4cm IED',...
            'FEA - 2.5cm IED',...
            'FEA - 4cm IED',...
            'Circuit Model - 2.5cm IED',...
            'Circuit Model - 4cm IED',},...
            'Location','northwest','NumColumns',1)
%     print -depsc p_leg_1
end

%%% Errors Bar Plots %%%

% a_meas_2cm = [a_left_2cm ; a_right_2cm];
% a_error_2cm = sqrt(mean((a_meas_2cm - a_sim_2cm).^2));
% p_meas_2cm = [p_left_2cm ; p_right_2cm];
% p_error_2cm = sqrt(mean((p_meas_2cm - p_sim_2cm).^2));
% a_meas_4cm = [a_left_4cm ; a_right_4cm];
% a_error_4cm = sqrt(mean((a_meas_4cm - a_sim_4cm).^2));
% p_meas_4cm = [p_left_4cm ; p_right_4cm];
% p_error_4cm = sqrt(mean((p_meas_4cm - p_sim_4cm).^2));

% figure
% bar(log10(f),[a_error_2cm; a_error_4cm]);
% set(gca,'Xtick',log10(f)); %// adjust manually; values in log scale
% set(gca,'Xticklabel',10.^get(gca,'Xtick')); %// use labels with linear values
% grid minor
% title('Magnitude Error')
% xlabel('frequency [Hz]')
% ylabel('absolute error [\Omega]')
% legend('magnitude - 2.5 cm','magnitude - 4 cm')
% 
% figure
% bar(log10(f),[p_error_2cm; p_error_4cm]);
% set(gca,'Xtick',log10(f)); %// adjust manually; values in log scale
% set(gca,'Xticklabel',10.^get(gca,'Xtick')); %// use labels with linear values
% grid minor
% title('Phase Error')
% xlabel('frequency [Hz]')
% ylabel('absolute error [Degrees]')
% legend('phase - 2.5 cm','phase - 4 cm')
    end
end