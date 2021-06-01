close all
clear all
clc
addpath('/home/alejandro/Documents/PhD Bioelectronics/Experiments/Manuscript1/Two_Subjects_Manuscrip1/8_more_subjects')
fds = fileDatastore('/home/alejandro/Documents/PhD Bioelectronics/Experiments/Manuscript1/Two_Subjects_Manuscrip1/8_more_subjects/*.xml', 'ReadFcn', @importdata);
fds_2sub = fileDatastore('*.xml', 'ReadFcn', @importdata);
fullFileNames = fds.Files;
fullFileNames_2sub = fds_2sub.Files;
numFiles = length(fullFileNames);
numFiles_2sub = length(fullFileNames_2sub);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for a = 1:10
    for b = 1:2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subject = a; % 
limb = b;    % 1 = Arm ; 2 = Leg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f=[976.563	1953.13	3906.25	7812.5	15625	31250	62500	125000	250000	500000	1000000];

if (subject == 1) && (limb == 1) % subject1 - Arms - wf = 0.4 [cm] 11-15
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames_2sub{2 + k});
        left_4cm(k+1)  = parseXML(fullFileNames_2sub{2 + 5*3 + k});
        right_2cm(k+1) = parseXML(fullFileNames_2sub{2 + 5*3*2*2 + k});
        right_4cm(k+1) = parseXML(fullFileNames_2sub{2 + 5*3 + 5*3*2*2 + k});
    end

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
    
elseif (subject == 1) && (limb == 2) % Legs - wf = 0.7 [cm]
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames_2sub{2 + 5*3*2 + k});
        left_4cm(k+1)  = parseXML(fullFileNames_2sub{2 + 5*3 + 5*3*2 + k});
        right_2cm(k+1) = parseXML(fullFileNames_2sub{2 + 5*3*2*2 + 5*3*2 + k});
        right_4cm(k+1) = parseXML(fullFileNames_2sub{2 + 5*3 + 5*3*2*2 + 5*3*2 + k});
    end

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
    
elseif (subject == 2) && (limb == 1) % Arms - wf = 0.4 [cm] 6-10
    for k = 0 : 4
        if k < 4
            left_2cm(k+1)  = parseXML(fullFileNames_2sub{2 + 5*2 + k});
            left_4cm(k+1)  = parseXML(fullFileNames_2sub{2 + 5*2 + 5*3 + k});
            right_2cm(k+1) = parseXML(fullFileNames_2sub{2 + 5*2 + 5*3*2*2 + k});
            right_4cm(k+1) = parseXML(fullFileNames_2sub{2 + 5*2 + 5*3 + 5*3*2*2 + k});
        elseif k == 4
            left_2cm(k+1)  = parseXML(fullFileNames_2sub{2 + k - 5});
            left_4cm(k+1)  = parseXML(fullFileNames_2sub{2 + 5*3 + k -5});
            right_2cm(k+1) = parseXML(fullFileNames_2sub{2 + 5*3*2*2 + k -5});
            right_4cm(k+1) = parseXML(fullFileNames_2sub{2 + 5*3 + 5*3*2*2 + k -5});
        end
    end
    
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
    
elseif (subject == 2) && (limb == 2) % Legs - wf = 1 [cm] 6-10
    for k = 0 : 4
        if k < 4
            left_2cm(k+1)  = parseXML(fullFileNames_2sub{2 + 5*2 + 5*3*2 + k});
            left_4cm(k+1)  = parseXML(fullFileNames_2sub{2 + 5*2 + 5*3 + 5*3*2 + k});
            right_2cm(k+1) = parseXML(fullFileNames_2sub{2 + 5*2 + 5*3*2*2 + 5*3*2 + k});
            right_4cm(k+1) = parseXML(fullFileNames_2sub{2 + 5*2 + 5*3 + 5*3*2*2 + 5*3*2 + k});
        elseif k == 4
            left_2cm(k+1)  = parseXML(fullFileNames_2sub{2 + 5*3*2 + k - 5});
            left_4cm(k+1)  = parseXML(fullFileNames_2sub{2 + 5*3 + 5*3*2 + k -5});
            right_2cm(k+1) = parseXML(fullFileNames_2sub{2 + 5*3*2*2 + 5*3*2 + k -5});
            right_4cm(k+1) = parseXML(fullFileNames_2sub{2 + 5*3 + 5*3*2*2 + 5*3*2 + k -5});
        end
    end
    
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
    
elseif (subject == 3) && (limb == 1) % Yukai - Legs - wf = 1 [cm]
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5});
    end
elseif (subject == 3) && (limb == 2) % Yukai - Legs - wf = 1 [cm]   
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5*2});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5 + 8*5*2});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5*2});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5 + 8*5*2});
    end
elseif (subject == 4) && (limb == 1) % Yukai - Legs - wf = 1 [cm]
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5});
    end
elseif (subject == 4) && (limb == 2) % Yukai - Legs - wf = 1 [cm]   
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5*2});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5 + 8*5*2});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5*2});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5 + 8*5*2});
    end
elseif (subject == 5) && (limb == 1) % Yukai - Legs - wf = 1 [cm]
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5});
    end
elseif (subject == 5) && (limb == 2) % Yukai - Legs - wf = 1 [cm]   
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5*2});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5 + 8*5*2});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5*2});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5 + 8*5*2});
    end
elseif (subject == 6) && (limb == 1) % Yukai - Legs - wf = 1 [cm]
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5});
    end
elseif (subject == 6) && (limb == 2) % Yukai - Legs - wf = 1 [cm]   
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5*2});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5 + 8*5*2});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5*2});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5 + 8*5*2});
    end
elseif (subject == 7) && (limb == 1) % Yukai - Legs - wf = 1 [cm]
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5});
    end
elseif (subject == 7) && (limb == 2) % Yukai - Legs - wf = 1 [cm]   
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5*2});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5 + 8*5*2});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5*2});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5 + 8*5*2});
    end
elseif (subject == 8) && (limb == 1) % Yukai - Legs - wf = 1 [cm]
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5});
    end
elseif (subject == 8) && (limb == 2) % Yukai - Legs - wf = 1 [cm]   
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5*2});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5 + 8*5*2});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5*2});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5 + 8*5*2});
    end
elseif (subject == 9) && (limb == 1) % Yukai - Legs - wf = 1 [cm]
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5});
    end
elseif (subject == 9) && (limb == 2) % Yukai - Legs - wf = 1 [cm]   
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5*2});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5 + 8*5*2});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5*2});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5 + 8*5*2});
    end
elseif (subject == 10) && (limb == 1) % Yukai - Legs - wf = 1 [cm]
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5});
    end
elseif (subject == 10) && (limb == 2) % Yukai - Legs - wf = 1 [cm]   
    for k = 0 : 4
        left_2cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5*2});
        left_4cm(k+1)  = parseXML(fullFileNames{subject-2 + 8*k + 8*5 + 8*5*2});
        right_2cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5*2});
        right_4cm(k+1) = parseXML(fullFileNames{subject-2 + 8*k + 8*5*4 + 8*5 + 8*5*2});
    end
end
 
a_left_2cm =   [str2double(left_2cm(1).Children(2).Children(10).Children(2).Children.Data) str2double(left_2cm(1).Children(2).Children(10).Children(4).Children.Data) str2double(left_2cm(1).Children(2).Children(10).Children(6).Children.Data) str2double(left_2cm(1).Children(2).Children(10).Children(8).Children.Data) str2double(left_2cm(1).Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm(1).Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm(1).Children(2).Children(10).Children(14).Children.Data) str2double(left_2cm(1).Children(2).Children(10).Children(16).Children.Data) str2double(left_2cm(1).Children(2).Children(10).Children(18).Children.Data) str2double(left_2cm(1).Children(2).Children(10).Children(20).Children.Data) str2double(left_2cm(1).Children(2).Children(10).Children(22).Children.Data)
                str2double(left_2cm(2).Children(2).Children(10).Children(2).Children.Data) str2double(left_2cm(2).Children(2).Children(10).Children(4).Children.Data) str2double(left_2cm(2).Children(2).Children(10).Children(6).Children.Data) str2double(left_2cm(2).Children(2).Children(10).Children(8).Children.Data) str2double(left_2cm(2).Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm(2).Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm(2).Children(2).Children(10).Children(14).Children.Data) str2double(left_2cm(2).Children(2).Children(10).Children(16).Children.Data) str2double(left_2cm(2).Children(2).Children(10).Children(18).Children.Data) str2double(left_2cm(2).Children(2).Children(10).Children(20).Children.Data) str2double(left_2cm(2).Children(2).Children(10).Children(22).Children.Data)
                str2double(left_2cm(3).Children(2).Children(10).Children(2).Children.Data) str2double(left_2cm(3).Children(2).Children(10).Children(4).Children.Data) str2double(left_2cm(3).Children(2).Children(10).Children(6).Children.Data) str2double(left_2cm(3).Children(2).Children(10).Children(8).Children.Data) str2double(left_2cm(3).Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm(3).Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm(3).Children(2).Children(10).Children(14).Children.Data) str2double(left_2cm(3).Children(2).Children(10).Children(16).Children.Data) str2double(left_2cm(3).Children(2).Children(10).Children(18).Children.Data) str2double(left_2cm(3).Children(2).Children(10).Children(20).Children.Data) str2double(left_2cm(3).Children(2).Children(10).Children(22).Children.Data)
                str2double(left_2cm(4).Children(2).Children(10).Children(2).Children.Data) str2double(left_2cm(4).Children(2).Children(10).Children(4).Children.Data) str2double(left_2cm(4).Children(2).Children(10).Children(6).Children.Data) str2double(left_2cm(4).Children(2).Children(10).Children(8).Children.Data) str2double(left_2cm(4).Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm(4).Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm(4).Children(2).Children(10).Children(14).Children.Data) str2double(left_2cm(4).Children(2).Children(10).Children(16).Children.Data) str2double(left_2cm(4).Children(2).Children(10).Children(18).Children.Data) str2double(left_2cm(4).Children(2).Children(10).Children(20).Children.Data) str2double(left_2cm(4).Children(2).Children(10).Children(22).Children.Data)
                str2double(left_2cm(5).Children(2).Children(10).Children(2).Children.Data) str2double(left_2cm(5).Children(2).Children(10).Children(4).Children.Data) str2double(left_2cm(5).Children(2).Children(10).Children(6).Children.Data) str2double(left_2cm(5).Children(2).Children(10).Children(8).Children.Data) str2double(left_2cm(5).Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm(5).Children(2).Children(10).Children(10).Children.Data) str2double(left_2cm(5).Children(2).Children(10).Children(14).Children.Data) str2double(left_2cm(5).Children(2).Children(10).Children(16).Children.Data) str2double(left_2cm(5).Children(2).Children(10).Children(18).Children.Data) str2double(left_2cm(5).Children(2).Children(10).Children(20).Children.Data) str2double(left_2cm(5).Children(2).Children(10).Children(22).Children.Data)];

a_right_2cm =  [str2double(right_2cm(1).Children(2).Children(10).Children(2).Children.Data) str2double(right_2cm(1).Children(2).Children(10).Children(4).Children.Data) str2double(right_2cm(1).Children(2).Children(10).Children(6).Children.Data) str2double(right_2cm(1).Children(2).Children(10).Children(8).Children.Data) str2double(right_2cm(1).Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm(1).Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm(1).Children(2).Children(10).Children(14).Children.Data) str2double(right_2cm(1).Children(2).Children(10).Children(16).Children.Data) str2double(right_2cm(1).Children(2).Children(10).Children(18).Children.Data) str2double(right_2cm(1).Children(2).Children(10).Children(20).Children.Data) str2double(right_2cm(1).Children(2).Children(10).Children(22).Children.Data)
                str2double(right_2cm(2).Children(2).Children(10).Children(2).Children.Data) str2double(right_2cm(2).Children(2).Children(10).Children(4).Children.Data) str2double(right_2cm(2).Children(2).Children(10).Children(6).Children.Data) str2double(right_2cm(2).Children(2).Children(10).Children(8).Children.Data) str2double(right_2cm(2).Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm(2).Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm(2).Children(2).Children(10).Children(14).Children.Data) str2double(right_2cm(2).Children(2).Children(10).Children(16).Children.Data) str2double(right_2cm(2).Children(2).Children(10).Children(18).Children.Data) str2double(right_2cm(2).Children(2).Children(10).Children(20).Children.Data) str2double(right_2cm(2).Children(2).Children(10).Children(22).Children.Data)
                str2double(right_2cm(3).Children(2).Children(10).Children(2).Children.Data) str2double(right_2cm(3).Children(2).Children(10).Children(4).Children.Data) str2double(right_2cm(3).Children(2).Children(10).Children(6).Children.Data) str2double(right_2cm(3).Children(2).Children(10).Children(8).Children.Data) str2double(right_2cm(3).Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm(3).Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm(3).Children(2).Children(10).Children(14).Children.Data) str2double(right_2cm(3).Children(2).Children(10).Children(16).Children.Data) str2double(right_2cm(3).Children(2).Children(10).Children(18).Children.Data) str2double(right_2cm(3).Children(2).Children(10).Children(20).Children.Data) str2double(right_2cm(3).Children(2).Children(10).Children(22).Children.Data)
                str2double(right_2cm(4).Children(2).Children(10).Children(2).Children.Data) str2double(right_2cm(4).Children(2).Children(10).Children(4).Children.Data) str2double(right_2cm(4).Children(2).Children(10).Children(6).Children.Data) str2double(right_2cm(4).Children(2).Children(10).Children(8).Children.Data) str2double(right_2cm(4).Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm(4).Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm(4).Children(2).Children(10).Children(14).Children.Data) str2double(right_2cm(4).Children(2).Children(10).Children(16).Children.Data) str2double(right_2cm(4).Children(2).Children(10).Children(18).Children.Data) str2double(right_2cm(4).Children(2).Children(10).Children(20).Children.Data) str2double(right_2cm(4).Children(2).Children(10).Children(22).Children.Data)
                str2double(right_2cm(5).Children(2).Children(10).Children(2).Children.Data) str2double(right_2cm(5).Children(2).Children(10).Children(4).Children.Data) str2double(right_2cm(5).Children(2).Children(10).Children(6).Children.Data) str2double(right_2cm(5).Children(2).Children(10).Children(8).Children.Data) str2double(right_2cm(5).Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm(5).Children(2).Children(10).Children(10).Children.Data) str2double(right_2cm(5).Children(2).Children(10).Children(14).Children.Data) str2double(right_2cm(5).Children(2).Children(10).Children(16).Children.Data) str2double(right_2cm(5).Children(2).Children(10).Children(18).Children.Data) str2double(right_2cm(5).Children(2).Children(10).Children(20).Children.Data) str2double(right_2cm(5).Children(2).Children(10).Children(22).Children.Data)];
                   
a_left_4cm =   [str2double(left_4cm(1).Children(2).Children(10).Children(2).Children.Data) str2double(left_4cm(1).Children(2).Children(10).Children(4).Children.Data) str2double(left_4cm(1).Children(2).Children(10).Children(6).Children.Data) str2double(left_4cm(1).Children(2).Children(10).Children(8).Children.Data) str2double(left_4cm(1).Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm(1).Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm(1).Children(2).Children(10).Children(14).Children.Data) str2double(left_4cm(1).Children(2).Children(10).Children(16).Children.Data) str2double(left_4cm(1).Children(2).Children(10).Children(18).Children.Data) str2double(left_4cm(1).Children(2).Children(10).Children(20).Children.Data) str2double(left_4cm(1).Children(2).Children(10).Children(22).Children.Data)
                str2double(left_4cm(2).Children(2).Children(10).Children(2).Children.Data) str2double(left_4cm(2).Children(2).Children(10).Children(4).Children.Data) str2double(left_4cm(2).Children(2).Children(10).Children(6).Children.Data) str2double(left_4cm(2).Children(2).Children(10).Children(8).Children.Data) str2double(left_4cm(2).Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm(2).Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm(2).Children(2).Children(10).Children(14).Children.Data) str2double(left_4cm(2).Children(2).Children(10).Children(16).Children.Data) str2double(left_4cm(2).Children(2).Children(10).Children(18).Children.Data) str2double(left_4cm(2).Children(2).Children(10).Children(20).Children.Data) str2double(left_4cm(2).Children(2).Children(10).Children(22).Children.Data)
                str2double(left_4cm(3).Children(2).Children(10).Children(2).Children.Data) str2double(left_4cm(3).Children(2).Children(10).Children(4).Children.Data) str2double(left_4cm(3).Children(2).Children(10).Children(6).Children.Data) str2double(left_4cm(3).Children(2).Children(10).Children(8).Children.Data) str2double(left_4cm(3).Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm(3).Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm(3).Children(2).Children(10).Children(14).Children.Data) str2double(left_4cm(3).Children(2).Children(10).Children(16).Children.Data) str2double(left_4cm(3).Children(2).Children(10).Children(18).Children.Data) str2double(left_4cm(3).Children(2).Children(10).Children(20).Children.Data) str2double(left_4cm(3).Children(2).Children(10).Children(22).Children.Data)
                str2double(left_4cm(4).Children(2).Children(10).Children(2).Children.Data) str2double(left_4cm(4).Children(2).Children(10).Children(4).Children.Data) str2double(left_4cm(4).Children(2).Children(10).Children(6).Children.Data) str2double(left_4cm(4).Children(2).Children(10).Children(8).Children.Data) str2double(left_4cm(4).Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm(4).Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm(4).Children(2).Children(10).Children(14).Children.Data) str2double(left_4cm(4).Children(2).Children(10).Children(16).Children.Data) str2double(left_4cm(4).Children(2).Children(10).Children(18).Children.Data) str2double(left_4cm(4).Children(2).Children(10).Children(20).Children.Data) str2double(left_4cm(4).Children(2).Children(10).Children(22).Children.Data)
                str2double(left_4cm(5).Children(2).Children(10).Children(2).Children.Data) str2double(left_4cm(5).Children(2).Children(10).Children(4).Children.Data) str2double(left_4cm(5).Children(2).Children(10).Children(6).Children.Data) str2double(left_4cm(5).Children(2).Children(10).Children(8).Children.Data) str2double(left_4cm(5).Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm(5).Children(2).Children(10).Children(10).Children.Data) str2double(left_4cm(5).Children(2).Children(10).Children(14).Children.Data) str2double(left_4cm(5).Children(2).Children(10).Children(16).Children.Data) str2double(left_4cm(5).Children(2).Children(10).Children(18).Children.Data) str2double(left_4cm(5).Children(2).Children(10).Children(20).Children.Data) str2double(left_4cm(5).Children(2).Children(10).Children(22).Children.Data)];
                   
a_right_4cm =  [str2double(right_4cm(1).Children(2).Children(10).Children(2).Children.Data) str2double(right_4cm(1).Children(2).Children(10).Children(4).Children.Data) str2double(right_4cm(1).Children(2).Children(10).Children(6).Children.Data) str2double(right_4cm(1).Children(2).Children(10).Children(8).Children.Data) str2double(right_4cm(1).Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm(1).Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm(1).Children(2).Children(10).Children(14).Children.Data) str2double(right_4cm(1).Children(2).Children(10).Children(16).Children.Data) str2double(right_4cm(1).Children(2).Children(10).Children(18).Children.Data) str2double(right_4cm(1).Children(2).Children(10).Children(20).Children.Data) str2double(right_4cm(1).Children(2).Children(10).Children(22).Children.Data)
                str2double(right_4cm(2).Children(2).Children(10).Children(2).Children.Data) str2double(right_4cm(2).Children(2).Children(10).Children(4).Children.Data) str2double(right_4cm(2).Children(2).Children(10).Children(6).Children.Data) str2double(right_4cm(2).Children(2).Children(10).Children(8).Children.Data) str2double(right_4cm(2).Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm(2).Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm(2).Children(2).Children(10).Children(14).Children.Data) str2double(right_4cm(2).Children(2).Children(10).Children(16).Children.Data) str2double(right_4cm(2).Children(2).Children(10).Children(18).Children.Data) str2double(right_4cm(2).Children(2).Children(10).Children(20).Children.Data) str2double(right_4cm(2).Children(2).Children(10).Children(22).Children.Data)
                str2double(right_4cm(3).Children(2).Children(10).Children(2).Children.Data) str2double(right_4cm(3).Children(2).Children(10).Children(4).Children.Data) str2double(right_4cm(3).Children(2).Children(10).Children(6).Children.Data) str2double(right_4cm(3).Children(2).Children(10).Children(8).Children.Data) str2double(right_4cm(3).Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm(3).Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm(3).Children(2).Children(10).Children(14).Children.Data) str2double(right_4cm(3).Children(2).Children(10).Children(16).Children.Data) str2double(right_4cm(3).Children(2).Children(10).Children(18).Children.Data) str2double(right_4cm(3).Children(2).Children(10).Children(20).Children.Data) str2double(right_4cm(3).Children(2).Children(10).Children(22).Children.Data)
                str2double(right_4cm(4).Children(2).Children(10).Children(2).Children.Data) str2double(right_4cm(4).Children(2).Children(10).Children(4).Children.Data) str2double(right_4cm(4).Children(2).Children(10).Children(6).Children.Data) str2double(right_4cm(4).Children(2).Children(10).Children(8).Children.Data) str2double(right_4cm(4).Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm(4).Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm(4).Children(2).Children(10).Children(14).Children.Data) str2double(right_4cm(4).Children(2).Children(10).Children(16).Children.Data) str2double(right_4cm(4).Children(2).Children(10).Children(18).Children.Data) str2double(right_4cm(4).Children(2).Children(10).Children(20).Children.Data) str2double(right_4cm(4).Children(2).Children(10).Children(22).Children.Data)
                str2double(right_4cm(5).Children(2).Children(10).Children(2).Children.Data) str2double(right_4cm(5).Children(2).Children(10).Children(4).Children.Data) str2double(right_4cm(5).Children(2).Children(10).Children(6).Children.Data) str2double(right_4cm(5).Children(2).Children(10).Children(8).Children.Data) str2double(right_4cm(5).Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm(5).Children(2).Children(10).Children(10).Children.Data) str2double(right_4cm(5).Children(2).Children(10).Children(14).Children.Data) str2double(right_4cm(5).Children(2).Children(10).Children(16).Children.Data) str2double(right_4cm(5).Children(2).Children(10).Children(18).Children.Data) str2double(right_4cm(5).Children(2).Children(10).Children(20).Children.Data) str2double(right_4cm(5).Children(2).Children(10).Children(22).Children.Data)];
                   
p_left_2cm =   [str2double(left_2cm(1).Children(2).Children(12).Children(2).Children.Data) str2double(left_2cm(1).Children(2).Children(12).Children(4).Children.Data) str2double(left_2cm(1).Children(2).Children(12).Children(6).Children.Data) str2double(left_2cm(1).Children(2).Children(12).Children(8).Children.Data) str2double(left_2cm(1).Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm(1).Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm(1).Children(2).Children(12).Children(14).Children.Data) str2double(left_2cm(1).Children(2).Children(12).Children(16).Children.Data) str2double(left_2cm(1).Children(2).Children(12).Children(18).Children.Data) str2double(left_2cm(1).Children(2).Children(12).Children(20).Children.Data) str2double(left_2cm(1).Children(2).Children(12).Children(22).Children.Data)
                str2double(left_2cm(2).Children(2).Children(12).Children(2).Children.Data) str2double(left_2cm(2).Children(2).Children(12).Children(4).Children.Data) str2double(left_2cm(2).Children(2).Children(12).Children(6).Children.Data) str2double(left_2cm(2).Children(2).Children(12).Children(8).Children.Data) str2double(left_2cm(2).Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm(2).Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm(2).Children(2).Children(12).Children(14).Children.Data) str2double(left_2cm(2).Children(2).Children(12).Children(16).Children.Data) str2double(left_2cm(2).Children(2).Children(12).Children(18).Children.Data) str2double(left_2cm(2).Children(2).Children(12).Children(20).Children.Data) str2double(left_2cm(2).Children(2).Children(12).Children(22).Children.Data)
                str2double(left_2cm(3).Children(2).Children(12).Children(2).Children.Data) str2double(left_2cm(3).Children(2).Children(12).Children(4).Children.Data) str2double(left_2cm(3).Children(2).Children(12).Children(6).Children.Data) str2double(left_2cm(3).Children(2).Children(12).Children(8).Children.Data) str2double(left_2cm(3).Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm(3).Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm(3).Children(2).Children(12).Children(14).Children.Data) str2double(left_2cm(3).Children(2).Children(12).Children(16).Children.Data) str2double(left_2cm(3).Children(2).Children(12).Children(18).Children.Data) str2double(left_2cm(3).Children(2).Children(12).Children(20).Children.Data) str2double(left_2cm(3).Children(2).Children(12).Children(22).Children.Data)
                str2double(left_2cm(4).Children(2).Children(12).Children(2).Children.Data) str2double(left_2cm(4).Children(2).Children(12).Children(4).Children.Data) str2double(left_2cm(4).Children(2).Children(12).Children(6).Children.Data) str2double(left_2cm(4).Children(2).Children(12).Children(8).Children.Data) str2double(left_2cm(4).Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm(4).Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm(4).Children(2).Children(12).Children(14).Children.Data) str2double(left_2cm(4).Children(2).Children(12).Children(16).Children.Data) str2double(left_2cm(4).Children(2).Children(12).Children(18).Children.Data) str2double(left_2cm(4).Children(2).Children(12).Children(20).Children.Data) str2double(left_2cm(4).Children(2).Children(12).Children(22).Children.Data)
                str2double(left_2cm(5).Children(2).Children(12).Children(2).Children.Data) str2double(left_2cm(5).Children(2).Children(12).Children(4).Children.Data) str2double(left_2cm(5).Children(2).Children(12).Children(6).Children.Data) str2double(left_2cm(5).Children(2).Children(12).Children(8).Children.Data) str2double(left_2cm(5).Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm(5).Children(2).Children(12).Children(12).Children.Data) str2double(left_2cm(5).Children(2).Children(12).Children(14).Children.Data) str2double(left_2cm(5).Children(2).Children(12).Children(16).Children.Data) str2double(left_2cm(5).Children(2).Children(12).Children(18).Children.Data) str2double(left_2cm(5).Children(2).Children(12).Children(20).Children.Data) str2double(left_2cm(5).Children(2).Children(12).Children(22).Children.Data)];

p_right_2cm =  [str2double(right_2cm(1).Children(2).Children(12).Children(2).Children.Data) str2double(right_2cm(1).Children(2).Children(12).Children(4).Children.Data) str2double(right_2cm(1).Children(2).Children(12).Children(6).Children.Data) str2double(right_2cm(1).Children(2).Children(12).Children(8).Children.Data) str2double(right_2cm(1).Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm(1).Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm(1).Children(2).Children(12).Children(14).Children.Data) str2double(right_2cm(1).Children(2).Children(12).Children(16).Children.Data) str2double(right_2cm(1).Children(2).Children(12).Children(18).Children.Data) str2double(right_2cm(1).Children(2).Children(12).Children(20).Children.Data) str2double(right_2cm(1).Children(2).Children(12).Children(22).Children.Data)
                str2double(right_2cm(2).Children(2).Children(12).Children(2).Children.Data) str2double(right_2cm(2).Children(2).Children(12).Children(4).Children.Data) str2double(right_2cm(2).Children(2).Children(12).Children(6).Children.Data) str2double(right_2cm(2).Children(2).Children(12).Children(8).Children.Data) str2double(right_2cm(2).Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm(2).Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm(2).Children(2).Children(12).Children(14).Children.Data) str2double(right_2cm(2).Children(2).Children(12).Children(16).Children.Data) str2double(right_2cm(2).Children(2).Children(12).Children(18).Children.Data) str2double(right_2cm(2).Children(2).Children(12).Children(20).Children.Data) str2double(right_2cm(2).Children(2).Children(12).Children(22).Children.Data)
                str2double(right_2cm(3).Children(2).Children(12).Children(2).Children.Data) str2double(right_2cm(3).Children(2).Children(12).Children(4).Children.Data) str2double(right_2cm(3).Children(2).Children(12).Children(6).Children.Data) str2double(right_2cm(3).Children(2).Children(12).Children(8).Children.Data) str2double(right_2cm(3).Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm(3).Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm(3).Children(2).Children(12).Children(14).Children.Data) str2double(right_2cm(3).Children(2).Children(12).Children(16).Children.Data) str2double(right_2cm(3).Children(2).Children(12).Children(18).Children.Data) str2double(right_2cm(3).Children(2).Children(12).Children(20).Children.Data) str2double(right_2cm(3).Children(2).Children(12).Children(22).Children.Data)
                str2double(right_2cm(4).Children(2).Children(12).Children(2).Children.Data) str2double(right_2cm(4).Children(2).Children(12).Children(4).Children.Data) str2double(right_2cm(4).Children(2).Children(12).Children(6).Children.Data) str2double(right_2cm(4).Children(2).Children(12).Children(8).Children.Data) str2double(right_2cm(4).Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm(4).Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm(4).Children(2).Children(12).Children(14).Children.Data) str2double(right_2cm(4).Children(2).Children(12).Children(16).Children.Data) str2double(right_2cm(4).Children(2).Children(12).Children(18).Children.Data) str2double(right_2cm(4).Children(2).Children(12).Children(20).Children.Data) str2double(right_2cm(4).Children(2).Children(12).Children(22).Children.Data)
                str2double(right_2cm(5).Children(2).Children(12).Children(2).Children.Data) str2double(right_2cm(5).Children(2).Children(12).Children(4).Children.Data) str2double(right_2cm(5).Children(2).Children(12).Children(6).Children.Data) str2double(right_2cm(5).Children(2).Children(12).Children(8).Children.Data) str2double(right_2cm(5).Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm(5).Children(2).Children(12).Children(12).Children.Data) str2double(right_2cm(5).Children(2).Children(12).Children(14).Children.Data) str2double(right_2cm(5).Children(2).Children(12).Children(16).Children.Data) str2double(right_2cm(5).Children(2).Children(12).Children(18).Children.Data) str2double(right_2cm(5).Children(2).Children(12).Children(20).Children.Data) str2double(right_2cm(5).Children(2).Children(12).Children(22).Children.Data)];
                   
p_left_4cm =   [str2double(left_4cm(1).Children(2).Children(12).Children(2).Children.Data) str2double(left_4cm(1).Children(2).Children(12).Children(4).Children.Data) str2double(left_4cm(1).Children(2).Children(12).Children(6).Children.Data) str2double(left_4cm(1).Children(2).Children(12).Children(8).Children.Data) str2double(left_4cm(1).Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm(1).Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm(1).Children(2).Children(12).Children(14).Children.Data) str2double(left_4cm(1).Children(2).Children(12).Children(16).Children.Data) str2double(left_4cm(1).Children(2).Children(12).Children(18).Children.Data) str2double(left_4cm(1).Children(2).Children(12).Children(20).Children.Data) str2double(left_4cm(1).Children(2).Children(12).Children(22).Children.Data)
                str2double(left_4cm(2).Children(2).Children(12).Children(2).Children.Data) str2double(left_4cm(2).Children(2).Children(12).Children(4).Children.Data) str2double(left_4cm(2).Children(2).Children(12).Children(6).Children.Data) str2double(left_4cm(2).Children(2).Children(12).Children(8).Children.Data) str2double(left_4cm(2).Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm(2).Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm(2).Children(2).Children(12).Children(14).Children.Data) str2double(left_4cm(2).Children(2).Children(12).Children(16).Children.Data) str2double(left_4cm(2).Children(2).Children(12).Children(18).Children.Data) str2double(left_4cm(2).Children(2).Children(12).Children(20).Children.Data) str2double(left_4cm(2).Children(2).Children(12).Children(22).Children.Data)
                str2double(left_4cm(3).Children(2).Children(12).Children(2).Children.Data) str2double(left_4cm(3).Children(2).Children(12).Children(4).Children.Data) str2double(left_4cm(3).Children(2).Children(12).Children(6).Children.Data) str2double(left_4cm(3).Children(2).Children(12).Children(8).Children.Data) str2double(left_4cm(3).Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm(3).Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm(3).Children(2).Children(12).Children(14).Children.Data) str2double(left_4cm(3).Children(2).Children(12).Children(16).Children.Data) str2double(left_4cm(3).Children(2).Children(12).Children(18).Children.Data) str2double(left_4cm(3).Children(2).Children(12).Children(20).Children.Data) str2double(left_4cm(3).Children(2).Children(12).Children(22).Children.Data)
                str2double(left_4cm(4).Children(2).Children(12).Children(2).Children.Data) str2double(left_4cm(4).Children(2).Children(12).Children(4).Children.Data) str2double(left_4cm(4).Children(2).Children(12).Children(6).Children.Data) str2double(left_4cm(4).Children(2).Children(12).Children(8).Children.Data) str2double(left_4cm(4).Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm(4).Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm(4).Children(2).Children(12).Children(14).Children.Data) str2double(left_4cm(4).Children(2).Children(12).Children(16).Children.Data) str2double(left_4cm(4).Children(2).Children(12).Children(18).Children.Data) str2double(left_4cm(4).Children(2).Children(12).Children(20).Children.Data) str2double(left_4cm(4).Children(2).Children(12).Children(22).Children.Data)
                str2double(left_4cm(5).Children(2).Children(12).Children(2).Children.Data) str2double(left_4cm(5).Children(2).Children(12).Children(4).Children.Data) str2double(left_4cm(5).Children(2).Children(12).Children(6).Children.Data) str2double(left_4cm(5).Children(2).Children(12).Children(8).Children.Data) str2double(left_4cm(5).Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm(5).Children(2).Children(12).Children(12).Children.Data) str2double(left_4cm(5).Children(2).Children(12).Children(14).Children.Data) str2double(left_4cm(5).Children(2).Children(12).Children(16).Children.Data) str2double(left_4cm(5).Children(2).Children(12).Children(18).Children.Data) str2double(left_4cm(5).Children(2).Children(12).Children(20).Children.Data) str2double(left_4cm(5).Children(2).Children(12).Children(22).Children.Data)];
                   
p_right_4cm =  [str2double(right_4cm(1).Children(2).Children(12).Children(2).Children.Data) str2double(right_4cm(1).Children(2).Children(12).Children(4).Children.Data) str2double(right_4cm(1).Children(2).Children(12).Children(6).Children.Data) str2double(right_4cm(1).Children(2).Children(12).Children(8).Children.Data) str2double(right_4cm(1).Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm(1).Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm(1).Children(2).Children(12).Children(14).Children.Data) str2double(right_4cm(1).Children(2).Children(12).Children(16).Children.Data) str2double(right_4cm(1).Children(2).Children(12).Children(18).Children.Data) str2double(right_4cm(1).Children(2).Children(12).Children(20).Children.Data) str2double(right_4cm(1).Children(2).Children(12).Children(22).Children.Data)
                str2double(right_4cm(2).Children(2).Children(12).Children(2).Children.Data) str2double(right_4cm(2).Children(2).Children(12).Children(4).Children.Data) str2double(right_4cm(2).Children(2).Children(12).Children(6).Children.Data) str2double(right_4cm(2).Children(2).Children(12).Children(8).Children.Data) str2double(right_4cm(2).Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm(2).Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm(2).Children(2).Children(12).Children(14).Children.Data) str2double(right_4cm(2).Children(2).Children(12).Children(16).Children.Data) str2double(right_4cm(2).Children(2).Children(12).Children(18).Children.Data) str2double(right_4cm(2).Children(2).Children(12).Children(20).Children.Data) str2double(right_4cm(2).Children(2).Children(12).Children(22).Children.Data)
                str2double(right_4cm(3).Children(2).Children(12).Children(2).Children.Data) str2double(right_4cm(3).Children(2).Children(12).Children(4).Children.Data) str2double(right_4cm(3).Children(2).Children(12).Children(6).Children.Data) str2double(right_4cm(3).Children(2).Children(12).Children(8).Children.Data) str2double(right_4cm(3).Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm(3).Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm(3).Children(2).Children(12).Children(14).Children.Data) str2double(right_4cm(3).Children(2).Children(12).Children(16).Children.Data) str2double(right_4cm(3).Children(2).Children(12).Children(18).Children.Data) str2double(right_4cm(3).Children(2).Children(12).Children(20).Children.Data) str2double(right_4cm(3).Children(2).Children(12).Children(22).Children.Data)
                str2double(right_4cm(4).Children(2).Children(12).Children(2).Children.Data) str2double(right_4cm(4).Children(2).Children(12).Children(4).Children.Data) str2double(right_4cm(4).Children(2).Children(12).Children(6).Children.Data) str2double(right_4cm(4).Children(2).Children(12).Children(8).Children.Data) str2double(right_4cm(4).Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm(4).Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm(4).Children(2).Children(12).Children(14).Children.Data) str2double(right_4cm(4).Children(2).Children(12).Children(16).Children.Data) str2double(right_4cm(4).Children(2).Children(12).Children(18).Children.Data) str2double(right_4cm(4).Children(2).Children(12).Children(20).Children.Data) str2double(right_4cm(4).Children(2).Children(12).Children(22).Children.Data)
                str2double(right_4cm(5).Children(2).Children(12).Children(2).Children.Data) str2double(right_4cm(5).Children(2).Children(12).Children(4).Children.Data) str2double(right_4cm(5).Children(2).Children(12).Children(6).Children.Data) str2double(right_4cm(5).Children(2).Children(12).Children(8).Children.Data) str2double(right_4cm(5).Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm(5).Children(2).Children(12).Children(12).Children.Data) str2double(right_4cm(5).Children(2).Children(12).Children(14).Children.Data) str2double(right_4cm(5).Children(2).Children(12).Children(16).Children.Data) str2double(right_4cm(5).Children(2).Children(12).Children(18).Children.Data) str2double(right_4cm(5).Children(2).Children(12).Children(20).Children.Data) str2double(right_4cm(5).Children(2).Children(12).Children(22).Children.Data)];
                   
ar_left_2cm =  [str2double(left_2cm(1).Children(2).Children(14).Children(2).Children.Data) str2double(left_2cm(1).Children(2).Children(14).Children(4).Children.Data) str2double(left_2cm(1).Children(2).Children(14).Children(6).Children.Data) str2double(left_2cm(1).Children(2).Children(14).Children(8).Children.Data) str2double(left_2cm(1).Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm(1).Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm(1).Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm(1).Children(2).Children(14).Children(16).Children.Data) str2double(left_2cm(1).Children(2).Children(14).Children(18).Children.Data) str2double(left_2cm(1).Children(2).Children(14).Children(20).Children.Data) str2double(left_2cm(1).Children(2).Children(14).Children(22).Children.Data)
                str2double(left_2cm(2).Children(2).Children(14).Children(2).Children.Data) str2double(left_2cm(2).Children(2).Children(14).Children(4).Children.Data) str2double(left_2cm(2).Children(2).Children(14).Children(6).Children.Data) str2double(left_2cm(2).Children(2).Children(14).Children(8).Children.Data) str2double(left_2cm(2).Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm(2).Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm(2).Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm(2).Children(2).Children(14).Children(16).Children.Data) str2double(left_2cm(2).Children(2).Children(14).Children(18).Children.Data) str2double(left_2cm(2).Children(2).Children(14).Children(20).Children.Data) str2double(left_2cm(2).Children(2).Children(14).Children(22).Children.Data)
                str2double(left_2cm(3).Children(2).Children(14).Children(2).Children.Data) str2double(left_2cm(3).Children(2).Children(14).Children(4).Children.Data) str2double(left_2cm(3).Children(2).Children(14).Children(6).Children.Data) str2double(left_2cm(3).Children(2).Children(14).Children(8).Children.Data) str2double(left_2cm(3).Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm(3).Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm(3).Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm(3).Children(2).Children(14).Children(16).Children.Data) str2double(left_2cm(3).Children(2).Children(14).Children(18).Children.Data) str2double(left_2cm(3).Children(2).Children(14).Children(20).Children.Data) str2double(left_2cm(3).Children(2).Children(14).Children(22).Children.Data)
                str2double(left_2cm(4).Children(2).Children(14).Children(2).Children.Data) str2double(left_2cm(4).Children(2).Children(14).Children(4).Children.Data) str2double(left_2cm(4).Children(2).Children(14).Children(6).Children.Data) str2double(left_2cm(4).Children(2).Children(14).Children(8).Children.Data) str2double(left_2cm(4).Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm(4).Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm(4).Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm(4).Children(2).Children(14).Children(16).Children.Data) str2double(left_2cm(4).Children(2).Children(14).Children(18).Children.Data) str2double(left_2cm(4).Children(2).Children(14).Children(20).Children.Data) str2double(left_2cm(4).Children(2).Children(14).Children(22).Children.Data)
                str2double(left_2cm(5).Children(2).Children(14).Children(2).Children.Data) str2double(left_2cm(5).Children(2).Children(14).Children(4).Children.Data) str2double(left_2cm(5).Children(2).Children(14).Children(6).Children.Data) str2double(left_2cm(5).Children(2).Children(14).Children(8).Children.Data) str2double(left_2cm(5).Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm(5).Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm(5).Children(2).Children(14).Children(14).Children.Data) str2double(left_2cm(5).Children(2).Children(14).Children(16).Children.Data) str2double(left_2cm(5).Children(2).Children(14).Children(18).Children.Data) str2double(left_2cm(5).Children(2).Children(14).Children(20).Children.Data) str2double(left_2cm(5).Children(2).Children(14).Children(22).Children.Data)];

ar_right_2cm = [str2double(right_2cm(1).Children(2).Children(14).Children(2).Children.Data) str2double(right_2cm(1).Children(2).Children(14).Children(4).Children.Data) str2double(right_2cm(1).Children(2).Children(14).Children(6).Children.Data) str2double(right_2cm(1).Children(2).Children(14).Children(8).Children.Data) str2double(right_2cm(1).Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm(1).Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm(1).Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm(1).Children(2).Children(14).Children(16).Children.Data) str2double(right_2cm(1).Children(2).Children(14).Children(18).Children.Data) str2double(right_2cm(1).Children(2).Children(14).Children(20).Children.Data) str2double(right_2cm(1).Children(2).Children(14).Children(22).Children.Data)
                str2double(right_2cm(2).Children(2).Children(14).Children(2).Children.Data) str2double(right_2cm(2).Children(2).Children(14).Children(4).Children.Data) str2double(right_2cm(2).Children(2).Children(14).Children(6).Children.Data) str2double(right_2cm(2).Children(2).Children(14).Children(8).Children.Data) str2double(right_2cm(2).Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm(2).Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm(2).Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm(2).Children(2).Children(14).Children(16).Children.Data) str2double(right_2cm(2).Children(2).Children(14).Children(18).Children.Data) str2double(right_2cm(2).Children(2).Children(14).Children(20).Children.Data) str2double(right_2cm(2).Children(2).Children(14).Children(22).Children.Data)
                str2double(right_2cm(3).Children(2).Children(14).Children(2).Children.Data) str2double(right_2cm(3).Children(2).Children(14).Children(4).Children.Data) str2double(right_2cm(3).Children(2).Children(14).Children(6).Children.Data) str2double(right_2cm(3).Children(2).Children(14).Children(8).Children.Data) str2double(right_2cm(3).Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm(3).Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm(3).Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm(3).Children(2).Children(14).Children(16).Children.Data) str2double(right_2cm(3).Children(2).Children(14).Children(18).Children.Data) str2double(right_2cm(3).Children(2).Children(14).Children(20).Children.Data) str2double(right_2cm(3).Children(2).Children(14).Children(22).Children.Data)
                str2double(right_2cm(4).Children(2).Children(14).Children(2).Children.Data) str2double(right_2cm(4).Children(2).Children(14).Children(4).Children.Data) str2double(right_2cm(4).Children(2).Children(14).Children(6).Children.Data) str2double(right_2cm(4).Children(2).Children(14).Children(8).Children.Data) str2double(right_2cm(4).Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm(4).Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm(4).Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm(4).Children(2).Children(14).Children(16).Children.Data) str2double(right_2cm(4).Children(2).Children(14).Children(18).Children.Data) str2double(right_2cm(4).Children(2).Children(14).Children(20).Children.Data) str2double(right_2cm(4).Children(2).Children(14).Children(22).Children.Data)
                str2double(right_2cm(5).Children(2).Children(14).Children(2).Children.Data) str2double(right_2cm(5).Children(2).Children(14).Children(4).Children.Data) str2double(right_2cm(5).Children(2).Children(14).Children(6).Children.Data) str2double(right_2cm(5).Children(2).Children(14).Children(8).Children.Data) str2double(right_2cm(5).Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm(5).Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm(5).Children(2).Children(14).Children(14).Children.Data) str2double(right_2cm(5).Children(2).Children(14).Children(16).Children.Data) str2double(right_2cm(5).Children(2).Children(14).Children(18).Children.Data) str2double(right_2cm(5).Children(2).Children(14).Children(20).Children.Data) str2double(right_2cm(5).Children(2).Children(14).Children(22).Children.Data)];
                   
ar_left_4cm =  [str2double(left_4cm(1).Children(2).Children(14).Children(2).Children.Data) str2double(left_4cm(1).Children(2).Children(14).Children(4).Children.Data) str2double(left_4cm(1).Children(2).Children(14).Children(6).Children.Data) str2double(left_4cm(1).Children(2).Children(14).Children(8).Children.Data) str2double(left_4cm(1).Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm(1).Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm(1).Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm(1).Children(2).Children(14).Children(16).Children.Data) str2double(left_4cm(1).Children(2).Children(14).Children(18).Children.Data) str2double(left_4cm(1).Children(2).Children(14).Children(20).Children.Data) str2double(left_4cm(1).Children(2).Children(14).Children(22).Children.Data)
                str2double(left_4cm(2).Children(2).Children(14).Children(2).Children.Data) str2double(left_4cm(2).Children(2).Children(14).Children(4).Children.Data) str2double(left_4cm(2).Children(2).Children(14).Children(6).Children.Data) str2double(left_4cm(2).Children(2).Children(14).Children(8).Children.Data) str2double(left_4cm(2).Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm(2).Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm(2).Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm(2).Children(2).Children(14).Children(16).Children.Data) str2double(left_4cm(2).Children(2).Children(14).Children(18).Children.Data) str2double(left_4cm(2).Children(2).Children(14).Children(20).Children.Data) str2double(left_4cm(2).Children(2).Children(14).Children(22).Children.Data)
                str2double(left_4cm(3).Children(2).Children(14).Children(2).Children.Data) str2double(left_4cm(3).Children(2).Children(14).Children(4).Children.Data) str2double(left_4cm(3).Children(2).Children(14).Children(6).Children.Data) str2double(left_4cm(3).Children(2).Children(14).Children(8).Children.Data) str2double(left_4cm(3).Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm(3).Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm(3).Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm(3).Children(2).Children(14).Children(16).Children.Data) str2double(left_4cm(3).Children(2).Children(14).Children(18).Children.Data) str2double(left_4cm(3).Children(2).Children(14).Children(20).Children.Data) str2double(left_4cm(3).Children(2).Children(14).Children(22).Children.Data)
                str2double(left_4cm(4).Children(2).Children(14).Children(2).Children.Data) str2double(left_4cm(4).Children(2).Children(14).Children(4).Children.Data) str2double(left_4cm(4).Children(2).Children(14).Children(6).Children.Data) str2double(left_4cm(4).Children(2).Children(14).Children(8).Children.Data) str2double(left_4cm(4).Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm(4).Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm(4).Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm(4).Children(2).Children(14).Children(16).Children.Data) str2double(left_4cm(4).Children(2).Children(14).Children(18).Children.Data) str2double(left_4cm(4).Children(2).Children(14).Children(20).Children.Data) str2double(left_4cm(4).Children(2).Children(14).Children(22).Children.Data)
                str2double(left_4cm(5).Children(2).Children(14).Children(2).Children.Data) str2double(left_4cm(5).Children(2).Children(14).Children(4).Children.Data) str2double(left_4cm(5).Children(2).Children(14).Children(6).Children.Data) str2double(left_4cm(5).Children(2).Children(14).Children(8).Children.Data) str2double(left_4cm(5).Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm(5).Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm(5).Children(2).Children(14).Children(14).Children.Data) str2double(left_4cm(5).Children(2).Children(14).Children(16).Children.Data) str2double(left_4cm(5).Children(2).Children(14).Children(18).Children.Data) str2double(left_4cm(5).Children(2).Children(14).Children(20).Children.Data) str2double(left_4cm(5).Children(2).Children(14).Children(22).Children.Data)];
                   
ar_right_4cm = [str2double(right_4cm(1).Children(2).Children(14).Children(2).Children.Data) str2double(right_4cm(1).Children(2).Children(14).Children(4).Children.Data) str2double(right_4cm(1).Children(2).Children(14).Children(6).Children.Data) str2double(right_4cm(1).Children(2).Children(14).Children(8).Children.Data) str2double(right_4cm(1).Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm(1).Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm(1).Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm(1).Children(2).Children(14).Children(16).Children.Data) str2double(right_4cm(1).Children(2).Children(14).Children(18).Children.Data) str2double(right_4cm(1).Children(2).Children(14).Children(20).Children.Data) str2double(right_4cm(1).Children(2).Children(14).Children(22).Children.Data)
                str2double(right_4cm(2).Children(2).Children(14).Children(2).Children.Data) str2double(right_4cm(2).Children(2).Children(14).Children(4).Children.Data) str2double(right_4cm(2).Children(2).Children(14).Children(6).Children.Data) str2double(right_4cm(2).Children(2).Children(14).Children(8).Children.Data) str2double(right_4cm(2).Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm(2).Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm(2).Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm(2).Children(2).Children(14).Children(16).Children.Data) str2double(right_4cm(2).Children(2).Children(14).Children(18).Children.Data) str2double(right_4cm(2).Children(2).Children(14).Children(20).Children.Data) str2double(right_4cm(2).Children(2).Children(14).Children(22).Children.Data)
                str2double(right_4cm(3).Children(2).Children(14).Children(2).Children.Data) str2double(right_4cm(3).Children(2).Children(14).Children(4).Children.Data) str2double(right_4cm(3).Children(2).Children(14).Children(6).Children.Data) str2double(right_4cm(3).Children(2).Children(14).Children(8).Children.Data) str2double(right_4cm(3).Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm(3).Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm(3).Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm(3).Children(2).Children(14).Children(16).Children.Data) str2double(right_4cm(3).Children(2).Children(14).Children(18).Children.Data) str2double(right_4cm(3).Children(2).Children(14).Children(20).Children.Data) str2double(right_4cm(3).Children(2).Children(14).Children(22).Children.Data)
                str2double(right_4cm(4).Children(2).Children(14).Children(2).Children.Data) str2double(right_4cm(4).Children(2).Children(14).Children(4).Children.Data) str2double(right_4cm(4).Children(2).Children(14).Children(6).Children.Data) str2double(right_4cm(4).Children(2).Children(14).Children(8).Children.Data) str2double(right_4cm(4).Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm(4).Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm(4).Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm(4).Children(2).Children(14).Children(16).Children.Data) str2double(right_4cm(4).Children(2).Children(14).Children(18).Children.Data) str2double(right_4cm(4).Children(2).Children(14).Children(20).Children.Data) str2double(right_4cm(4).Children(2).Children(14).Children(22).Children.Data)
                str2double(right_4cm(5).Children(2).Children(14).Children(2).Children.Data) str2double(right_4cm(5).Children(2).Children(14).Children(4).Children.Data) str2double(right_4cm(5).Children(2).Children(14).Children(6).Children.Data) str2double(right_4cm(5).Children(2).Children(14).Children(8).Children.Data) str2double(right_4cm(5).Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm(5).Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm(5).Children(2).Children(14).Children(14).Children.Data) str2double(right_4cm(5).Children(2).Children(14).Children(16).Children.Data) str2double(right_4cm(5).Children(2).Children(14).Children(18).Children.Data) str2double(right_4cm(5).Children(2).Children(14).Children(20).Children.Data) str2double(right_4cm(5).Children(2).Children(14).Children(22).Children.Data)];
                   
pr_left_2cm =  [str2double(left_2cm(1).Children(2).Children(16).Children(2).Children.Data) str2double(left_2cm(1).Children(2).Children(16).Children(4).Children.Data) str2double(left_2cm(1).Children(2).Children(16).Children(6).Children.Data) str2double(left_2cm(1).Children(2).Children(16).Children(8).Children.Data) str2double(left_2cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(1).Children(2).Children(16).Children(18).Children.Data) str2double(left_2cm(1).Children(2).Children(16).Children(20).Children.Data) str2double(left_2cm(1).Children(2).Children(16).Children(22).Children.Data)
                str2double(left_2cm(2).Children(2).Children(16).Children(2).Children.Data) str2double(left_2cm(2).Children(2).Children(16).Children(4).Children.Data) str2double(left_2cm(2).Children(2).Children(16).Children(6).Children.Data) str2double(left_2cm(2).Children(2).Children(16).Children(8).Children.Data) str2double(left_2cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(2).Children(2).Children(16).Children(18).Children.Data) str2double(left_2cm(2).Children(2).Children(16).Children(20).Children.Data) str2double(left_2cm(2).Children(2).Children(16).Children(22).Children.Data)
                str2double(left_2cm(3).Children(2).Children(16).Children(2).Children.Data) str2double(left_2cm(3).Children(2).Children(16).Children(4).Children.Data) str2double(left_2cm(3).Children(2).Children(16).Children(6).Children.Data) str2double(left_2cm(3).Children(2).Children(16).Children(8).Children.Data) str2double(left_2cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(3).Children(2).Children(16).Children(18).Children.Data) str2double(left_2cm(3).Children(2).Children(16).Children(20).Children.Data) str2double(left_2cm(3).Children(2).Children(16).Children(22).Children.Data)
                str2double(left_2cm(4).Children(2).Children(16).Children(2).Children.Data) str2double(left_2cm(4).Children(2).Children(16).Children(4).Children.Data) str2double(left_2cm(4).Children(2).Children(16).Children(6).Children.Data) str2double(left_2cm(4).Children(2).Children(16).Children(8).Children.Data) str2double(left_2cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(4).Children(2).Children(16).Children(18).Children.Data) str2double(left_2cm(4).Children(2).Children(16).Children(20).Children.Data) str2double(left_2cm(4).Children(2).Children(16).Children(22).Children.Data)
                str2double(left_2cm(5).Children(2).Children(16).Children(2).Children.Data) str2double(left_2cm(5).Children(2).Children(16).Children(4).Children.Data) str2double(left_2cm(5).Children(2).Children(16).Children(6).Children.Data) str2double(left_2cm(5).Children(2).Children(16).Children(8).Children.Data) str2double(left_2cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(left_2cm(5).Children(2).Children(16).Children(18).Children.Data) str2double(left_2cm(5).Children(2).Children(16).Children(20).Children.Data) str2double(left_2cm(5).Children(2).Children(16).Children(22).Children.Data)];

pr_right_2cm = [str2double(right_2cm(1).Children(2).Children(16).Children(2).Children.Data) str2double(right_2cm(1).Children(2).Children(16).Children(4).Children.Data) str2double(right_2cm(1).Children(2).Children(16).Children(6).Children.Data) str2double(right_2cm(1).Children(2).Children(16).Children(8).Children.Data) str2double(right_2cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(1).Children(2).Children(16).Children(18).Children.Data) str2double(right_2cm(1).Children(2).Children(16).Children(20).Children.Data) str2double(right_2cm(1).Children(2).Children(16).Children(22).Children.Data)
                str2double(right_2cm(2).Children(2).Children(16).Children(2).Children.Data) str2double(right_2cm(2).Children(2).Children(16).Children(4).Children.Data) str2double(right_2cm(2).Children(2).Children(16).Children(6).Children.Data) str2double(right_2cm(2).Children(2).Children(16).Children(8).Children.Data) str2double(right_2cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(2).Children(2).Children(16).Children(18).Children.Data) str2double(right_2cm(2).Children(2).Children(16).Children(20).Children.Data) str2double(right_2cm(2).Children(2).Children(16).Children(22).Children.Data)
                str2double(right_2cm(3).Children(2).Children(16).Children(2).Children.Data) str2double(right_2cm(3).Children(2).Children(16).Children(4).Children.Data) str2double(right_2cm(3).Children(2).Children(16).Children(6).Children.Data) str2double(right_2cm(3).Children(2).Children(16).Children(8).Children.Data) str2double(right_2cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(3).Children(2).Children(16).Children(18).Children.Data) str2double(right_2cm(3).Children(2).Children(16).Children(20).Children.Data) str2double(right_2cm(3).Children(2).Children(16).Children(22).Children.Data)
                str2double(right_2cm(4).Children(2).Children(16).Children(2).Children.Data) str2double(right_2cm(4).Children(2).Children(16).Children(4).Children.Data) str2double(right_2cm(4).Children(2).Children(16).Children(6).Children.Data) str2double(right_2cm(4).Children(2).Children(16).Children(8).Children.Data) str2double(right_2cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(4).Children(2).Children(16).Children(18).Children.Data) str2double(right_2cm(4).Children(2).Children(16).Children(20).Children.Data) str2double(right_2cm(4).Children(2).Children(16).Children(22).Children.Data)
                str2double(right_2cm(5).Children(2).Children(16).Children(2).Children.Data) str2double(right_2cm(5).Children(2).Children(16).Children(4).Children.Data) str2double(right_2cm(5).Children(2).Children(16).Children(6).Children.Data) str2double(right_2cm(5).Children(2).Children(16).Children(8).Children.Data) str2double(right_2cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(right_2cm(5).Children(2).Children(16).Children(18).Children.Data) str2double(right_2cm(5).Children(2).Children(16).Children(20).Children.Data) str2double(right_2cm(5).Children(2).Children(16).Children(22).Children.Data)];
                   
pr_left_4cm =  [str2double(left_4cm(1).Children(2).Children(16).Children(2).Children.Data) str2double(left_4cm(1).Children(2).Children(16).Children(4).Children.Data) str2double(left_4cm(1).Children(2).Children(16).Children(6).Children.Data) str2double(left_4cm(1).Children(2).Children(16).Children(8).Children.Data) str2double(left_4cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(1).Children(2).Children(16).Children(18).Children.Data) str2double(left_4cm(1).Children(2).Children(16).Children(20).Children.Data) str2double(left_4cm(1).Children(2).Children(16).Children(22).Children.Data)
                str2double(left_4cm(2).Children(2).Children(16).Children(2).Children.Data) str2double(left_4cm(2).Children(2).Children(16).Children(4).Children.Data) str2double(left_4cm(2).Children(2).Children(16).Children(6).Children.Data) str2double(left_4cm(2).Children(2).Children(16).Children(8).Children.Data) str2double(left_4cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(2).Children(2).Children(16).Children(18).Children.Data) str2double(left_4cm(2).Children(2).Children(16).Children(20).Children.Data) str2double(left_4cm(2).Children(2).Children(16).Children(22).Children.Data)
                str2double(left_4cm(3).Children(2).Children(16).Children(2).Children.Data) str2double(left_4cm(3).Children(2).Children(16).Children(4).Children.Data) str2double(left_4cm(3).Children(2).Children(16).Children(6).Children.Data) str2double(left_4cm(3).Children(2).Children(16).Children(8).Children.Data) str2double(left_4cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(3).Children(2).Children(16).Children(18).Children.Data) str2double(left_4cm(3).Children(2).Children(16).Children(20).Children.Data) str2double(left_4cm(3).Children(2).Children(16).Children(22).Children.Data)
                str2double(left_4cm(4).Children(2).Children(16).Children(2).Children.Data) str2double(left_4cm(4).Children(2).Children(16).Children(4).Children.Data) str2double(left_4cm(4).Children(2).Children(16).Children(6).Children.Data) str2double(left_4cm(4).Children(2).Children(16).Children(8).Children.Data) str2double(left_4cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(4).Children(2).Children(16).Children(18).Children.Data) str2double(left_4cm(4).Children(2).Children(16).Children(20).Children.Data) str2double(left_4cm(4).Children(2).Children(16).Children(22).Children.Data)
                str2double(left_4cm(5).Children(2).Children(16).Children(2).Children.Data) str2double(left_4cm(5).Children(2).Children(16).Children(4).Children.Data) str2double(left_4cm(5).Children(2).Children(16).Children(6).Children.Data) str2double(left_4cm(5).Children(2).Children(16).Children(8).Children.Data) str2double(left_4cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(left_4cm(5).Children(2).Children(16).Children(18).Children.Data) str2double(left_4cm(5).Children(2).Children(16).Children(20).Children.Data) str2double(left_4cm(5).Children(2).Children(16).Children(22).Children.Data)];
                   
pr_right_4cm = [str2double(right_4cm(1).Children(2).Children(16).Children(2).Children.Data) str2double(right_4cm(1).Children(2).Children(16).Children(4).Children.Data) str2double(right_4cm(1).Children(2).Children(16).Children(6).Children.Data) str2double(right_4cm(1).Children(2).Children(16).Children(8).Children.Data) str2double(right_4cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(1).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(1).Children(2).Children(16).Children(18).Children.Data) str2double(right_4cm(1).Children(2).Children(16).Children(20).Children.Data) str2double(right_4cm(1).Children(2).Children(16).Children(22).Children.Data)
                str2double(right_4cm(2).Children(2).Children(16).Children(2).Children.Data) str2double(right_4cm(2).Children(2).Children(16).Children(4).Children.Data) str2double(right_4cm(2).Children(2).Children(16).Children(6).Children.Data) str2double(right_4cm(2).Children(2).Children(16).Children(8).Children.Data) str2double(right_4cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(2).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(2).Children(2).Children(16).Children(18).Children.Data) str2double(right_4cm(2).Children(2).Children(16).Children(20).Children.Data) str2double(right_4cm(2).Children(2).Children(16).Children(22).Children.Data)
                str2double(right_4cm(3).Children(2).Children(16).Children(2).Children.Data) str2double(right_4cm(3).Children(2).Children(16).Children(4).Children.Data) str2double(right_4cm(3).Children(2).Children(16).Children(6).Children.Data) str2double(right_4cm(3).Children(2).Children(16).Children(8).Children.Data) str2double(right_4cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(3).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(3).Children(2).Children(16).Children(18).Children.Data) str2double(right_4cm(3).Children(2).Children(16).Children(20).Children.Data) str2double(right_4cm(3).Children(2).Children(16).Children(22).Children.Data)
                str2double(right_4cm(4).Children(2).Children(16).Children(2).Children.Data) str2double(right_4cm(4).Children(2).Children(16).Children(4).Children.Data) str2double(right_4cm(4).Children(2).Children(16).Children(6).Children.Data) str2double(right_4cm(4).Children(2).Children(16).Children(8).Children.Data) str2double(right_4cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(4).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(4).Children(2).Children(16).Children(18).Children.Data) str2double(right_4cm(4).Children(2).Children(16).Children(20).Children.Data) str2double(right_4cm(4).Children(2).Children(16).Children(22).Children.Data)
                str2double(right_4cm(5).Children(2).Children(16).Children(2).Children.Data) str2double(right_4cm(5).Children(2).Children(16).Children(4).Children.Data) str2double(right_4cm(5).Children(2).Children(16).Children(6).Children.Data) str2double(right_4cm(5).Children(2).Children(16).Children(8).Children.Data) str2double(right_4cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(5).Children(2).Children(16).Children(16).Children.Data) str2double(right_4cm(5).Children(2).Children(16).Children(18).Children.Data) str2double(right_4cm(5).Children(2).Children(16).Children(20).Children.Data) str2double(right_4cm(5).Children(2).Children(16).Children(22).Children.Data)];

a_2cm(:,subject,limb) = mean([a_left_2cm ; a_right_2cm]);
a_4cm(:,subject,limb) = mean([a_left_4cm ; a_right_4cm]);
p_2cm(:,subject,limb) = mean([p_left_2cm ; p_right_2cm]);
p_4cm(:,subject,limb) = mean([p_left_4cm ; p_right_4cm]);
ar_2cm(:,subject,limb) = sqrt(mean([ar_left_2cm.^2 ; ar_right_2cm.^2]));
ar_4cm(:,subject,limb) = sqrt(mean([ar_left_4cm.^2 ; ar_right_4cm.^2]));
pr_2cm(:,subject,limb) = sqrt(mean([pr_left_2cm.^2 ; pr_right_2cm.^2]));
pr_4cm(:,subject,limb) = sqrt(mean([pr_left_4cm.^2 ; pr_right_4cm.^2]));

if subject == 10
    rmse_FEA_pha_2cm(limb) = sqrt(mean((mean(p_2cm(:,:,limb),2) - p_sim_2cm).^2));
    rmse_FEA_pha_4cm(limb) = sqrt(mean((mean(p_4cm(:,:,limb),2) - p_sim_4cm).^2));
    rmse_FEA_mag_2cm(limb) = sqrt(mean((mean(a_2cm(:,:,limb),2) - a_sim_2cm).^2));
    rmse_FEA_mag_4cm(limb) = sqrt(mean((mean(a_4cm(:,:,limb),2) - a_sim_4cm).^2));
    rmse_circ_pha_2cm(limb) = sqrt(mean((mean(p_2cm(:,:,limb),2) - p_cadence_2cm).^2));
    rmse_circ_pha_4cm(limb) = sqrt(mean((mean(p_4cm(:,:,limb),2) - p_cadence_4cm).^2));
    rmse_circ_mag_2cm(limb) = sqrt(mean((mean(a_2cm(:,:,limb),2) - a_cadence_2cm).^2));
    rmse_circ_mag_4cm(limb) = sqrt(mean((mean(a_4cm(:,:,limb),2) - a_cadence_4cm).^2));
end

if limb == 1 && subject == 2
    % simulated magnitude
    figure(1)
    hAx=axes;
    hAx.XScale='log';
    hold on
    set(gca,'FontSize',20);
    set(gcf, 'WindowState', 'maximized')
    plot(log10(f),a_sim_2cm,'*--r','LineWidth',2,'MarkerSize',16)
    plot(log10(f),a_sim_4cm,'*--b','LineWidth',2,'MarkerSize',16)
    plot(log10(f),a_cadence_2cm,'s:r','LineWidth',2,'MarkerSize',16)
    plot(log10(f),a_cadence_4cm,'s:b','LineWidth',2,'MarkerSize',16)
    
    % simulated phase
    figure(2)
    hAx=axes;
    hAx.XScale='log';
    hold on
    set(gca,'FontSize',20);
    set(gcf, 'WindowState', 'maximized')
    plot(log10(f),p_sim_2cm,'*--r','LineWidth',2,'MarkerSize',16)
    plot(log10(f),p_sim_4cm,'*--b','LineWidth',2,'MarkerSize',16)
    plot(log10(f),p_cadence_2cm,'s:r','LineWidth',2,'MarkerSize',16)
    plot(log10(f),p_cadence_4cm,'s:b','LineWidth',2,'MarkerSize',16)
elseif limb == 1 && subject == 10
    a_2cm_4cm(:,1,:) = a_2cm(:,:,limb);
    a_2cm_4cm(:,2,:) = a_4cm(:,:,limb);
    p_2cm_4cm(:,1,:) = p_2cm(:,:,limb);
    p_2cm_4cm(:,2,:) = p_4cm(:,:,limb);
    
    % boxplot magnitude
    figure(1)
    box_a = boxplot2(a_2cm_4cm,log10(f));
    cmap = get(0, 'defaultaxescolororder');
    for ii = 1:2
        structfun(@(x) set(x(3-ii,:), 'color', cmap(ii,:), ...
            'markeredgecolor', cmap(ii,:)), box_a);
    end
    ylim([0 130])
    xticks(log10(f))
    xticklabels(num2str(f'/1e3,'%.0f'))
    grid minor
    xlabel('Frequency [kHz]')
    ylabel('Magnitude [Ohm]')
    legend({'FEA - 2.5cm IED',...
            'FEA - 4cm IED',...
            'Circuit Model - 2.5cm IED',...
            'Circuit Model - 4cm IED',...
            'Measurements - 2.5cm IED',...
            'Measurements - 4cm IED'},...
            'Location','northeast','NumColumns',1)
%      print -depsc a_arm_1
    
    % boxplot phase
    figure(2)
    box_a = boxplot2(p_2cm_4cm,log10(f));
    cmap = get(0, 'defaultaxescolororder');
    for ii = 1:2
        structfun(@(x) set(x(3-ii,:), 'color', cmap(ii,:), ...
            'markeredgecolor', cmap(ii,:)), box_a);
    end
    xticks(log10(f))
    xticklabels(num2str(f'/1e3,'%.0f'))
    grid minor
    xlabel('Frequency [kHz]')
    ylabel('Phase [Degrees]')
    legend({'FEA - 2.5cm IED',...
            'FEA - 4cm IED',...
            'Circuit Model - 2.5cm IED',...
            'Circuit Model - 4cm IED',...
            'Measurements - 2.5cm IED',...
            'Measurements - 4cm IED'},...
            'Location','southeast','NumColumns',1)
%      print -depsc p_arm_1
end


if limb == 2 && subject == 2
    % simulated magnitude
    figure(3)
    hAx=axes;
    hAx.XScale='log';
    hold on
    set(gca,'FontSize',20);
    set(gcf, 'WindowState', 'maximized')
    plot(log10(f),a_sim_2cm,'*--r','LineWidth',2,'MarkerSize',16)
    plot(log10(f),a_sim_4cm,'*--b','LineWidth',2,'MarkerSize',16)
    plot(log10(f),a_cadence_2cm,'s:r','LineWidth',2,'MarkerSize',16)
    plot(log10(f),a_cadence_4cm,'s:b','LineWidth',2,'MarkerSize',16)
    % simulated phase
    figure(4)
    hAx=axes;
    hAx.XScale='log';
    hold on
    set(gca,'FontSize',20);
    set(gcf, 'WindowState', 'maximized')
    plot(log10(f),p_sim_2cm,'*--r','LineWidth',2,'MarkerSize',16)
    plot(log10(f),p_sim_4cm,'*--b','LineWidth',2,'MarkerSize',16)
    plot(log10(f),p_cadence_2cm,'s:r','LineWidth',2,'MarkerSize',16)
    plot(log10(f),p_cadence_4cm,'s:b','LineWidth',2,'MarkerSize',16)
elseif limb == 2 && subject == 10
    a_2cm_4cm(:,1,:) = a_2cm(:,:,limb);
    a_2cm_4cm(:,2,:) = a_4cm(:,:,limb);
    p_2cm_4cm(:,1,:) = p_2cm(:,:,limb);
    p_2cm_4cm(:,2,:) = p_4cm(:,:,limb);
    
    % boxplot magnitude
    figure(3)
    box_a = boxplot2(a_2cm_4cm,log10(f));
    cmap = get(0, 'defaultaxescolororder');
    for ii = 1:2
        structfun(@(x) set(x(3-ii,:), 'color', cmap(ii,:), ...
            'markeredgecolor', cmap(ii,:)), box_a);
    end
    xticks(log10(f))
    xticklabels(num2str(f'/1e3,'%.0f'))
    grid minor
    xlabel('Frequency [kHz]')
    ylabel('Magnitude [Ohm]')
    legend({'FEA - 2.5cm IED',...
            'FEA - 4cm IED',...
            'Circuit Model - 2.5cm IED',...
            'Circuit Model - 4cm IED',...
            'Measurements - 2.5cm IED',...
            'Measurements - 4cm IED'},...
            'Location','northeast','NumColumns',1)
%      print -depsc a_leg_1
    
    % boxplot phase
    figure(4)
    box_a = boxplot2(p_2cm_4cm,log10(f));
    cmap = get(0, 'defaultaxescolororder');
    for ii = 1:2
        structfun(@(x) set(x(3-ii,:), 'color', cmap(ii,:), ...
            'markeredgecolor', cmap(ii,:)), box_a);
    end
    xticks(log10(f))
    xticklabels(num2str(f'/1e3,'%.0f'))
    grid minor
    xlabel('Frequency [kHz]')
    ylabel('Phase [Degrees]')
    legend({'FEA - 2.5cm IED',...
            'FEA - 4cm IED',...
            'Circuit Model - 2.5cm IED',...
            'Circuit Model - 4cm IED',...
            'Measurements - 2.5cm IED',...
            'Measurements - 4cm IED'},...
            'Location','southeast','NumColumns',1)
%      print -depsc p_leg_1
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Mean - 2.5 cm IED %%%
% 
% if limb == 1 && subject == 1
%     figure(1)
%     hAx=axes;
%     hAx.XScale='log';
%     % xlim([f(4) f(8)]);
%     % ylim([-10 70]);
%     hold on
%     set(gca,'FontSize',28);
%     set(gcf, 'WindowState', 'maximized');
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#D95319')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#EDB120')
% elseif limb == 1 && subject == 3
%     figure(1)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 1 && subject == 4
%     figure(1)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 1 && subject == 5
%     figure(1)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 1 && subject == 6
%     figure(1)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 1 && subject == 7
%     figure(1)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 1 && subject == 8
%     figure(1)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 1 && subject == 9
%     figure(1)
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 1 && subject == 10
%     figure(1)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% 
% elseif limb == 1 && subject == 2
%     figure(1)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
%     plot(f,a_sim_2cm,'*--r','LineWidth',2,'MarkerSize',16)
%     plot(f,a_sim_4cm,'*--b','LineWidth',2,'MarkerSize',16)
%     plot(f,a_cadence_2cm,'s:r','LineWidth',2,'MarkerSize',16)
%     plot(f,a_cadence_4cm,'s:b','LineWidth',2,'MarkerSize',16)
%     
% %     plot(f,a_sim_2cm,'*--g','LineWidth',2,'MarkerSize',16)
% %     plot(f,a_sim_4cm,'*--y','LineWidth',2,'MarkerSize',16)
% %     plot(f,a_cadence_2cm,'s:g','LineWidth',2,'MarkerSize',16)
% %     plot(f,a_cadence_4cm,'s:y','LineWidth',2,'MarkerSize',16)
%     grid minor
%     xlabel('Frequency [Hz]')
%     ylabel('Magnitude [Ohm]')
%     legend({'Subject 1 - Measured Mean - 2.5cm IED',...
%             'Subject 1 - Measured Mean - 4cm IED',...
%             'Subject 2 - Measured Mean - 2.5cm IED',...
%             'Subject 2 - Measured Mean - 4cm IED',...
%             'FEA - 2.5cm IED',...
%             'FEA - 4cm IED',...
%             'Circuit Model - 2.5cm IED',...
%             'Circuit Model - 4cm IED',},...
%             'Location','southwest','NumColumns',1)
% %      print -depsc a_arm_1
% 
% elseif limb == 2 && subject == 1
%     figure(3)
%     hAx=axes;
%     hAx.XScale='log';
%     % xlim([f(4) f(8)]);
%     ylim([-15 60]);
%     hold on
%     set(gca,'FontSize',28);
%     set(gcf, 'WindowState', 'maximized');
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#D95319')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#EDB120')
% elseif limb == 2 && subject == 3
%     figure(3)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 4
%     figure(3)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 5
%     figure(3)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 6
%     figure(3)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 7
%     figure(3)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 8
%     figure(3)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 9
%     figure(3)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 10
%     figure(3)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 2
%     figure(3)
%     hold on
%     errorbar(f,mean([a_left_2cm ; a_right_2cm]),std([a_left_2cm ; a_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([a_left_4cm ; a_right_4cm]),std([a_left_4cm ; a_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
%     plot(f,a_sim_2cm,'*--r','LineWidth',2,'MarkerSize',16)
%     plot(f,a_sim_4cm,'*--b','LineWidth',2,'MarkerSize',16)
%     plot(f,a_cadence_2cm,'s:r','LineWidth',2,'MarkerSize',16)
%     plot(f,a_cadence_4cm,'s:b','LineWidth',2,'MarkerSize',16)
%     
% %     plot(f,a_sim_2cm,'*--g','LineWidth',2,'MarkerSize',16)
% %     plot(f,a_sim_4cm,'*--y','LineWidth',2,'MarkerSize',16)
% %     plot(f,a_cadence_2cm,'s:g','LineWidth',2,'MarkerSize',16)
% %     plot(f,a_cadence_4cm,'s:y','LineWidth',2,'MarkerSize',16)
%     grid minor
%     xlabel('Frequency [Hz]')
%     ylabel('Magnitude [Ohm]')
%     legend({'Subject 1 - Measured Mean - 2.5cm IED',...
%             'Subject 1 - Measured Mean - 4cm IED',...
%             'Subject 2 - Measured Mean - 2.5cm IED',...
%             'Subject 2 - Measured Mean - 4cm IED',...
%             'FEA - 2.5cm IED',...
%             'FEA - 4cm IED',...
%             'Circuit Model - 2.5cm IED',...
%             'Circuit Model - 4cm IED',},...
%             'Location','southwest','NumColumns',1)
% %      print -depsc a_leg_1
% end
% 
% 
% if limb == 1 && subject == 1
%     figure(2)
%     hAx=axes;
%     hAx.XScale='log';
%     % xlim([f(4) f(8)]);
%     % ylim([-10 60]);
%     set(gca,'FontSize',28);
%     set(gcf, 'WindowState', 'maximized');
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([p_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#D95319')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([p_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#EDB120')
% elseif limb == 1 && subject == 3
%     figure(2)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 1 && subject == 4
%     figure(2)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 1 && subject == 5
%     figure(2)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 1 && subject == 6
%     figure(2)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 1 && subject == 7
%     figure(2)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 1 && subject == 8
%     figure(2)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 1 && subject == 9
%     figure(2)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 1 && subject == 10
%     figure(2)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 1 && subject == 2
%     figure(2)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([p_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([p_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
%     plot(f,p_sim_2cm,'*--r','LineWidth',2,'MarkerSize',16)
%     plot(f,p_sim_4cm,'*--b','LineWidth',2,'MarkerSize',16)
%     plot(f,p_cadence_2cm,'s:r','LineWidth',2,'MarkerSize',16)
%     plot(f,p_cadence_4cm,'s:b','LineWidth',2,'MarkerSize',16)
%     
% %     plot(f,p_sim_2cm,'*--g','LineWidth',2,'MarkerSize',16)
% %     plot(f,p_sim_4cm,'*--y','LineWidth',2,'MarkerSize',16)
% %     plot(f,p_cadence_2cm,'s:g','LineWidth',2,'MarkerSize',16)
% %     plot(f,p_cadence_4cm,'s:y','LineWidth',2,'MarkerSize',16)
%     
%     grid minor
%     xlabel('Frequency [Hz]')
%     ylabel('Phase [Degrees]')
%     legend({'Subject 1 - Measured Mean - 2.5cm IED',...
%             'Subject 1 - Measured Mean - 4cm IED',...
%             'Subject 2 - Measured Mean - 2.5cm IED',...
%             'Subject 2 - Measured Mean - 4cm IED',...
%             'FEA - 2.5cm IED',...
%             'FEA - 4cm IED',...
%             'Circuit Model - 2.5cm IED',...
%             'Circuit Model - 4cm IED',},...
%             'Location','northwest','NumColumns',1)
% %      print -depsc p_arm_1
%  
% elseif limb == 2 && subject == 1
%     figure(4)
%     hAx=axes;
%     hAx.XScale='log';
%     % xlim([f(4) f(8)]);
%     % ylim([-10 60]);
%     set(gca,'FontSize',28);
%     set(gcf, 'WindowState', 'maximized');
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([p_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#D95319')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([p_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#EDB120')
% elseif limb == 2 && subject == 3
%     figure(4)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 4
%     figure(4)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 5
%     figure(4)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 6
%     figure(4)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 7
%     figure(4)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 8
%     figure(4)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 9
%     figure(4)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 10
%     figure(4)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([a_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([a_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
% elseif limb == 2 && subject == 2
%     figure(4)
%     hold on
%     errorbar(f,mean([p_left_2cm ; p_right_2cm]),std([p_left_2cm ; p_right_2cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#7E2F8E')
%     errorbar(f,mean([p_left_4cm ; p_right_4cm]),std([p_left_4cm ; p_right_4cm]),'o-','LineWidth',2,'MarkerSize',16,'Color','#77AC30')
%     plot(f,p_sim_2cm,'*--r','LineWidth',2,'MarkerSize',16)
%     plot(f,p_sim_4cm,'*--b','LineWidth',2,'MarkerSize',16)
%     plot(f,p_cadence_2cm,'s:r','LineWidth',2,'MarkerSize',16)
%     plot(f,p_cadence_4cm,'s:b','LineWidth',2,'MarkerSize',16)
%     
% %     plot(f,p_sim_2cm,'*--g','LineWidth',2,'MarkerSize',16)
% %     plot(f,p_sim_4cm,'*--y','LineWidth',2,'MarkerSize',16)
% %     plot(f,p_cadence_2cm,'s:g','LineWidth',2,'MarkerSize',16)
% %     plot(f,p_cadence_4cm,'s:y','LineWidth',2,'MarkerSize',16)
%     
%     grid minor
%     xlabel('Frequency [Hz]')
%     ylabel('Phase [Degrees]')
%     legend({'Subject 1 - Measured Mean - 2.5cm IED',...
%             'Subject 1 - Measured Mean - 4cm IED',...
%             'Subject 2 - Measured Mean - 2.5cm IED',...
%             'Subject 2 - Measured Mean - 4cm IED',...
%             'FEA - 2.5cm IED',...
%             'FEA - 4cm IED',...
%             'Circuit Model - 2.5cm IED',...
%             'Circuit Model - 4cm IED',},...
%             'Location','northwest','NumColumns',1)
% %     print -depsc p_leg_1
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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