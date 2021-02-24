clc
clear
close all

% bioz_cont = csvread('contributions_tissue_impedance.csv',5);
bioz_cont_arm = csvread('sel_tissue_impedance_arm_2.csv',5);
bioz_cont_leg = csvread('sel_tissue_impedance_leg_2.csv',5);

freq = bioz_cont_arm(11:end,1);

stackData = zeros(11,5,2);
% cont_1 = real(bioz_cont(11:end,2:6));
% cont_2 = real(bioz_cont(11:end,17:21));
% stackData(:,:,1) = cont_1;
% stackData(:,:,2) = cont_2;
fc = hsv(10);

% plotBarStackGroups(stackData,freq)

figure('WindowState','maximized')
% set(gca,'FontSize',30);
% xlim([log10(freq(5)-2.5e3) log10(freq(9)+5e4)]);
hold all
stem(log10(freq),bioz_cont_arm(11:end,2),'o','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[255/255 204/255 204/255],'LineWidth',2);
stem(log10(freq),bioz_cont_arm(11:end,3),'o','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[255/255 255/255 0/255],'LineWidth',2);
stem(log10(freq),bioz_cont_arm(11:end,4),'o','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[255/255 0/255 0/255],'LineWidth',2);
stem(log10(freq),bioz_cont_arm(11:end,5),'o','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[160/255 160/255 160/255],'LineWidth',2);
stem(log10(freq),bioz_cont_arm(11:end,6),'o','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[204/255 102/255 0/255],'LineWidth',2);

stem(log10(freq),bioz_cont_arm(11:end,7),'d','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[255/255 102/255 102/255],'LineWidth',2);
stem(log10(freq),bioz_cont_arm(11:end,8),'d','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[153/255 153/255 0/255],'LineWidth',2);
stem(log10(freq),bioz_cont_arm(11:end,9),'d','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[153/255 0/255 0/255],'LineWidth',2);
stem(log10(freq),bioz_cont_arm(11:end,10),'d','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[96/255 96/255 96/255],'LineWidth',2);
stem(log10(freq),bioz_cont_arm(11:end,11),'d','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[102/255 51/255 0/255],'LineWidth',2);

% stem(log10(freq),bioz_cont(11:end,k),'o','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fc(k-11,:));
xticks(log10(freq))
xticklabels(string(freq))
grid minor
title('Percentage Contribution of Tissues - UpperArm')
xlabel('Frequency [Hz]')
ylabel('Percentage [%]')
legend('Skin - 2.5 cm','Fat - 2.5 cm','Muscle - 2.5 cm','Bone - 2.5 cm','Marrow - 2.5 cm','Skin - 4 cm','Fat - 4 cm','Muscle - 4 cm','Bone - 4 cm','Marrow - 4 cm','Location','northeastoutside')
% print -depsc cont_arm

figure('WindowState','maximized')
% set(gca,'FontSize',30);
% xlim([log10(freq(5)-2.5e3) log10(freq(9)+5e4)]);
hold all
stem(log10(freq),bioz_cont_leg(11:end,2),'o','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[255/255 204/255 204/255],'LineWidth',2);
stem(log10(freq),bioz_cont_leg(11:end,3),'o','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[255/255 255/255 0/255],'LineWidth',2);
stem(log10(freq),bioz_cont_leg(11:end,4),'o','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[255/255 0/255 0/255],'LineWidth',2);
stem(log10(freq),bioz_cont_leg(11:end,5),'o','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[160/255 160/255 160/255],'LineWidth',2);
stem(log10(freq),bioz_cont_leg(11:end,6),'o','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[204/255 102/255 0/255],'LineWidth',2);

stem(log10(freq),bioz_cont_leg(11:end,7),'d','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[255/255 102/255 102/255],'LineWidth',2);
stem(log10(freq),bioz_cont_leg(11:end,8),'d','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[153/255 153/255 0/255],'LineWidth',2);
stem(log10(freq),bioz_cont_leg(11:end,9),'d','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[153/255 0/255 0/255],'LineWidth',2);
stem(log10(freq),bioz_cont_leg(11:end,10),'d','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[96/255 96/255 96/255],'LineWidth',2);
stem(log10(freq),bioz_cont_leg(11:end,11),'d','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[102/255 51/255 0/255],'LineWidth',2);

% stem(log10(freq),bioz_cont(11:end,k),'o','MarkerSize',20,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fc(k-11,:));
xticks(log10(freq))
xticklabels(string(freq))
grid minor
title('Percentage Contribution of Tissues - Lower Leg')
xlabel('frequency [Hz]')
ylabel('Percentage [%]')
legend('Skin - 2.5 cm','Fat - 2.5 cm','Muscle - 2.5 cm','Bone - 2.5 cm','Marrow - 2.5 cm','Skin - 4 cm','Fat - 4 cm','Muscle - 4 cm','Bone - 4 cm','Marrow - 4 cm','Location','northeastoutside')
% print -depsc cont_leg
