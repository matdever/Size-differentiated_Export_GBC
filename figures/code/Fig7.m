clear

%% This figure plots the histogram of the relative biomass vs vertical velocity at 
% DOY 25 in Summer/Winter for xi = 2 and xi = 4

%% Import w for each size class in winter and summer

% Load data previously extracted from big database
load dataset/Fig7_data.mat

%% Set colors for each class
color0025 = [215,48,39]/256;
color1 = [39,100,25]/256;
color5 = [49,54,149]/256;

%% compute biomass per simulated particle

% Set of marameters
wmax = 5; wmin = 0.025;
slope = 3; B0 = 1;
% Compute total biomass based on the paramters
Btot = B0/(wmax.^((3-slope)/2))*(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));

% For Junge slope = 2
slope = 2;
% Compute biomass of 1 largest particle necessary to keep total biomass constant
B0 = Btot*(wmax.^((3-slope)/2))/(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));
B0025_2 = B0*(0.025/wmax)^((3-slope)/2);
B1_2 = B0*(1/wmax)^((3-slope)/2);
B5_2 = B0*(5/wmax)^((3-slope)/2);

% For Junge slope = 4
slope = 4;
% Compute biomass of 1 largest particle necessary to keep total biomass constant
B0 = Btot*(wmax.^((3-slope)/2))/(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));
B0025_4 = B0*(0.025/wmax)^((3-slope)/2);
B1_4 = B0*(1/wmax)^((3-slope)/2);
B5_4 = B0*(5/wmax)^((3-slope)/2);

%%
% Calculate total biomass to normalize by
totalB = sum(N0025s.*B0025_2)+sum(N1s.*B1_2)+sum(N5s.*B5_2);
totalBw = sum(N0025w.*B0025_2)+sum(N1w.*B1_2)+sum(N5w.*B5_2);
totalB4 = sum(N0025s.*B0025_4)+sum(N1s.*B1_4)+sum(N5s.*B5_4);
totalB4w = sum(N0025w.*B0025_4)+sum(N1w.*B1_4)+sum(N5w.*B5_4);
totalF = sum(N0025s(E<=0).*B0025_2.*E(E<=0))+sum(N1s(E<=0).*B1_2.*E(E<=0))+sum(N5s(E<=0).*B5_2.*E(E<=0));
totalF4 = sum(N0025s(E<=0).*B0025_4.*E(E<=0))+sum(N1s(E<=0).*B1_4.*E(E<=0))+sum(N5s(E<=0).*B5_4.*E(E<=0));
totalFw = sum(N0025w(Ew<=0).*B0025_2.*Ew(Ew<=0))+sum(N1w(Ew<=0).*B1_2.*Ew(Ew<=0))+sum(N5w(Ew<=0).*B5_2.*Ew(Ew<=0));
totalF4w = sum(N0025w(Ew<=0).*B0025_4.*Ew(Ew<=0))+sum(N1w(Ew<=0).*B1_4.*Ew(Ew<=0))+sum(N5w(Ew<=0).*B5_4.*Ew(Ew<=0));

figure('position',[15 150 2119 1361])
subplot(2,2,1)
% Plot PDFs
h1 = plot(E,N0025s.*B0025_2/totalB*100,'-','color',color0025,'linewidth',2);
hold on
h2 = plot(E,N1s.*B1_2/totalB*100,'-','color',color1,'linewidth',2);
h3 = plot(E,N5s.*B5_2/totalB*100,'-','color',color5,'linewidth',2);
% Shade negative flux
patch([E(E<=0) 0 min(E(E<=0)) min(E(E<=0))],[N0025s(E<=0).*B0025_2/totalB*100 1e-2 1e-2 N0025s(1).*B0025_2/totalB*100],color0025,'edgecolor','none','facealpha',.3)
patch([E(E<=0) 0 min(E(E<=0)) min(E(E<=0))],[N1s(E<=0).*B1_2/totalB*100 1e-2 1e-2 N1s(1).*B1_2/totalB*100],color1,'edgecolor','none','facealpha',.3)
patch([E(E<=0) 0 min(E(E<=0)) min(E(E<=0))],[N5s(E<=0).*B5_2/totalB*100 1e-2 1e-2 N5s(1).*B5_2/totalB*100],color5,'edgecolor','none','facealpha',.3)

set(gca,'Yscale','log','FontSize',28)
%xlabel('w [m/day]')
ylabel('Relative Biomass [%]')
grid on; axis tight; xlim([-10 10])
ylim([1e-2 100])
title('Summer','fontsize',50,'position',[-6 11e+01 0]);
text(-9.5,6e1,'(a)   \xi = 2','fontweight','b','FontSize',28)
HH = legend([h1 h2 h3],'0.025 m/day','1 m/day','5 m/day','location','southeast');
set(HH,'fontsize',30,'linewidth',2)
%[HH,icons,plots,txt] = legend([h1 h2 h3],'0.025 m/day','1 m/day','5 m/day','location','southeast');
%set(plots,'linewidth',15)


axes('position',[.36 .82 .15 .14])
bar(1,sum(N0025s(E<=0).*B0025_2.*E(E<=0)/totalF*100),'facecolor',color0025,'edgecolor',color0025)
text(1,sum(N0025s(E<=0).*B0025_2.*E(E<=0)/totalF*100),[num2str(sum(N0025s(E<=0).*B0025_2.*E(E<=0)/totalF*100),'%2.1f'),'%'],'horizontalalignment','center','verticalalignment','top','fontsize',17,'fontweight','b')
hold on
bar(2,sum(N1s(E<=0).*B1_2.*E(E<=0)/totalF*100),'facecolor',color1,'edgecolor',color1);
text(2,sum(N1s(E<=0).*B1_2.*E(E<=0)/totalF*100),[num2str(sum(N1s(E<=0).*B1_2.*E(E<=0)/totalF*100),'%2.1f'),'%'],'horizontalalignment','center','verticalalignment','top','fontsize',17,'fontweight','b')
bar(3,sum(N5s(E<=0).*B5_2.*E(E<=0)/totalF*100),'facecolor',color5,'edgecolor',color5);
text(3,sum(N5s(E<=0).*B5_2.*E(E<=0)/totalF*100),[num2str(sum(N5s(E<=0).*B5_2.*E(E<=0)/totalF*100),'%2.1f'),'%'],'horizontalalignment','center','verticalalignment','top','fontsize',17,'fontweight','b','color','w')
set(gca,'YScale','log','FontSize',28,'xtick',[1:3],'XTickLabel',{},'yaxislocation','right')
grid on; box on
ylim([1e-2 100])
title('Relative Downward Flux','FontSize',28,'position',[2 200 0])

subplot(2,2,3)
h1 = plot(E,N0025s.*B0025_4/totalB4*100,'-','color',color0025,'linewidth',2);
hold on
plot(E,N1s.*B1_4/totalB4*100,'-','color',color1,'linewidth',2)
plot(E,N5s.*B5_4/totalB4*100,'-','color',color5,'linewidth',2)
% Shade negative flux
patch([E(E<=0) 0 min(E(E<=0)) min(E(E<=0))],[N0025s(E<=0).*B0025_4/totalB4*100 1e-2 1e-2 N0025s(1).*B0025_4/totalB4*100],color0025,'edgecolor','none','facealpha',.3)
patch([E(E<=0) 0 min(E(E<=0)) min(E(E<=0))],[N1s(E<=0).*B1_4/totalB4*100 1e-2 1e-2 N1s(1).*B1_4/totalB4*100],color1,'edgecolor','none','facealpha',.3)
patch([E(E<=0) 0 min(E(E<=0)) min(E(E<=0))],[N5s(E<=0).*B5_4/totalB4*100 1e-2 1e-2 N5s(1).*B5_4/totalB4*100],color5,'edgecolor','none','facealpha',.3)

set(gca,'Yscale','log','FontSize',28)
xlabel('w_t_o_t [m/day]')
ylabel('Relative Biomass [%]')
grid on; axis tight; xlim([-10 10])
ylim([1e-2 100])
text(-9.5,6e1,'(c)   \xi = 4','fontweight','b','FontSize',28)

axes('position',[.36 .34 .15 .14])
bar(1,sum(N0025s(E<=0).*B0025_4.*E(E<=0))/totalF4*100,'facecolor',color0025,'edgecolor',color0025)
text(1,sum(N0025s(E<=0).*B0025_4.*E(E<=0))/totalF4*100,[num2str(sum(N0025s(E<=0).*B0025_4.*E(E<=0))/totalF4*100,'%2.1f'),'%'],'horizontalalignment','center','verticalalignment','top','fontsize',17,'fontweight','b')
hold on
bar(2,sum(N1s(E<=0).*B1_4.*E(E<=0))/totalF4*100,'facecolor',color1,'edgecolor',color1);
text(2,sum(N1s(E<=0).*B1_4.*E(E<=0))/totalF4*100,[num2str(sum(N1s(E<=0).*B1_4.*E(E<=0))/totalF4*100,'%2.1f'),'%'],'horizontalalignment','center','verticalalignment','top','fontsize',17,'fontweight','b')
bar(3,sum(N5s(E<=0).*B5_4.*E(E<=0))/totalF4*100,'facecolor',color5,'edgecolor',color5);
text(3,sum(N5s(E<=0).*B5_4.*E(E<=0))/totalF4*100,[num2str(sum(N5s(E<=0).*B5_4.*E(E<=0))/totalF4*100,'%2.1f'),'%'],'horizontalalignment','center','verticalalignment','top','fontsize',17,'fontweight','b','color','w')
set(gca,'YScale','log','FontSize',28,'xtick',[1:3],'XTickLabel',{},'yaxislocation','right')
grid on; box on
ylim([1e-2 100])
title('Relative Downward Flux','FontSize',28,'position',[2 200 0])

subplot(2,2,2)
h1 = plot(Ew,N0025w.*B0025_2/totalBw*100,'-','color',color0025,'linewidth',2);
hold on
plot(Ew,N1w.*B1_2/totalBw*100,'-','color',color1,'linewidth',2)
plot(Ew,N5w.*B5_2/totalBw*100,'-','color',color5,'linewidth',2)
% Shade negative flux
patch([Ew(Ew<=0) 0 min(Ew(Ew<=0)) min(Ew(Ew<=0))],[N0025w(Ew<=0).*B0025_2/totalBw*100 1e-2 1e-2 N0025w(1).*B0025_2/totalBw*100],color0025,'edgecolor','none','facealpha',.3)
patch([Ew(Ew<=0) 0 min(Ew(Ew<=0)) min(Ew(Ew<=0))],[N1w(Ew<=0).*B1_2/totalBw*100 1e-2 1e-2 N1w(1).*B1_2/totalBw*100],color1,'edgecolor','none','facealpha',.3)
patch([Ew(Ew<=0) 0 min(Ew(Ew<=0)) min(Ew(Ew<=0))],[N5w(Ew<=0).*B5_2/totalBw*100 1e-2 1e-2 N5w(1).*B5_2/totalBw*100],color5,'edgecolor','none','facealpha',.3)

set(gca,'Yscale','log','FontSize',28,'yticklabel',{})
%xlabel('w [m/day]')
%ylabel('Relative Biomass [%]')
grid on; axis tight; xlim([-150 150])
ylim([1e-2 100])
title('Winter','fontsize',50,'position',[-90 11e+01 0]);
text(-143,6e1,'(b)   \xi = 2','fontweight','b','FontSize',28)


axes('position',[.8 .82 .15 .14])
bar(1,sum(N0025w(Ew<=0).*B0025_2.*Ew(Ew<=0))/totalFw*100,'facecolor',color0025,'edgecolor',color0025)
text(1,sum(N0025w(Ew<=0).*B0025_2.*Ew(Ew<=0))/totalFw*100,[num2str(sum(N0025w(Ew<=0).*B0025_2.*Ew(Ew<=0))/totalFw*100,'%2.1f'),'%'],'horizontalalignment','center','verticalalignment','top','fontsize',17,'fontweight','b')
hold on
bar(2,sum(N1w(Ew<=0).*B1_2.*Ew(Ew<=0))/totalFw*100,'facecolor',color1,'edgecolor',color1);
text(2,sum(N1w(Ew<=0).*B1_2.*Ew(Ew<=0))/totalFw*100,[num2str(sum(N1w(Ew<=0).*B1_2.*Ew(Ew<=0))/totalFw*100,'%2.1f'),'%'],'horizontalalignment','center','verticalalignment','top','fontsize',17,'fontweight','b')
bar(3,sum(N5w(Ew<=0).*B5_2.*Ew(Ew<=0))/totalFw*100,'facecolor',color5,'edgecolor',color5);
text(3,sum(N5w(Ew<=0).*B5_2.*Ew(Ew<=0))/totalFw*100,[num2str(sum(N5w(Ew<=0).*B5_2.*Ew(Ew<=0))/totalFw*100,'%2.1f'),'%'],'horizontalalignment','center','verticalalignment','top','fontsize',17,'fontweight','b','color','w')
set(gca,'YScale','log','FontSize',28,'xtick',[1:3],'XTickLabel',{},'yaxislocation','right')
grid on; box on
ylim([1e-2 100])
title('Relative Downward Flux','FontSize',28,'position',[2 200 0])

subplot(2,2,4)
h1 = plot(Ew,N0025w.*B0025_4/totalB4w*100,'-','color',color0025,'linewidth',2);
hold on
plot(Ew,N1w.*B1_4/totalB4w*100,'-','color',color1,'linewidth',2)
plot(Ew,N5w.*B5_4/totalB4w*100,'-','color',color5,'linewidth',2)

% Shade negative flux
patch([Ew(Ew<=0) 0 min(Ew(Ew<=0)) min(Ew(Ew<=0))],[N0025w(Ew<=0).*B0025_4/totalB4w*100 1e-2 1e-2 N0025w(1).*B0025_4/totalB4w*100],color0025,'edgecolor','none','facealpha',.3)
patch([Ew(Ew<=0) 0 min(Ew(Ew<=0)) min(Ew(Ew<=0))],[N1w(Ew<=0).*B1_4/totalB4w*100 1e-2 1e-2 N1w(1).*B1_4/totalB4w*100],color1,'edgecolor','none','facealpha',.3)
patch([Ew(Ew<=0) 0 min(Ew(Ew<=0)) min(Ew(Ew<=0))],[N5w(Ew<=0).*B5_4/totalB4w*100 1e-2 1e-2 N5w(1).*B5_4/totalB4w*100],color5,'edgecolor','none','facealpha',.3)

set(gca,'Yscale','log','FontSize',28,'yticklabel',{})
xlabel('w_t_o_t [m/day]')
%ylabel('Relative Biomass [%]')
grid on; axis tight; xlim([-150 150])
ylim([1e-2 100])
text(-143,6e1,'(d)   \xi = 4','fontweight','b','FontSize',28)


axes('position',[.8 .34 .15 .14])
bar(1,sum(N0025w(Ew<=0).*B0025_4.*Ew(Ew<=0))/totalF4w*100,'facecolor',color0025,'edgecolor',color0025)
text(1,sum(N0025w(Ew<=0).*B0025_4.*Ew(Ew<=0))/totalF4w*100,[num2str(sum(N0025w(Ew<=0).*B0025_4.*Ew(Ew<=0))/totalF4w*100,'%2.1f'),'%'],'horizontalalignment','center','verticalalignment','top','fontsize',17,'fontweight','b')
hold on
bar(2,sum(N1w(Ew<=0).*B1_4.*Ew(Ew<=0))/totalF4w*100,'facecolor',color1,'edgecolor',color1);
text(2,sum(N1w(Ew<=0).*B1_4.*Ew(Ew<=0))/totalF4w*100,[num2str(sum(N1w(Ew<=0).*B1_4.*Ew(Ew<=0))/totalF4w*100,'%2.1f'),'%'],'horizontalalignment','center','verticalalignment','top','fontsize',17,'fontweight','b')
bar(3,sum(N5w(Ew<=0).*B5_4.*Ew(Ew<=0))/totalF4w*100,'facecolor',color5,'edgecolor',color5);
text(3,sum(N5w(Ew<=0).*B5_4.*Ew(Ew<=0))/totalF4w*100,[num2str(sum(N5w(Ew<=0).*B5_4.*Ew(Ew<=0))/totalF4w*100,'%2.1f'),'%'],'horizontalalignment','center','verticalalignment','top','fontsize',17,'fontweight','b','color','w')
set(gca,'YScale','log','FontSize',28,'xtick',[1:3],'XTickLabel',{},'yaxislocation','right')
grid on; box on
ylim([1e-2 100])
title('Relative Downward Flux','FontSize',28,'position',[2 200 0])

set(gcf,'color','w')
export_fig -r300 Fig7_revised.png
close
%%


% Calculate total biomass to normalize by
totalB = sum(TS.N0025s.*B0025_2)+sum(TS.N1s.*B1_2)+sum(TS.N5s.*B5_2);
totalBw = sum(TS.N0025w.*B0025_2)+sum(TS.N1w.*B1_2)+sum(TS.N5w.*B5_2);
totalB4 = sum(TS.N0025s.*B0025_4)+sum(TS.N1s.*B1_4)+sum(TS.N5s.*B5_4);
totalB4w = sum(TS.N0025w.*B0025_4)+sum(TS.N1w.*B1_4)+sum(TS.N5w.*B5_4);

totalF = sum(TS.N0025s(E<=0,:).*B0025_2.*E(E<=0)')+sum(TS.N1s(E<=0,:).*B1_2.*E(E<=0)')+sum(TS.N5s(E<=0,:).*B5_2.*E(E<=0)');
totalF4 = sum(TS.N0025s(E<=0,:).*B0025_4.*E(E<=0)')+sum(TS.N1s(E<=0,:).*B1_4.*E(E<=0)')+sum(TS.N5s(E<=0,:).*B5_4.*E(E<=0)');

totalFw = sum(TS.N0025w(Ew<=0,:).*B0025_2.*Ew(Ew<=0)')+sum(TS.N1w(Ew<=0,:).*B1_2.*Ew(Ew<=0)')+sum(TS.N5w(Ew<=0,:).*B5_2.*Ew(Ew<=0)');
totalF4w = sum(TS.N0025w(Ew<=0,:).*B0025_4.*Ew(Ew<=0)')+sum(TS.N1w(Ew<=0,:).*B1_4.*Ew(Ew<=0)')+sum(TS.N5w(Ew<=0,:).*B5_4.*Ew(Ew<=0)');


figure('position',[15 150 2119 1361])
subplot(2,2,1)
h1 = plot(0:.125:24,sum(TS.N0025s(E<=0,:).*B0025_2.*E(E<=0)')./totalF*100,'color',color0025,'linewidth',2);
hold on
h2 = plot(0:.125:24,sum(TS.N1s(E<=0,:).*B1_2.*E(E<=0)')./totalF*100,'color',color1,'linewidth',2);
h3 = plot(0:.125:24,sum(TS.N5s(E<=0,:).*B5_2.*E(E<=0)')./totalF*100,'color',color5,'linewidth',2);
set(gca,'fontsize',18)
ylabel('Relative Biomass [%]')
grid on; axis tight;
ylim([0 100])
title('Summer','fontsize',30);
text(.5,95,'(a)   \xi = 2','fontweight','b','fontsize',20)
HH = legend([h1 h2 h3],'0.025 m/day','1 m/day','5 m/day','location','best');
set(HH,'fontsize',18,'linewidth',2)


subplot(2,2,2)
h1 = plot(0:.125:24,sum(TS.N0025w(Ew<=0,:).*B0025_2.*Ew(Ew<=0)')./totalFw*100,'color',color0025,'linewidth',2);
hold on
h2 = plot(0:.125:24,sum(TS.N1w(Ew<=0,:).*B1_2.*Ew(Ew<=0)')./totalFw*100,'color',color1,'linewidth',2);
h3 = plot(0:.125:24,sum(TS.N5w(Ew<=0,:).*B5_2.*Ew(Ew<=0)')./totalFw*100,'color',color5,'linewidth',2);
set(gca,'fontsize',18)
ylabel('Relative Biomass [%]')
grid on; axis tight;
ylim([0 100])
title('Winter','fontsize',30);
text(.5,95,'(b)   \xi = 2','fontweight','b','fontsize',20)

subplot(2,2,3)
h1 = plot(0:.125:24,sum(TS.N0025s(E<=0,:).*B0025_4.*E(E<=0)')./totalF4*100,'color',color0025,'linewidth',2);
hold on
h2 = plot(0:.125:24,sum(TS.N1s(E<=0,:).*B1_4.*E(E<=0)')./totalF4*100,'color',color1,'linewidth',2);
h3 = plot(0:.125:24,sum(TS.N5s(E<=0,:).*B5_4.*E(E<=0)')./totalF4*100,'color',color5,'linewidth',2);
set(gca,'fontsize',18)
ylabel('Relative Biomass [%]')
xlabel('Day')
grid on; axis tight;
ylim([0 100])
text(.5,95,'(c)   \xi = 4','fontweight','b','fontsize',20)


subplot(2,2,4)
h1 = plot(0:.125:24,sum(TS.N0025w(Ew<=0,:).*B0025_4.*Ew(Ew<=0)')./totalF4w*100,'color',color0025,'linewidth',2);
hold on
h2 = plot(0:.125:24,sum(TS.N1w(Ew<=0,:).*B1_4.*Ew(Ew<=0)')./totalF4w*100,'color',color1,'linewidth',2);
h3 = plot(0:.125:24,sum(TS.N5w(Ew<=0,:).*B5_4.*Ew(Ew<=0)')./totalF4w*100,'color',color5,'linewidth',2);
set(gca,'fontsize',18)
ylabel('Relative Biomass [%]')
xlabel('Day')
grid on; axis tight;
ylim([0 100])
title('Winter','fontsize',30);
text(.5,95,'(d)   \xi = 4','fontweight','b','fontsize',20)



set(gcf,'color','w')
export_fig -r300 Fig7b_revised.png
close