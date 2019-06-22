clear

%% This figure plots the histogram of the biomass vs vertical velocity at 
% DOY 25 in Summer/Winter with and without remineralization for slope = 4

%% Import w for each size class in winter and summer
tt = 25;
[~,~,wt0025s] = importw('/Volumes/mdever/particlespeed/doy25/w0025mday_summer.csv');
[~,~,wt1s] = importw('/Volumes/mdever/particlespeed/doy25/w1mday_summer.csv');
[~,~,wt5s] = importw('/Volumes/mdever/particlespeed/doy25/w5mday_summer.csv');

[~,~,wt0025w] = importw('/Volumes/mdever/particlespeed/doy25/w0025mday_winter.csv');
[~,~,wt1w] = importw('/Volumes/mdever/particlespeed/doy25/w1mday_winter.csv');
[~,~,wt5w] = importw('/Volumes/mdever/particlespeed/doy25/w5mday_winter.csv');

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

% For Junge slope = 4
slope = 4;
% Compute biomass of 1 largest particle necessary to keep total biomass constant
B0 = Btot*(wmax.^((3-slope)/2))/(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));
B0025_4 = B0*(0.025/wmax)^((3-slope)/2);
B1_4 = B0*(1/wmax)^((3-slope)/2);
B5_4 = B0*(5/wmax)^((3-slope)/2);

% Compute histograms in summer
[N0025s,E] = histcounts(wt0025s*86400,-10.05:.1:10.05);
[N1s,~] = histcounts(wt1s*86400,-10.05:.1:10.05);
[N5s,~] = histcounts(wt5s*86400,-10.05:.1:10.05);
E = (E(1:end-1)+E(2:end))/2;

N0025s(N0025s==0) = 1e-2; % to lok better in log plots
N1s(N1s==0) = 1e-2;
N5s(N5s==0) = 1e-2;

% Compute histograms in winter
[N0025w,Ew] = histcounts(wt0025w*86400,-155:10:155);
[N1w,~] = histcounts(wt1w*86400,-155:10:155);
[N5w,~] = histcounts(wt5w*86400,-155:10:155);
Ew = (Ew(1:end-1)+Ew(2:end))/2;


N0025w(N0025w==0) = 1e-2;
N1w(N1w==0) = 1e-2;
N5w(N5w==0) = 1e-2;

%% Plot
% Calculate total biomass to normalize by
totalB = sum(N0025s.*B0025_4)+sum(N1s.*B1_4)+sum(N5s.*B5_4);
totalBw = sum(N0025w.*B0025_4)+sum(N1w.*B1_4)+sum(N5w.*B5_4);
totalF = sum(N0025s(E<=0).*B0025_4.*E(E<=0))+sum(N1s(E<=0).*B1_4.*E(E<=0))+sum(N5s(E<=0).*B5_4.*E(E<=0));
totalFw = sum(N0025w(Ew<=0).*B0025_4.*Ew(Ew<=0))+sum(N1w(Ew<=0).*B1_4.*Ew(Ew<=0))+sum(N5w(Ew<=0).*B5_4.*Ew(Ew<=0));

figure('position',[15 150 2119 1361])
subplot(2,2,1)
h1 = plot(E,N0025s.*B0025_4/totalB*100,'-','color',color0025,'linewidth',2);
hold on
h2 = plot(E,N1s.*B1_4/totalB*100,'-','color',color1,'linewidth',2);
h3 = plot(E,N5s.*B5_4/totalB*100,'-','color',color5,'linewidth',2);
% Shade negative flux
patch([E(E<=0) 0 min(E(E<=0)) min(E(E<=0))],[N0025s(E<=0).*B0025_4/totalB*100 1e-2 1e-2 N0025s(1).*B0025_4/totalB*100],color0025,'edgecolor','none','facealpha',.3)
patch([E(E<=0) 0 min(E(E<=0)) min(E(E<=0))],[N1s(E<=0).*B1_4/totalB*100 1e-2 1e-2 N1s(1).*B1_4/totalB*100],color1,'edgecolor','none','facealpha',.3)
patch([E(E<=0) 0 min(E(E<=0)) min(E(E<=0))],[N5s(E<=0).*B5_4/totalB*100 1e-2 1e-2 N5s(1).*B5_4/totalB*100],color5,'edgecolor','none','facealpha',.3)

set(gca,'Yscale','log','FontSize',40)
%xlabel('w [m/day]')
ylabel('Relative Biomass [%]')
grid on; axis tight; xlim([-10 10])
ylim([1e-2 100])
title('Summer','fontsize',50,'position',[-6 11e+01 0]);
text(-9.5,6e1,'(a)   Remin OFF','fontweight','b','fontsize',40)
HH = legend([h1 h2 h3],'0.025 m/day','1 m/day','5 m/day','location','southeast');
set(HH,'fontsize',25)

axes('position',[.36 .82 .15 .14])
bar(1,sum(N0025s(E<=0).*B0025_4.*E(E<=0))/totalF*100,'FaceColor',color0025,'EdgeColor',color0025)
hold on
bar(2,sum(N1s(E<=0).*B1_4.*E(E<=0))/totalF*100,'facecolor',color1,'edgecolor',color1);
bar(3,sum(N5s(E<=0).*B5_4.*E(E<=0))/totalF*100,'facecolor',color5,'edgecolor',color5);
set(gca,'YScale','log','fontsize',35,'xtick',[1:3],'XTickLabel',{},'yaxislocation','right')
grid on; box on
ylim([1e-2 100])
title('Relative Downward Flux','fontsize',40,'position',[2 200 0])


subplot(2,2,2)
h1 = plot(Ew,N0025w.*B0025_4/totalBw*100,'-','color',color0025,'linewidth',2);
hold on
plot(Ew,N1w.*B1_4/totalBw*100,'-','color',color1,'linewidth',2)
plot(Ew,N5w.*B5_4/totalBw*100,'-','color',color5,'linewidth',2)
% Shade negative flux
patch([Ew(Ew<=0) 0 min(Ew(Ew<=0)) min(Ew(Ew<=0))],[N0025w(Ew<=0).*B0025_4/totalBw*100 1e-2 1e-2 N0025w(1).*B0025_4/totalBw*100],color0025,'edgecolor','none','facealpha',.3)
patch([Ew(Ew<=0) 0 min(Ew(Ew<=0)) min(Ew(Ew<=0))],[N1w(Ew<=0).*B1_4/totalBw*100 1e-2 1e-2 N1w(1).*B1_4/totalBw*100],color1,'edgecolor','none','facealpha',.3)
patch([Ew(Ew<=0) 0 min(Ew(Ew<=0)) min(Ew(Ew<=0))],[N5w(Ew<=0).*B5_4/totalBw*100 1e-2 1e-2 N5w(1).*B5_4/totalBw*100],color5,'edgecolor','none','facealpha',.3)

set(gca,'Yscale','log','FontSize',40,'yticklabel',{})
%xlabel('w [m/day]')
%ylabel('Relative Biomass [%]')
grid on; axis tight; xlim([-150 150])
ylim([1e-2 100])
title('Winter','fontsize',50,'position',[-90 11e+01 0]);
text(-143,6e1,'(b)   Remin OFF','fontweight','b','fontsize',40)

axes('position',[.8 .82 .15 .14])
bar(1,sum(N0025w(Ew<=0).*B0025_4.*Ew(Ew<=0))/totalFw*100,'FaceColor',color0025,'EdgeColor',color0025)
hold on
bar(2,sum(N1w(Ew<=0).*B1_4.*Ew(Ew<=0))/totalFw*100,'facecolor',color1,'edgecolor',color1);
bar(3,sum(N5w(Ew<=0).*B5_4.*Ew(Ew<=0))/totalFw*100,'facecolor',color5,'edgecolor',color5);
set(gca,'YScale','log','fontsize',35,'xtick',[1:3],'XTickLabel',{},'yaxislocation','right')
grid on; box on
ylim([1e-2 100])
title('Relative Downward Flux','fontsize',40,'position',[2 200 0])

%%
clearvars -except total*
tt = 25;
[~,~,wt0025s] = importw('/Volumes/mdever/particlespeed/doy25/w0025mday_remin_summer.csv');
[~,~,wt1s] = importw('/Volumes/mdever/particlespeed/doy25/w1mday_remin_summer.csv');
[~,~,wt5s] = importw('/Volumes/mdever/particlespeed/doy25/w5mday_remin_summer.csv');

[~,~,wt0025w] = importw('/Volumes/mdever/particlespeed/doy25/w0025mday_remin_winter.csv');
[~,~,wt1w] = importw('/Volumes/mdever/particlespeed/doy25/w1mday_remin_winter.csv');
[~,~,wt5w] = importw('/Volumes/mdever/particlespeed/doy25/w5mday_remin_winter.csv');

%% Set colors for each class
color0025 = [215,48,39]/256;
color1 = [39,100,25]/256;
color5 = [49,54,149]/256;

%% compute biomass per simulated particle
wmax = 5; wmin = 0.025;
slope = 3; B0 = 1;
Btot = B0/(wmax.^((3-slope)/2))*(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));

slope = 4;
B0 = Btot*(wmax.^((3-slope)/2))/(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));
B0025_4 = B0*(0.025/wmax)^((3-slope)/2)*exp(-0.13*tt);
B1_4 = B0*(1/wmax)^((3-slope)/2)*exp(-0.13*tt);
B5_4 = B0*(5/wmax)^((3-slope)/2)*exp(-0.13*tt);

% Compute histograms
[N0025s,E] = histcounts(wt0025s*86400,-10.05:.1:10.05);
[N1s,~] = histcounts(wt1s*86400,-10.05:.1:10.05);
[N5s,~] = histcounts(wt5s*86400,-10.05:.1:10.05);
E = (E(1:end-1)+E(2:end))/2;

[N0025w,Ew] = histcounts(wt0025w*86400,-155:10:155);
[N1w,~] = histcounts(wt1w*86400,-155:10:155);
[N5w,~] = histcounts(wt5w*86400,-155:10:155);
Ew = (Ew(1:end-1)+Ew(2:end))/2;

N0025s(N0025s==0) = 1e-2;
N1s(N1s==0) = 1e-2;
N5s(N5s==0) = 1e-2;

N0025w(N0025w==0) = 1e-2;
N1w(N1w==0) = 1e-2;
N5w(N5w==0) = 1e-2;
%%
% Calculate total biomass to normalize by
% totalB = sum(N0025s.*B0025_4)+sum(N1s.*B1_4)+sum(N5s.*B5_4);
% totalBw = sum(N0025w.*B0025_4)+sum(N1w.*B1_4)+sum(N5w.*B5_4);
% totalF = sum(N0025s(E<=0).*B0025_4.*E(E<=0))+sum(N1s(E<=0).*B1_4.*E(E<=0))+sum(N5s(E<=0).*B5_4.*E(E<=0));
% totalFw = sum(N0025w(Ew<=0).*B0025_4.*Ew(Ew<=0))+sum(N1w(Ew<=0).*B1_4.*Ew(Ew<=0))+sum(N5w(Ew<=0).*B5_4.*Ew(Ew<=0));
% 

subplot(2,2,3)
h1 = plot(E,N0025s.*B0025_4/totalB*100,'-','color',color0025,'linewidth',2);
hold on
plot(E,N1s.*B1_4/totalB*100,'-','color',color1,'linewidth',2)
plot(E,N5s.*B5_4/totalB*100,'-','color',color5,'linewidth',2)
% Shade negative flux
patch([E(E<=0) 0 min(E(E<=0)) min(E(E<=0))],[N0025s(E<=0).*B0025_4/totalB*100 1e-2 1e-2 N0025s(1).*B0025_4/totalB*100],color0025,'edgecolor','none','facealpha',.3)
patch([E(E<=0) 0 min(E(E<=0)) min(E(E<=0))],[N1s(E<=0).*B1_4/totalB*100 1e-2 1e-2 N1s(1).*B1_4/totalB*100],color1,'edgecolor','none','facealpha',.3)
patch([E(E<=0) 0 min(E(E<=0)) min(E(E<=0))],[N5s(E<=0).*B5_4/totalB*100 1e-2 1e-2 N5s(1).*B5_4/totalB*100],color5,'edgecolor','none','facealpha',.3)

set(gca,'Yscale','log','FontSize',40)
xlabel('w [m/day]')
ylabel('Relative Biomass [%]')
grid on; axis tight; xlim([-10 10])
ylim([1e-2 100])
text(-9.5,6e1,'(c)   Remin ON','fontweight','b','fontsize',40)

axes('position',[.36 .34 .15 .14])
bar(1,sum(N0025s(E<=0).*B0025_4.*E(E<=0))/totalF*100,'FaceColor',color0025,'EdgeColor',color0025)
hold on
bar(2,sum(N1s(E<=0).*B1_4.*E(E<=0))/totalF*100,'facecolor',color1,'edgecolor',color1);
bar(3,sum(N5s(E<=0).*B5_4.*E(E<=0))/totalF*100,'facecolor',color5,'edgecolor',color5);
set(gca,'YScale','log','fontsize',35,'xtick',[1:3],'XTickLabel',{},'yaxislocation','right')
grid on; box on
ylim([1e-2 100])
title('Relative Downward Flux','fontsize',40,'position',[2 400 0])


subplot(2,2,4)
h1 = plot(Ew,N0025w.*B0025_4/totalBw*100,'-','color',color0025,'linewidth',2);
hold on
plot(Ew,N1w.*B1_4/totalBw*100,'-','color',color1,'linewidth',2)
plot(Ew,N5w.*B5_4/totalBw*100,'-','color',color5,'linewidth',2)
% Shade negative flux
patch([Ew(Ew<=0) 0 min(Ew(Ew<=0)) min(Ew(Ew<=0))],[N0025w(Ew<=0).*B0025_4/totalBw*100 1e-2 1e-2 N0025w(1).*B0025_4/totalBw*100],color0025,'edgecolor','none','facealpha',.3)
patch([Ew(Ew<=0) 0 min(Ew(Ew<=0)) min(Ew(Ew<=0))],[N1w(Ew<=0).*B1_4/totalBw*100 1e-2 1e-2 N1w(1).*B1_4/totalBw*100],color1,'edgecolor','none','facealpha',.3)
patch([Ew(Ew<=0) 0 min(Ew(Ew<=0)) min(Ew(Ew<=0))],[N5w(Ew<=0).*B5_4/totalBw*100 1e-2 1e-2 N5w(1).*B5_4/totalBw*100],color5,'edgecolor','none','facealpha',.3)

set(gca,'Yscale','log','FontSize',40,'yticklabel',{})
xlabel('w [m/day]')
%ylabel('Relative Biomass [%]')
grid on; axis tight; xlim([-150 150])
ylim([1e-2 100])
text(-143,6e1,'(d)   Remin ON','fontweight','b','fontsize',40)

axes('position',[.8 .34 .15 .14])
bar(1,sum(N0025w(Ew<=0).*B0025_4.*Ew(Ew<=0))/totalFw*100,'FaceColor',color0025,'EdgeColor',color0025)
hold on
bar(2,sum(N1w(Ew<=0).*B1_4.*Ew(Ew<=0))/totalFw*100,'facecolor',color1,'edgecolor',color1);
bar(3,sum(N5w(Ew<=0).*B5_4.*Ew(Ew<=0))/totalFw*100,'facecolor',color5,'edgecolor',color5);
set(gca,'YScale','log','fontsize',35,'xtick',[1:3],'XTickLabel',{},'yaxislocation','right')
grid on; box on
ylim([1e-2 100])
title('Relative Downward Flux','fontsize',40,'position',[2 400 0])

set(gcf,'color','w')
export_fig -r300 Fig8.png
%%
function [w,wsink,wtotal] = importw(filename)

%% Initialize variables.
delimiter = ',';
startRow = 1;
endRow = inf;

%% Format for each line of text:
formatSpec = '%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');


%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
w = dataArray{:, 1};
wsink = dataArray{:, 2};
wtotal = w+wsink;
end

