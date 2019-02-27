clear

% this code import the number of simulated particles located below 100 m in
% Papa_meso and Papa_submeso

%% Start with Papa_summer

% total number of particles
N0 = 1971717;
% start day
t0 = 60.5;
% Time vector.
time_meso = 60.5:0.125:88.375;

path = '/Volumes/mahadevanlab/mathieudever/particle_backup/Papa_summer/';
for ii = 1:5
    if ii ==1
        continue
    elseif ii == 2
        continue
    elseif ii == 3
        filenamedown = [path,'05mday/exportTS05mdaydown.txt'];
        filenameup = [path,'05mday/exportTS05mdayup.txt'];
    elseif ii == 4
        filenamedown = [path,'1mday/exportTS1mdaydown.txt'];
        filenameup = [path,'1mday/exportTS1mdayup.txt'];
    elseif ii == 5
        filenamedown = [path,'5mday/exportTS5mdaydown.txt'];
        filenameup = [path,'5mday/exportTS5mdayup.txt'];
    end
    
    % Get downward flux
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(filenamedown,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    flux_meso_down(:,ii) = dataArray{:, 1};
    clearvars filenamedown delimiter formatSpec fileID dataArray ans
    
    % Get upward flux
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(filenameup,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    flux_meso_up(:,ii) = dataArray{:, 1};
    
    clearvars filenameup delimiter formatSpec fileID dataArray ans
end; clear ii path

flux_meso_down(flux_meso_down==0) = 1e-10;
flux_meso_up(flux_meso_up==0) = 1e-10;
%% Compute export anomaly based on expected export computed from sinking
% velocity and distance from 100m horizon

grav_flux_meso = zeros(size(flux_meso_up));
wscounter = 2;
for ws = [0.025 0.5 1 5]
    for zz = 15:25
        ind = find(t0+zz/ws<=time_meso,1,'first');
        grav_flux_meso(ind,wscounter) = N0./11./0.125;
    end; clear zz ind
    grav_flux_meso(:,wscounter) = smooth(grav_flux_meso(:,wscounter),8);
    a =find(grav_flux_meso(:,wscounter)~=0,1,'first');
    b = find(grav_flux_meso(:,wscounter)~=0,1,'last');
    grav_flux_meso(a:b,wscounter) = max(grav_flux_meso(a:b,wscounter));
    wscounter = wscounter +1;
end; clear ws wscounter a b

grav_flux_meso(grav_flux_meso==0) = 1e-10;

%% Then do Papa_submeso
N0 = 1971717;
t0 = 50;
time_submeso = 50:0.125:80;

path = '/Volumes/mahadevanlab/mathieudever/particle_backup/Papa_winter/';
for ii = 1:5
    if ii ==1
        continue
    elseif ii == 2
        filenamedown = [path,'0025mday/exportTS0025mdaydown.txt'];
        filenameup = [path,'0025mday/exportTS0025mdayup.txt'];
    elseif ii == 3
        filenamedown = [path,'05mday/exportTS05mdaydown.txt'];
        filenameup = [path,'05mday/exportTS05mdayup.txt'];
    elseif ii == 4
        filenamedown = [path,'1mday/exportTS1mdaydown.txt'];
        filenameup = [path,'1mday/exportTS1mdayup.txt'];
    elseif ii == 5
        filenamedown = [path,'5mday/exportTS5mdaydown.txt'];
        filenameup = [path,'5mday/exportTS5mdayup.txt'];
    end
    
    % Get downward flux
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(filenamedown,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    flux_submeso_down(:,ii) = dataArray{:, 1};
    clearvars filenamedown delimiter formatSpec fileID dataArray ans
    
    % Get upward flux
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(filenameup,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    flux_submeso_up(:,ii) = dataArray{:, 1};
    
    clearvars filenameup delimiter formatSpec fileID dataArray ans
end; clear ii path

flux_submeso_up(flux_submeso_up==0) = 1e-10;
flux_submeso_down(flux_submeso_down==0) = 1e-10;

%% Compute export anomaly based on expected export computed from sinking
% velocity and distance from 100m horizon

grav_flux_submeso = zeros(size(flux_submeso_up));
wscounter = 2;
for ws = [0.025 0.5 1 5]
    for zz = 15:25
        ind = find(t0+zz/ws<=time_submeso,1,'first');
        grav_flux_submeso(ind,wscounter) = N0./11./0.125;
    end; clear zz ind
    wscounter = wscounter +1;
end; clear ws wscounter

grav_flux_submeso(grav_flux_submeso==0) = 1e-10;
%% figure

cmap = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

%% MESOSCALES
fullfigure
subplot(3,2,1)
for ii = 1:5
    plot(time_meso+105,smooth(flux_meso_up(:,ii)./unique(diff(time_meso)),12),'color',cmap(ii,:),'linewidth',1);
    hold on; grid on
    theh(ii) = line([0 0],[0 0],'color',cmap(ii,:),'linewidth',2);
end; clear ii
set(gca,'yscale','log','xticklabel',{},'fontsize',16,'xlim',[min(time_meso)+105 max(time_meso)+105])
HH = ylabel('\langleparticle flux\rangle [particles/day]','fontsize',18,'Interpreter','tex');
POS = get(HH,'position');
POS(2) = 1;
set(HH,'position',POS); clear POS
%plot(time_meso,smooth(grav_flux_meso(:,ii),8),'color',cmap(ii,:),'linestyle','--')
set(gca,'ylim',[10 12e5],'ytick',[10,100,1000,10000,100000]);
text(60.75+105,5e5,'(a)','fontweight','b','fontsize',18)
text(87.5+105,5e5,'Upward flux','fontweight','b','fontsize',18,'HorizontalAlignment','right')
text(180,5e6,'Papa\_summer','HorizontalAlignment','center','FontWeight','b','FontSize',30)

YLIM = get(gca,'ylim');
POS = get(gca,'position');
subplot(3,2,3)
for ii = 1:5
    plot(time_meso+105,smooth(flux_meso_down(:,ii)./unique(diff(time_meso)),12),'color',cmap(ii,:),'linewidth',2)
    hold on; grid on; axis ij
    plot(time_meso+105,grav_flux_meso(:,ii),'color',cmap(ii,:),'linestyle','--','Linewidth',2)
end; clear ii
H1 = line([0 0],[0 0],'color','k','linewidth',1);
H2 = line([0 0],[0 0],'color','k','linewidth',2);
H3 = line([0 0],[0 0],'color','k','linewidth',2,'linestyle','--');
POS2 = get(gca,'position');
POS2(2) = POS(2)-.005-POS2(4);
set(gca,'ylim',[10 12e5],'ytick',[10,100,1000,10000,100000],'YTickLabel',{'','10^2','10^3','10^4','10^5'},'xlim',[min(time_meso)+105 max(time_meso)+105]);
text(87.5+105,5e5,'Downward flux','fontweight','b','fontsize',16,'HorizontalAlignment','right')
set(gca,'yscale','log','ylim',[10 12e5],'position',POS2,'fontsize',16)
xlabel('Day of year','fontsize',18)

% LEGENDS
subplot(3,2,1)
H = legend(fliplr(theh(2:5)),'w_s = 5  ','w_s = 1  ','w_s = 0.5  ','w_s = 0.025 ','location','north','Orientation','horizontal');
set(H,'fontsize',18,'EdgeColor','k','position',[0.2029 0.86 0.1891 0.0265]);clear H
set(gca,'position',POS)
subplot(3,2,1)
H = legend([H1 H2 H3],'Upward/Downward flux','Net flux','Gravitational flux','location','north','Orientation','horizontal');
set(H,'fontsize',18,'EdgeColor','k','position',[0.1862 0.83 0.2224 0.0208]);clear H


% SUBMESOSCALES
subplot(3,2,2)
for ii = 1:5
    plot(time_submeso,smooth(flux_submeso_up(:,ii)./unique(diff(time_submeso)),12),'color',cmap(ii,:),'linewidth',1)
    hold on; grid on
end; clear ii
set(gca,'yscale','log','xticklabel',{},'fontsize',16,'xlim',[min(time_submeso) max(time_submeso)])
%HH = ylabel('\langleparticle flux\rangle [particles/day]','fontsize',18,'Interpreter','tex');
%POS = get(HH,'position');
%POS(2) = 1;
%set(HH,'position',POS); clear POS
POS = get(gca,'Position');
POS(1) = 0.48;
set(gca,'position',POS);
set(gca,'ylim',[10 12e5],'ytick',[10,100,1000,10000,100000],'yticklabel',{});
text(50.25,5e5,'(b)','fontweight','b','fontsize',18)
text(79,5e5,'Upward flux','fontweight','b','fontsize',18,'HorizontalAlignment','right')
text(65,5e6,'Papa\_winter','HorizontalAlignment','center','FontWeight','b','FontSize',30)

YLIM = get(gca,'ylim');
POS = get(gca,'position');
subplot(3,2,4)
for ii = 1:5
    plot(time_submeso,smooth(flux_submeso_down(:,ii)./unique(diff(time_submeso)),12),'color',cmap(ii,:),'linewidth',1)
    hold on; grid on; axis ij
    blop = flux_submeso_down(:,ii)-flux_submeso_up(:,ii);
    plot(time_submeso,smooth(blop./unique(diff(time_submeso)),12),'color',cmap(ii,:),'linewidth',2)
end; clear ii
POS2 = get(gca,'position');
POS2(2) = POS(2)-.005-POS2(4);
POS2(1) = POS(1);
set(gca,'ylim',[10 12e5],'ytick',[10,100,1000,10000,100000],'YTickLabel',{''},'xlim',[min(time_submeso) max(time_submeso)]);
text(79,5e5,'Downward flux','fontweight','b','fontsize',16,'HorizontalAlignment','right')
set(gca,'yscale','log','ylim',[10 12e5],'position',POS2,'fontsize',16)
xlabel('Day of year','fontsize',18)

set(gcf,'color','w')
export_fig -r300 Fig6_particle_export_flux.png
close all