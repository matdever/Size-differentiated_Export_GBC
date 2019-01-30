clear

% this code import the number of simulated particles located below 100 m in
% Papa_meso and Papa_submeso

%% Start with Papa_meso

path = '/Volumes/garuda/Mathieu/particle_tracking/MATLAB/Papa500_';
for ii = 1:5
    if ii ==1
        filename = [path,'0mday/exportTS0mday.txt'];
    elseif ii == 2
        filename = [path,'0025mday/exportTS0025mday.txt'];
    elseif ii == 3
        filename = [path,'05mday/exportTS05mday.txt'];
    elseif ii == 4
        filename = [path,'1mday/exportTS1mday.txt'];
    elseif ii == 5
        filename = [path,'5mday/exportTS5mday.txt'];
    end
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    
    export_meso(:,ii) = dataArray{:, 1};
    
    clearvars filename delimiter formatSpec fileID dataArray ans
end; clear ii path

%% Compute export anomaly based on expected export computed from sinking
% velocity and distance from 100m horizon

N0 = 1971717;
t0 = 60.5;
%ws = 5;

time_meso = 60.5:0.125:88.375;
Nmeso = zeros(size(export_meso));

wscounter = 1;
for ws = [0 0.025 0.5 1 5]
    counter = 1;
    for zz = 15:25
        if zz == 15
            Nmeso(time_meso<t0+zz/ws,wscounter) = 0;
        end
        if zz == 25
            Nmeso(time_meso>t0+zz/ws,wscounter) = 1;
        end
        Nmeso(time_meso>=t0+zz/ws & time_meso<t0+(zz+1)/ws,wscounter) = counter/11;
        counter = counter+1;
    end
    wscounter = wscounter +1;
end; clear ws t0 counter wscounter zz


%
%
% figure(2)
% plot(time,(data./Ntot-N)*100,'--','LineWidth',2)
% legend('w_s = 0','w_s = 0.025','w_s = 0.5','w_s = 1','w_s = 5')
% title('Export anomaly')

%% Then do Papa_submeso

path = '/Volumes/garuda/Mathieu/particle_tracking/MATLAB/Papa500steroid_';
for ii = 1:5
    if ii ==1
        filename = [path,'0mday/exportTS0mday.txt'];
    elseif ii == 2
        filename = [path,'0025mday/exportTS0025mday.txt'];
    elseif ii == 3
        filename = [path,'05mday/exportTS05mday.txt'];
    elseif ii == 4
        filename = [path,'1mday/exportTS1mday.txt'];
    elseif ii == 5
        filename = [path,'5mday/exportTS5mday.txt'];
    end
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    
    export_submeso(:,ii) = dataArray{:, 1};
    
    clearvars filename delimiter formatSpec fileID dataArray ans
end; clear ii path

    export_submeso(229:end,:) = [];


%% Compute export anomaly based on expected export computed from sinking
% velocity and distance from 100m horizon

N0 = 1971717;
t0 = 50;

time_submeso = 50:0.125:80;
Nsubmeso = zeros(size(export_submeso));

wscounter = 1;
for ws = [0 0.025 0.5 1 5]
    counter = 1;
    for zz = 15:25
        if zz == 15
            Nsubmeso(time_submeso<t0+zz/ws,wscounter) = 0;
        end
        if zz == 25
            Nsubmeso(time_submeso>t0+zz/ws,wscounter) = 1;
        end
        Nsubmeso(time_submeso>=t0+zz/ws & time_submeso<t0+(zz+1)/ws,wscounter) = counter/11;
        counter = counter+1;
    end
    wscounter = wscounter +1;
end; clear ws t0 counter wscounter zz

    Nsubmeso(229:end,:) = [];
    time_submeso(229:end) = [];

%% figure

cmap = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

% non-caucher overwriting
Nmeso(120:201,4) = NaN;
Nmeso(24:41,5) = NaN;
Nsubmeso(120:201,4) = NaN;
Nsubmeso(24:41,5) = NaN;

fullfigure
subplot(3,2,1)
for ii = 5:-1:2
h(5-ii+1) = plot(time_meso+105,export_meso(:,ii)/N0*100,'linewidth',2,'color',cmap(ii,:));
hold on
plot(time_meso(isnan(Nmeso(:,ii))==0)+105,Nmeso(isnan(Nmeso(:,ii))==0,ii)*100,'--','color',[.3 .3 .3],'linewidth',3,'color',cmap(ii,:));
end; clear ii
XLIM = get(gca,'xlim');
g = plot(time_meso,Nmeso(:,1)*100,'--','color',[.3 .3 .3],'linewidth',2);
set(gca,'xlim',XLIM); 
H = legend([h g],'w_s = 5  ','w_s = 1  ','w_s = 0.5  ','w_s = 0.025  ','expected','location','northoutside','Orientation','horizontal');
set(H,'fontsize',28);clear H
clear g h XLIM
text(180,150,'Papa\_summer','HorizontalAlignment','center','FontWeight','b','FontSize',30)
grid on
ylabel('Exported Particles [%]')
xlabel('Day of Year')
set(gca,'fontsize',28)

subplot(3,2,2)
for ii = 5:-1:2
h(5-ii+1) = plot(time_submeso,export_submeso(:,ii)/N0*100,'linewidth',2,'color',cmap(ii,:));
hold on
plot(time_submeso(isnan(Nsubmeso(:,ii))==0),Nsubmeso((isnan(Nsubmeso(:,ii))==0),ii)*100,'--','color',[.3 .3 .3],'linewidth',3,'color',cmap(ii,:));
end; clear ii
XLIM = get(gca,'xlim');
g = plot(time_submeso,Nsubmeso(:,1)*100,'--','color',[.3 .3 .3],'linewidth',2);
set(gca,'xlim',XLIM); 
H = legend([h g],'w_s = 5  ','w_s = 1  ','w_s = 0.5  ','w_s = 0.025  ','expected','location','northoutside','Orientation','horizontal');
set(H,'fontsize',28);clear H
clear g h XLIM
text(65,150,'Papa\_winter','HorizontalAlignment','center','FontWeight','b','FontSize',30)
grid on
ylabel('Exported Particles [%]')
xlabel('Day of Year')
set(gca,'fontsize',28)

% 
% subplot(3,2,3)
% for ii = 5:-1:1
% h(5-ii+1) = plot(time_meso+105,smooth(export_meso(:,ii)/N0*100-Nmeso(:,ii)*100,20),'linewidth',2,'color',cmap(ii,:));
% hold on
% end; clear ii
% grid on
% ylabel('Exported Particles Anomaly [%]')
% xlabel('Day of Year')
% set(gca,'fontsize',21)
% H = legend(h,'w_s = 5','w_s = 1','w_s = 0.5','w_s = 0.025','w_s = 0','location','northoutside','Orientation','horizontal');
% set(H,'fontsize',21)
% clear h H
% 
% subplot(3,2,4)
% for ii = 5:-1:1
% h(5-ii+1) = plot(time_submeso,smooth(export_submeso(:,ii)/N0*100-Nsubmeso(:,ii)*100,10),'linewidth',2,'color',cmap(ii,:));
% hold on
% end; clear ii
% grid on
% ylabel('Exported Particles Anomaly [%]')
% xlabel('Day of Year')
% set(gca,'fontsize',21)
% H = legend(h,'w_s = 5','w_s = 1','w_s = 0.5','w_s = 0.025','w_s = 0','location','northoutside','Orientation','horizontal');
% set(H,'fontsize',21)
% clear h H

set(gcf,'color','w')
export_fig -r200 Fig6_particle_export.png
close all