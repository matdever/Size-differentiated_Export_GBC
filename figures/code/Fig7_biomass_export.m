clear

% this code import the number of simulated particles located below 100 m in
% Papa_meso and Papa_submeso

%% Start with Papa_meso

time_meso = 60.5:0.125:88.375;

path = '/Volumes/mathieudever/particle_backup/Papa500/';
for ii = 3:5
    if ii ==1
        filename = [path,'0mday/exportTS0mday.txt'];
    elseif ii == 2
        filename = [path,'0025mday_remin/exportTS0025mday_remin.txt'];
    elseif ii == 3
        filename = [path,'05mday_remin/exportTS05mday_remin.txt'];
    elseif ii == 4
        filename = [path,'1mday_remin/exportTS1mday_remin.txt'];
    elseif ii == 5
        filename = [path,'5mday_remin/exportTS5mday_remin.txt'];
    end
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    
    export_meso_remin(:,ii) = dataArray{:, 1};
    
    clearvars filename delimiter formatSpec fileID dataArray ans
end; clear ii path

time_meso = 60.5:0.125:88.375;

path = '/Volumes/mathieudever/particle_backup/Papa500/';
for ii = 3:5
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

%% Then do Papa_submeso

time_submeso = 50:0.125:80;

path = '/Volumes/mathieudever/particle_backup/Papa500steroid/Papa500steroid_';
for ii = 2:5
    if ii ==1
        filename = [path,'0mday_remin/exportTS0mday_remin.txt'];
    elseif ii == 2
        filename = [path,'0025mday_remin/exportTS0025mday_remin.txt'];
    elseif ii == 3
        filename = [path,'05mday_remin/exportTS05mday_remin.txt'];
    elseif ii == 4
        filename = [path,'1mday_remin/exportTS1mday_remin.txt'];
    elseif ii == 5
        filename = [path,'5mday_remin/exportTS5mday_remin.txt'];
    end
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    
    export_submeso_remin(:,ii) = dataArray{:, 1};
    
    clearvars filename delimiter formatSpec fileID dataArray ans
end; clear ii path


time_submeso = 50:0.125:80;

path = '/Volumes/mathieudever/particle_backup/Papa500steroid/Papa500steroid_';
for ii = 2:5
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
%% Figure
ws = [0 0.025 0.5 1 5];
ws0 = max(ws);
%ws0_meso = max(ws)/86400 *... % initial sinking velocity
%    exp(-2/3*0.13*(time_meso-time_meso(1)));
%ws0_submeso = max(ws)/86400 *... % initial sinking velocity
%    exp(-2/3*0.13*(time_submeso-time_submeso(1)));
B0 = 1; %arbitrary biomass unit

for ii = 1:length(ws)
    ws_meso(ii,:) = ws(ii)./86400 .*... % initial sinking velocity
        exp(-2/3*0.13*(time_meso-time_meso(1)));
end; clear ii

for ii = 1:length(ws)
    ws_submeso(ii,:) = ws(ii)./86400 .*... % initial sinking velocity
        exp(-2/3*0.13*(time_submeso-time_submeso(1)));
end; clear ii

cmap = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

%% figure with log scales

fullfigure

counter = 0;
for slope = [4 3 2]
    
    for ii = 1:length(ws)
        B_meso(:,ii) = export_meso(:,ii).*B0.*(ws(ii)/ws0)^((3-slope)/2);
        B_meso_remin(:,ii) = export_meso_remin(:,ii).*... % Number of model particles exported
            B0*(ws(ii)/ws0).^(3/2).*... % biomass of 1 model particle at t = 0
            (ws_meso(ii,:)'./ws_meso(ii,1)).^(3/2).*... % Change in biomass vs time for 1 model particle
            (ws(ii)/ws0).^(-slope/2); % Convert to real particles using the spectral slope
        
        
        B_submeso(:,ii) = export_submeso(:,ii).*B0.*(ws(ii)/ws0)^((3-slope)/2);
        B_submeso_remin(:,ii) = export_submeso_remin(:,ii).*... % Number of model particles exported
            B0*(ws(ii)/ws0).^(3/2).*... % biomass of 1 model particle at t = 0
            (ws_submeso(ii,:)'./ws_submeso(ii,1)).^(3/2).*... % Change in biomass vs time for 1 model particle
            (ws(ii)/ws0).^(-slope/2); % Convert to real particles using the spectral slope
        
    end; clear ii
    
    subplot(3,2,counter*2+1)
    for ii = 5:-1:2
        h1(5-ii+1) = plot(time_meso+105,B_meso(:,ii),'LineWidth',2,'Color',cmap(ii,:),'linestyle','--');
        hold on
        h2(5-ii+1) = plot(time_meso+105,B_meso_remin(:,ii),'LineWidth',2,'Color',cmap(ii,:));
        h3(5-ii+1) = plot(time_meso+105,B_meso_remin(:,ii),'LineWidth',2,'Color',cmap(ii,:));
    end; clear ii
    H = legend(h2,'w_s = 5','w_s = 1','w_s = 0.5','w_s = 0.025','w_s = 0','location','northoutside','Orientation','horizontal');
    set(H,'fontsize',28);clear H
    set(gca,'yscale','log')
    grid on
    xlabel('Day of Year')
    ylabel({'Exported Biomass','[normalized units]'})
    XLIM = get(gca,'xlim');
    YLIM = get(gca,'ylim');    
    text(0.07,0.8,['\xi = ',num2str(slope)],'fontsize',32,'HorizontalAlignment','center','edgecolor','k','Units','normalized')
    if counter == 0
        text(0.5,1.4,'Papa\_summer','HorizontalAlignment','center','FontWeight','b','FontSize',30,'Units','normalized')
    end
    set(gca,'fontsize',28)
    set(gca,'ytick',[1 1e2 1e4 1e6])
    
    if counter == 0
        a=axes('position',get(gca,'position'),'visible','off');
        H = legend(a,[h1(1) h3(1)],'remin OFF','remin ON','Location','south');
        set(H,'fontsize',28,'Position',[0.28 0.73 0.0497 0.0258]);clear H
    end
    
    subplot(3,2,counter*2+2)
    for ii = 5:-1:2
        h1(5-ii+1) = plot(time_submeso,B_submeso(:,ii),'LineWidth',2,'Color',cmap(ii,:),'linestyle','--');
        hold on
        h2(5-ii+1) = plot(time_submeso,B_submeso_remin(:,ii),'LineWidth',2,'Color',cmap(ii,:));
        h3(5-ii+1) = plot(time_submeso,B_submeso_remin(:,ii),'LineWidth',2,'Color',cmap(ii,:));
    end; clear ii
    H = legend(h2,'w_s = 5','w_s = 1','w_s = 0.5','w_s = 0.025','w_s = 0','location','northoutside','Orientation','horizontal');
    set(H,'fontsize',28);clear H
    set(gca,'yscale','log')
    grid on
    xlabel('Day of Year')
    ylabel({'Exported Biomass','[normalized units]'})
    text(0.9,0.2,['\xi = ',num2str(slope)],'fontsize',32,'HorizontalAlignment','center','edgecolor','k','Units','normalized')
    if counter == 0
        text(0.5,1.4,'Papa\_winter','HorizontalAlignment','center','FontWeight','b','FontSize',30,'Units','normalized')
    end
    set(gca,'fontsize',28)
    set(gca,'ytick',[1 1e2 1e4 1e6])
    if counter == 0
        a=axes('position',get(gca,'position'),'visible','off');
        H = legend(a,[h1(1) h3(1)],'remin OFF','remin ON','Location','South');
        set(H,'fontsize',28);clear H
    end
    counter = counter +1;
end

set(gcf,'color','w')
export_fig -r200 Fig7_biomass_export.png
close all
