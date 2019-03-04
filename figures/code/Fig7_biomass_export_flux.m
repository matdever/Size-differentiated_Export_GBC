clear

% this code import the number of simulated particles located below 100 m in
% Papa_meso and Papa_submeso

%% Start with Papa_summer

% arbitrary biomass unit
B0 = 1;
% Sinking-velocity classes
ws = [0 0.025 0.5 1 5];
% Reference sinking-velocity
ws0 = max(ws);

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
end; clear ii

for ii = 1:5
    if ii ==1
        continue
    elseif ii == 2
        continue
    elseif ii == 3
        continue
    elseif ii == 4
        filenamedown = [path,'1mday_remin/exportTS1mday_remindown.txt'];
        filenameup = [path,'1mday_remin/exportTS1mday_reminup.txt'];
    elseif ii == 5
        filenamedown = [path,'5mday_remin/exportTS5mday_remindown.txt'];
        filenameup = [path,'5mday_remin/exportTS5mday_reminup.txt'];
    end
    
    % Get downward flux
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(filenamedown,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    flux_meso_remin_down(:,ii) = dataArray{:, 1};
    clearvars filenamedown delimiter formatSpec fileID dataArray ans
    
    % Get upward flux
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(filenameup,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    flux_meso_remin_up(:,ii) = dataArray{:, 1};
    
    clearvars filenameup delimiter formatSpec fileID dataArray ans
end; clear ii path

flux_meso = flux_meso_down-flux_meso_up;
flux_meso_remin = flux_meso_remin_down-flux_meso_remin_up;
clear *_up *_down


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
end; clear ii

for ii = 1:5
    if ii ==1
        continue
    elseif ii == 2
        filenamedown = [path,'0025mday_remin/exportTS0025mday_remindown.txt'];
        filenameup = [path,'0025mday_remin/exportTS0025mday_reminup.txt'];
    elseif ii == 3
        filenamedown = [path,'05mday_remin/exportTS05mday_remindown.txt'];
        filenameup = [path,'05mday_remin/exportTS05mday_reminup.txt'];
    elseif ii == 4
        filenamedown = [path,'1mday_remin/exportTS1mday_remindown.txt'];
        filenameup = [path,'1mday_remin/exportTS1mday_reminup.txt'];
    elseif ii == 5
        filenamedown = [path,'5mday_remin/exportTS5mday_remindown.txt'];
        filenameup = [path,'5mday_remin/exportTS5mday_reminup.txt'];
    end
    
    % Get downward flux
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(filenamedown,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    flux_submeso_remin_down(:,ii) = dataArray{:, 1};
    clearvars filenamedown delimiter formatSpec fileID dataArray ans
    
    % Get upward flux
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(filenameup,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    flux_submeso_remin_up(:,ii) = dataArray{:, 1};
    
    clearvars filenameup delimiter formatSpec fileID dataArray ans
end; clear ii path


flux_submeso = flux_submeso_down-flux_submeso_up;
flux_submeso_remin = flux_submeso_remin_down-flux_submeso_remin_up;

clear *_up *_down
%%

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

%%
fullfigure
counter = 0;
for slope = [4 3 2]
    
    for ii = 1:length(ws)
        B_meso(:,ii) = flux_meso(:,ii).*B0.*(ws(ii)/ws0)^((3-slope)/2);
        B_meso_remin(:,ii) = flux_meso_remin(:,ii).*... % Number of model particles exported
            B0*(ws(ii)/ws0).^(3/2).*... % biomass of 1 model particle at t = 0
            (ws_meso(ii,:)'./ws_meso(ii,1)).^(3/2).*... % Change in biomass vs time for 1 model particle
            (ws(ii)/ws0).^(-slope/2); % Convert to real particles using the spectral slope
        
        
        B_submeso(:,ii) = flux_submeso(:,ii).*B0.*(ws(ii)/ws0)^((3-slope)/2);
        B_submeso_remin(:,ii) = flux_submeso_remin(:,ii).*... % Number of model particles exported
            B0*(ws(ii)/ws0).^(3/2).*... % biomass of 1 model particle at t = 0
            (ws_submeso(ii,:)'./ws_submeso(ii,1)).^(3/2).*... % Change in biomass vs time for 1 model particle
            (ws(ii)/ws0).^(-slope/2); % Convert to real particles using the spectral slope
        
        if ii == 5
            B_meso([1:find(B_meso(:,ii)>0,1,'first')-1 find(B_meso(:,ii)>0,1,'last')+1:end],ii) = 1e-10;
            B_meso_remin([1:find(B_meso_remin(:,ii)>0,1,'first')-1 find(B_meso_remin(:,ii)>0,1,'last')+1:end],ii) = 1e-10;
        else
            B_meso(1:find(B_meso(:,ii)>0,1,'first')-1,ii) = 1e-10;
            B_meso_remin(1:find(B_meso_remin(:,ii)>1e-9,1,'first')-1,ii) = 1e-10;  
        end
        
        B_submeso(1:find(B_submeso(:,ii)>0,1,'first')-1,ii) = 1e-10;
        B_submeso_remin(1:find(B_submeso_remin(:,ii)>0,1,'first')-1,ii) = 1e-10;
        
        B_meso(:,ii) = smooth(B_meso(:,ii),8);
        B_meso_remin(:,ii) = smooth(B_meso_remin(:,ii),8);
        B_submeso(:,ii) = smooth(B_submeso(:,ii),8);
        B_submeso_remin(:,ii) = smooth(B_submeso_remin(:,ii),8);
    end; clear ii
    
    subplot(3,2,counter*2+1)
    for ii = 5:-1:2
        h1(5-ii+1) = plot(time_meso+105,B_meso(:,ii)./unique(diff(time_meso)),'LineWidth',2,'Color',cmap(ii,:),'linestyle','--');
        hold on
        h2(5-ii+1) = plot(time_meso+105,B_meso_remin(:,ii)./unique(diff(time_meso)),'LineWidth',2,'Color',cmap(ii,:));
        h3(5-ii+1) = plot(time_meso+105,B_meso_remin(:,ii)./unique(diff(time_meso)),'LineWidth',2,'Color',cmap(ii,:));
    end; clear ii
    H = legend(h2,'w_s = 5','w_s = 1','w_s = 0.5','w_s = 0.025','w_s = 0','location','northoutside','Orientation','horizontal');
    set(H,'fontsize',28);clear H
    set(gca,'yscale','log','ylim',[1e0 3e6])
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
        h1(5-ii+1) = plot(time_submeso,B_submeso(:,ii)./unique(diff(time_submeso)),'LineWidth',2,'Color',cmap(ii,:),'linestyle','--');
        hold on
        h2(5-ii+1) = plot(time_submeso,B_submeso_remin(:,ii)./unique(diff(time_submeso)),'LineWidth',2,'Color',cmap(ii,:));
        h3(5-ii+1) = plot(time_submeso,B_submeso_remin(:,ii)./unique(diff(time_submeso)),'LineWidth',2,'Color',cmap(ii,:));
    end; clear ii
    H = legend(h2,'w_s = 5','w_s = 1','w_s = 0.5','w_s = 0.025','w_s = 0','location','northoutside','Orientation','horizontal');
    set(H,'fontsize',28);clear H
    set(gca,'yscale','log','ylim',[1e0 3e6])
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