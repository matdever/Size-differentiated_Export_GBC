clear

%speed = {'0025','05','1','5'};
speed = {'0025','1','5'};
fullfigure
for ii = 1:length(speed)
    subplot(2,1,1)
    
    %% SUMMER
    % Get data
    filename = ['/Volumes/home/mdever/allparticles/wtotal',num2str(speed{ii}),'mday_summer.csv'];
    M = csvread(filename);
    %M(M(:,2)==max(M(:,2)),2) = NaN;
    
    % Plots data
    plot(M(:,1)*86400,M(:,2)./sum(M(:,2))*100,'r')
    hold on
    set(gca,'fontsize',18)
    xlabel('Vertical velocity [m/day]')
    ylabel('Particles PDF [%]')
    
    % computes stats
    mu = sum(M(:,1).*M(:,2))/sum(M(:,2));
    std = sqrt(sum(M(:,2).*(M(:,1)-mu).^2)/(sum(M(:,2))-1));
    s = (sum(M(:,2).*(M(:,1)-mu).^3)/sum(M(:,2)))/(std^3);
    k = (sum(M(:,2).*(M(:,1)-mu).^4)/sum(M(:,2)))/(std^4) - 3;
    
    xlim([-20 20])
    %ylim([0 60])
    
    % display stats
    %YLIM = get(gca,'ylim');
    %text(-28,0.85*(YLIM(2)-YLIM(1))+YLIM(1),{['\mu = ',num2str(mu*86400,'%2.2f')],...
    %    ['\sigma = ',num2str(std*86400,'%2.2f')],...
    %    ['s = ',num2str(s*86400,'%2.2e')],...
    %    ['\kappa = ',num2str(k*86400,'%2.2e')]},...
    %    'fontsize',21,'color','r')
    
    
end

%% GET MODEL PDF
istart = 200;
ifinish = 400;

thefolder = pwd;
cd /Volumes/mahadevanlab/mathieudever/PSOM_backup/PSOM_outputs/Papa_summer/highrez
timesteps = 24000:1000:36000;
counter = 0;
thew = [];
for tt = timesteps
    counter = counter +1;
    
    w = ncread(['full_',num2str(tt,'%05d'),'.cdf'],'w');
    w(1,:,:) = NaN; w(end,:,:) = NaN;
    
    w = permute(w,[2 1 3]);
    
    %% PLOT HISTOGRAMS
    
    blop = w(istart:ifinish,:,16:end)*86400;
    thew = cat(1,thew,blop(:));
    
end
cd(thefolder)
yyaxis right
histogram(thew,-100:0.1:100,'facecolor',[.8 .8 .8],'edgecolor',[.8 .8 .8],'Normalization','probability')
set(gca,'ycolor',[.3 .3 .3])
ylabel('Model PDF [%]')
title('Summer','fontsize',21)
grid on



for ii = 1:length(speed)
    subplot(2,1,2)
    %% WINTER
    % Get data
    filename = ['/Volumes/home/mdever/allparticles/wtotal',num2str(speed{ii}),'mday_winter.csv'];
    M = csvread(filename);
    [~,I] = max(M(:,2));
    M(I,:) = [];
    
    % Plots data
    plot(M(:,1)*86400,M(:,2)./sum(M(:,2))*100,'b')
    hold on
    set(gca,'fontsize',18)
    ylabel('Particles PDF [%]')
    xlim([-20 20])
    
    % computes stats
    mu = sum(M(:,1).*M(:,2))/sum(M(:,2));
    std = sqrt(sum(M(:,2).*(M(:,1)-mu).^2)/(sum(M(:,2))-1));
    s = (sum(M(:,2).*(M(:,1)-mu).^3)/sum(M(:,2)))/(std^3);
    k = (sum(M(:,2).*(M(:,1)-mu).^4)/sum(M(:,2)))/(std^4) - 3;
    
    % display stats
%     YLIM = get(gca,'ylim');
%     text(14,0.85*(YLIM(2)-YLIM(1))+YLIM(1),{['\mu = ',num2str(mu*86400,'%2.2f')],...
%         ['\sigma = ',num2str(std*86400,'%2.2f')],...
%         ['s = ',num2str(s*86400,'%2.2e')],...
%         ['\kappa = ',num2str(k*86400,'%2.2e')]},...
%         'fontsize',21,'color','b')
    
    %legend('Summer','Winter','location','northoutside','orientation','horizontal')
    
end; clear ii filename

%% GET MODEL PDF
istart = 200;
ifinish = 400;

thefolder = pwd;
cd /Volumes/mahadevanlab/mathieudever/PSOM_backup/PSOM_outputs/Papa_winter/highrez
timesteps = 32000:1000:56000;
counter = 0;
thew = [];
for tt = timesteps
    counter = counter +1;
    
    w = ncread(['full_',num2str(tt,'%05d'),'.cdf'],'w');
    w(1,:,:) = NaN; w(end,:,:) = NaN;
    
    w = permute(w,[2 1 3]);
    
    %% PLOT HISTOGRAMS
    
    blop = w(istart:ifinish,:,16:end)*86400;
    thew = cat(1,thew,blop(:));
    
end
cd(thefolder)
yyaxis right
%histogram(thew,-100:0.1:100,'facecolor',[.8 .8 .8],'edgecolor',[.8 .8 .8],'Normalization','probability')
[N,E] = histcounts(thew,-100:0.1:100,'Normalization','probability');
N(1001) = (N(1000)+N(1002))/2;
bar((E(2:end)+E(1:end-1))/2,N*100,E(2)-E(1),'facecolor',[.8 .8 .8],'edgecolor',[.8 .8 .8])
grid on
set(gca,'ycolor',[.3 .3 .3])
ylabel('Model PDF [%]')
title('Winter','fontsize',21)












%
%
%
%
%
%
%
%
%
%
%
% %%
% clear
%
% speed = {'0025','05','1','5'};
% fullfigure
% for ii = 1:length(speed)
%     subplot(2,2,ii)
%
%     %% SUMMER
%     % Get data
%     filename = ['wtotal',num2str(speed{ii}),'mday_remin_summer.csv'];
%     M = csvread(filename);
%     %M(M(:,2)==max(M(:,2)),2) = NaN;
%
%     % Plots data
%     plot(M(:,1)*86400,M(:,2)./sum(M(:,2))*100,'r')
%     hold on
%     set(gca,'fontsize',18)
%         xlabel('Vertical velocity [m/day]')
%     ylabel('PDF [%]')
%
%     % computes stats
%     mu = sum(M(:,1).*M(:,2))/sum(M(:,2));
%     std = sqrt(sum(M(:,2).*(M(:,1)-mu).^2)/(sum(M(:,2))-1));
%     s = (sum(M(:,2).*(M(:,1)-mu).^3)/sum(M(:,2)))/(std^3);
%     k = (sum(M(:,2).*(M(:,1)-mu).^4)/sum(M(:,2)))/(std^4) - 3;
%
%     xlim([-30 30])
%         ylim([0 60])
%
%     % display stats
%     YLIM = get(gca,'ylim');
%     text(-28,0.85*(YLIM(2)-YLIM(1))+YLIM(1),{['\mu = ',num2str(mu*86400,'%2.2f')],...
%         ['\sigma = ',num2str(std*86400,'%2.2f')],...
%         ['s = ',num2str(s*86400,'%2.2e')],...
%         ['\kappa = ',num2str(k*86400,'%2.2e')]},...
%         'fontsize',21,'color','r')
%
%     %% WINTER
%     % Get data
%     filename = ['wtotal',num2str(speed{ii}),'mday_remin_winter.csv'];
%     M = csvread(filename);
%     [~,I] = max(M(:,2));
%     M(I,:) = [];
%
%     % Plots data
%     yyaxis right
%     plot(M(:,1)*86400,M(:,2)./sum(M(:,2))*100,'b')
%     hold on
%     set(gca,'fontsize',18)
%     ylabel('PDF [%]')
%         ylim([0 1])
%
%     % computes stats
%     mu = sum(M(:,1).*M(:,2))/sum(M(:,2));
%     std = sqrt(sum(M(:,2).*(M(:,1)-mu).^2)/(sum(M(:,2))-1));
%     s = (sum(M(:,2).*(M(:,1)-mu).^3)/sum(M(:,2)))/(std^3);
%     k = (sum(M(:,2).*(M(:,1)-mu).^4)/sum(M(:,2)))/(std^4) - 3;
%
%     % display stats
%     YLIM = get(gca,'ylim');
%     text(14,0.85*(YLIM(2)-YLIM(1))+YLIM(1),{['\mu = ',num2str(mu*86400,'%2.2f')],...
%         ['\sigma = ',num2str(std*86400,'%2.2f')],...
%         ['s = ',num2str(s*86400,'%2.2e')],...
%         ['\kappa = ',num2str(k*86400,'%2.2e')]},...
%         'fontsize',21,'color','b')
%     legend('Summer','Winter','location','northoutside','orientation','horizontal')
%
% end; clear ii filename


set(gcf,'color','w')
export_fig Fig6.png
%close all