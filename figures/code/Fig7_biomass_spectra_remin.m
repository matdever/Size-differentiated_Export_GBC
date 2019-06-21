clear

cmap = get(gca,'ColorOrder');
close all

%speed = {'0025','05','1','5'};
speed = {'0025','1','5'};
%numspeed = [0.025,0.5,1,5];
numspeed = [0.025,1,5];
styles = {'-',':'};
fullfigure
for ii = 2:length(speed)
    subplot(2,1,1)
    
    for xi = [2 4]
        %% SUMMER
        % Get data
        filename = ['/Volumes/home/mdever/sousML/wtotal',speed{ii},'mday_remin_summer_xi',num2str(xi),'.csv'];
        M = csvread(filename);
        %M(M(:,2)==max(M(:,2)),2) = NaN;
        
        % Plots data
        h1(xi/2)=plot(M(:,1)*86400,M(:,2),styles{xi/2},'linewidth',2,'color',cmap(ii,:));
        hold on
        set(gca,'fontsize',18,'yscale','log')
        xlabel('Vertical velocity [m/day]')
        ylabel('Biomass [arbitrary unit]')
        grid on
        
        % computes stats
        mu = sum(M(:,1).*M(:,2))/sum(M(:,2));
        std = sqrt(sum(M(:,2).*(M(:,1)-mu).^2)/(sum(M(:,2))-1));
        s = (sum(M(:,2).*(M(:,1)-mu).^3)/sum(M(:,2)))/(std^3);
        k = (sum(M(:,2).*(M(:,1)-mu).^4)/sum(M(:,2)))/(std^4) - 3;
    end
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
title('Summer','fontsize',21)
legend([h1],'\xi = 2','\xi = 4')

%% WINTER
for ii = 1:length(speed)
    subplot(2,1,2)
    
    for xi = [2 4]
        
        %% SUMMER
        % Get data
        filename = ['/Volumes/home/mdever/sousML/wtotal',speed{ii},'mday_remin_winter_xi',num2str(xi),'.csv'];
        M = csvread(filename);
        
        % Plots data
        h1(ii) = plot(M(:,1)*86400,M(:,2),styles{xi/2},'linewidth',2,'color',cmap(ii,:));
        hold on
        set(gca,'fontsize',18,'yscale','log')
        xlabel('Vertical velocity [m/day]')
        ylabel(' Biomass [arbitrary unit]')
        grid on
        
        % computes stats
        mu = sum(M(:,1).*M(:,2))/sum(M(:,2));
        std = sqrt(sum(M(:,2).*(M(:,1)-mu).^2)/(sum(M(:,2))-1));
        s = (sum(M(:,2).*(M(:,1)-mu).^3)/sum(M(:,2)))/(std^3);
        k = (sum(M(:,2).*(M(:,1)-mu).^4)/sum(M(:,2)))/(std^4) - 3;
    end
    %xlim([-20 20])
    %ylim([0 60])
    
    % display stats
    %YLIM = get(gca,'ylim');
    %text(-28,0.85*(YLIM(2)-YLIM(1))+YLIM(1),{['\mu = ',num2str(mu*86400,'%2.2f')],...
    %    ['\sigma = ',num2str(std*86400,'%2.2f')],...
    %    ['s = ',num2str(s*86400,'%2.2e')],...
    %    ['\kappa = ',num2str(k*86400,'%2.2e')]},...
    %    'fontsize',21,'color','r')
end
title('Winter','fontsize',21)
%legend([h2 h3 h4],'\xi = 2','\xi = 3','\xi = 4')
legend([h1],speed)

set(gcf,'color','w')
export_fig biomass_spectra_remin.png
close all