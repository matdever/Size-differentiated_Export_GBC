clear

cmap = get(gca,'ColorOrder');
%close all

%speed = {'0025','05','1','5'};
speed = {'0025','1','5'};
%numspeed = [0.025,0.5,1,5];
numspeed = [0.025,1,5];
fullfigure
for ii = 1:length(speed)
    subplot(2,1,1)
    
    %% SUMMER
    % Get data
    filename = ['/Volumes/home/mdever/sousML/wtotal',num2str(speed{ii}),'mday_summer.csv'];
    M = csvread(filename);
    %M(M(:,2)==max(M(:,2)),2) = NaN;
    
    % convert Particles to biomass
    wmax = 5; wmin = 0.025;
    slope = 3; B0 = 1;
    Btot = B0/(wmax.^((3-slope)/2))*(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));

    slope = 2;
    B0 = Btot*(wmax.^((3-slope)/2))/(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));
    B2 = B0*(numspeed(ii)/wmax)^((3-slope)/2);
    
    slope = 3;
    B0 = Btot*(wmax.^((3-slope)/2))/(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));
    B3 = B0*(numspeed(ii)/wmax)^((3-slope)/2);
    
    slope = 4;
    B0 = Btot*(wmax.^((3-slope)/2))/(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));
    B4 = B0*(numspeed(ii)/wmax)^((3-slope)/2);
    
    % Plots data
    h1(ii)=plot(M(:,1)*86400,M(:,2)*B2,'-','linewidth',2,'color',cmap(ii,:));
    disp(['Total biomass for speed ',speed{ii},'is ',num2str(sum(M(:,2)*B2)),' and xi = 2'])
    disp(['Total biomass for speed ',speed{ii},'is ',num2str(sum(M(:,2)*B4)),' and xi = 4'])
    hold on
    %plot(M(:,1)*86400,M(:,2)*B3,'--','linewidth',2,'color',cmap(ii,:))
    plot(M(:,1)*86400,M(:,2)*B4,':','linewidth',2,'color',cmap(ii,:))
    set(gca,'fontsize',18,'yscale','log')
    xlabel('Vertical velocity [m/day]')
    ylabel(' Biomass [arbitrary unit]')
    grid on
    
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
title('Summer','fontsize',21)
legend(h1,speed)

%% WINTER
for ii = 1:length(speed)
    subplot(2,1,2)
    
    % Get data
    filename = ['/Volumes/home/mdever/sousML/wtotal',num2str(speed{ii}),'mday_winter.csv'];
    M = csvread(filename);
    [~,I] = max(M(:,2));
    M(I,:) = [];
    
    % convert Particles to biomass
    wmax = 5; wmin = 0.025;
    slope = 3; B0 = 1;
    Btot = B0/(wmax.^((3-slope)/2))*(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));

    slope = 2;
    B0 = Btot*(wmax.^((3-slope)/2))/(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));
    B2 = B0*(numspeed(ii)/wmax)^((3-slope)/2);
    
    slope = 3;
    B0 = Btot*(wmax.^((3-slope)/2))/(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));
    B3 = B0*(numspeed(ii)/wmax)^((3-slope)/2);
    
    slope = 4;
    B0 = Btot*(wmax.^((3-slope)/2))/(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));
    B4 = B0*(numspeed(ii)/wmax)^((3-slope)/2);
    
    % Plots data
    h2 = plot(M(:,1)*86400,M(:,2)*B2,'-','linewidth',2,'color',cmap(ii,:));
    hold on
    %h3 = plot(M(:,1)*86400,M(:,2)*B3,'--','linewidth',2,'color',cmap(ii,:));
    h4 = plot(M(:,1)*86400,M(:,2)*B4,':','linewidth',2,'color',cmap(ii,:));
    set(gca,'fontsize',18,'yscale','log')
    xlabel('Vertical velocity [m/day]')
    ylabel(' Biomass [arbitrary unit]')
    grid on
    
    disp(['Total biomass for speed ',speed{ii},'is ',num2str(sum(M(:,2)*B2)),' and xi = 2'])
    disp(['Total biomass for speed ',speed{ii},'is ',num2str(sum(M(:,2)*B4)),' and xi = 4'])
    
    
    % computes stats
    mu = sum(M(:,1).*M(:,2))/sum(M(:,2));
    std = sqrt(sum(M(:,2).*(M(:,1)-mu).^2)/(sum(M(:,2))-1));
    s = (sum(M(:,2).*(M(:,1)-mu).^3)/sum(M(:,2)))/(std^3);
    k = (sum(M(:,2).*(M(:,1)-mu).^4)/sum(M(:,2)))/(std^4) - 3;
    
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
legend([h2 h4],'\xi = 2','\xi = 4')

set(gcf,'color','w')
export_fig biomass_spectra.png
close all