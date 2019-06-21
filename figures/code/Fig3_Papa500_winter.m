clear
% Create a blue-to-red colorbar for line plot
cmap = jet(16); cmap = [cmap(1,:);cmap(5,:);cmap(12,:);cmap(16,:)];

%% this code creates the figures that synthetize the forcing and density
% field for the Papa500 xperiment

% Load glider data
load('/Users/mathieudever/Documents/EXPORTS/Data/station_Papa_glider/glider_non_gridded.mat')
glider.N2 = 9.81./nanmean(glider.density(:))*(glider.density(3:end,:)-glider.density(1:end-2,:))./(glider.depth(3:end,:)-glider.depth(1:end-2,:));
glider.time = datevec(nanmean(glider.time));
clearvars -except glider cmap


%% STEP 1 - Compute the TKE vs time. cannot use the full resolution, as it is domain integrated
thefolder = pwd;

cd /Volumes/mahadevanlab/mathieudever/PSOM_backup/PSOM_outputs/Papa_winter
% Surface view
depth = 50; % depth to plot in m
lat = 50;   % Latitude of mid-domain (for f)
f = gsw_f(lat);
g = 9.81;

jstart = 201;
jfinish = 402;

%% Step 2 - load the forcing time series

% Heat Fluxes
fileID = fopen('/Users/mathieudever/Google Drive/Publications/size_differentiated_sinking/figures/code/Papa500_winter/ini_heatfluxes.in','r');
dataArray = textscan(fileID, '%f%[^\n\r]');
fclose(fileID);
HF = dataArray{:, 1};
clearvars fileID dataArray ans;

% Windstress
fileID = fopen('/Users/mathieudever/Google Drive/Publications/size_differentiated_sinking/figures/code/Papa500_winter/ini_windstress.in','r');
dataArray = textscan(fileID, '%f%f%[^\n\r]','delimiter',',');
fclose(fileID);
WS = [dataArray{1:end-1}];
clearvars fileID dataArray ans;

% Applies the scaling imposed in the model
HF(1:5*86400/216) = 0;
for ii = 5*86400/216:10*86400/216
    HF(ii) = HF(ii)*(ii-(5*86400/216))/(5*86400/216);
end
WS(1:5*86400/216,:) = 0;
for ii = 5*86400/216:10*86400/216
    WS(ii,:) = WS(ii,:)*(ii-(5*86400/216))/(5*86400/216);
end

% PLOT TIME SERIES OF FORCINGS
fullfigure
subplot(2,4,1:3)
h1 = plot([0:36000]*216/86400,HF(1:36001),'k','linewidth',3);
ylabel('Net Heat Flux [W/m^2]')
xlabel('Day of Year');
set(gca,'fontsize',30); axis tight
ylim([-60 60]); xlim([0 80])
line([50 50],get(gca,'ylim'),'color','k')
grid on
text(25,50,' Model spin-up ','HorizontalAlignment','center','FontSize',26,'FontWeight','b','EdgeColor','k','BackgroundColor','w')
text(65,50,' Particle tracking ','HorizontalAlignment','center','FontSize',26,'FontWeight','b','EdgeColor','k','BackgroundColor','w')
set(gca,'ytick',-60:30:60)
box on
text(78,55,'(a)','fontsize',30,'fontweight','b')

hold on
yyaxis right
h2 = plot([0:36000]*216/86400,WS(1:36001,1),'--k','linewidth',3);
ylabel('Zonal Wind Stress [N/m^2]')
set(gca,'ytick',-0.2:0.1:0.2,'ylim',[-0.2 0.2],'ycolor','k')

line(get(gca,'xlim'),[0 0],'linestyle',':','color','k','linewidth',2)

H = legend([h1 h2],'Surface Heat Flux','Zonal Wind Stress');
set(H,'location','northwest','fontsize',32)

axes('position',get(gca,'position'),...
    'color','none',...
    'xlim',get(gca,'xlim'),...
    'xtick',[0 31 59],...
    'xticklabel',{'January','February','March'},...
    'ylim',get(gca,'ylim'),...
    'ytick',[],...
    'yticklabel',{},...
    'tickdir','out',...
    'fontsize',30,...
    'xaxislocation','top')

%% Step 3 - Plot initial conditions and snapshots

loop_counter = 1;
for ii =  [0 16000 32000 56000]
    
    rho = ncread(['full_',num2str(ii,'%05d'),'.cdf'],'rho');
    rho(1,:,:) = NaN; rho(end,:,:) = NaN;
    
    x = ncread(['full_',num2str(ii,'%05d'),'.cdf'],'xc')*1000;
    y = ncread(['full_',num2str(ii,'%05d'),'.cdf'],'yc')*1000;
    Z = ncread(['full_',num2str(ii,'%05d'),'.cdf'],'zc');
    z = squeeze(Z(1,1,:));
    
    [X,Y] = meshgrid(x,y);
    [X2,Y2,~] = meshgrid(x,y,squeeze(Z(1,1,:)));
    X2 = permute(X2,[2 1 3]);
    Y2 = permute(Y2,[2 1 3]);
    depth_level = find(abs(abs(z)-depth) == min(abs(abs(z)-depth)));
    
    M2x = -g./nanmean(rho(:)).*...
        (rho(3:end,:,:)-rho(1:end-2,:,:))./(X2(3:end,:,:)-X2(1:end-2,:,:));
    M2y = -g./nanmean(rho(:)).*...
        (rho(:,3:end,:)-rho(:,1:end-2,:))./(Y2(:,3:end,:)-Y2(:,1:end-2,:));
    M2 = sqrt(M2x(:,2:end-1,:).^2+M2y(2:end-1,:,:).^2);
    
    N2 = -9.81./nanmean(rho(:)).*...
        (rho(:,:,3:end)-rho(:,:,1:end-2))./(Z(:,:,3:end)-Z(:,:,1:end-2));
    
    %%
    if loop_counter<4
        subplot(2,4,loop_counter+4)
        pcolor(x(2:end-1)/1000,y(2:end-1)/1000,M2(:,:,depth_level)'); shading interp
        CLIM = get(gca,'clim');
        hold on
        contour(x/1000,y/1000,rho(:,:,depth_level)',1025:.01:1026,'k');
        set(gca,'clim',CLIM)
        box on
        set(gca,'fontsize',30)
        xlabel('x [km]')
        ylabel('y [km]')
        caxis([0 3e-7])
        ylim([y(jstart)/1000 y(jfinish)/1000])
        title(['t = ',num2str((ii-8000)/86400*108 + 20),' days'])
        cmocean('amp')
        if loop_counter==1
            text(4,195,'(b)','fontsize',30,'fontweight','b','backgroundcolor','w')
            title(['doy ',num2str((ii)/86400*216)])
        elseif loop_counter==2
            text(4,195,'(c)','fontsize',30,'fontweight','b','backgroundcolor','w')
            title(['doy ',num2str((ii-8000)/86400*108 + 20)])
            set(gca,'yticklabels',{})
            ylabel('')            
        elseif loop_counter==3
            text(4,195,'(d)','fontsize',30,'fontweight','b','backgroundcolor','w')
            title(['doy ',num2str((ii-8000)/86400*108 + 20)])
            set(gca,'yticklabels',{})
            ylabel('')
            colorbar
            POS = get(gca,'pos');
            POS (3) = 0.1566;
            set(gca,'pos',POS); clear POS
        end
    end
    
    subplot(2,4,[4 8])
    
     if loop_counter<3
    h(loop_counter) = plot(squeeze(nanmean(nanmean(N2(:,jstart:jfinish,:)))),z(2:end-1),'linewidth',2,'color',cmap(loop_counter,:));
    grid on; hold on
    set(gca,'ylim',[-300 0],'xlim',[-0.5e-5 3.5e-4],'fontsize',30,'xaxislocation','top')
    xlabel('N^2 [s^-^2]'); ylabel('Depth [m]')
    text(2.75e-4,-6,'(e)','fontsize',30,'fontweight','b')
    
     elseif loop_counter == 3
         h(loop_counter) = plot(squeeze(nanmean(nanmean(N2(:,jstart:jfinish,:)))),z(2:end-1),'color',cmap(loop_counter,:),'linewidth',2);
% 
%     ind = find(glider.time(:,2)==2 & glider.time(:,3)>=10);        
%         N2data = interp1(-glider.depth(2:end-1,1),nanmean(glider.N2(:,ind),2),z(2:end-1));
%         %theotherh(1) = plot(nanmean(glider.N2(:,ind),2),-glider.depth(2:end-1,1),'--k');
%         theotherh(1) = plot(N2data,z(2:end-1),'color',[0.4940 0.1840 0.5560],'LineStyle','--','linewidth',2);
%         blop = squeeze(nanmean(nanmean(N2(:,jstart:jfinish,:))));
%         good = find(isnan(N2data)==0 & isnan(blop) == 0);%  & z(2:end-1)>-300);
%         r6 = corrcoef(N2data(good),blop(good));
        
    elseif loop_counter == 4
        h(loop_counter) = plot(squeeze(nanmean(nanmean(N2(:,jstart:jfinish,:)))),z(2:end-1),'color',cmap(loop_counter,:),'linewidth',2);

        ind = find(glider.time(:,2)==3);% & glider.time(:,3)==15);        
        N2data = interp1(-glider.depth(2:end-1,1),nanmean(glider.N2(:,ind),2),z(2:end-1));
        %theotherh(1) = plot(nanmean(glider.N2(:,ind),2),-glider.depth(2:end-1,1),'--k');
        theotherh(1) = plot(N2data,z(2:end-1),'color',cmap(loop_counter,:),'LineStyle','--','linewidth',2);
        blop = squeeze(nanmean(nanmean(N2(:,jstart:jfinish,:))));
        good = find(isnan(N2data)==0 & isnan(blop) == 0);% & z(2:end-1)>-300);
        r7 = corrcoef(N2data(good),blop(good));
        
    end
    
    loop_counter = loop_counter + 1;
end
H = legend([h theotherh],'doy 0','doy 30','doy 50','doy 80',['Glider (March; r = ',num2str(r7(1,2),'%1.2f'),')']);
set(H,'fontsize',32,'location','southeast')

POS = get(gca,'position');
set(gca,'position',[POS(1)+0.05 POS(2) POS(3) POS(4)])

cd(thefolder)
set(gcf,'color','w')
%export_fig('Fig3_Papa500_winter.png','-r200')
%close all