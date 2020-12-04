clear
close all

%% Load gradients computed from glider data for comparison curve
load('dataset/Fig5_data.mat')
glider.M2time = datevec(glider.M2time);
clearvars -except glider cmap


lat = 50;   % Latitude of mid-domain (for f)
f = gsw_f(lat);
g = 9.81;

istart = 200;
ifinish = 400;

%fullfigure

%%
thefolder = pwd;
fullfigure
for plt = 1:2
    if plt == 1
        cd /Volumes/mahadevanlab/mathieudever/PSOM_backup/PSOM_outputs/Papa_summer/highrez
        timesteps = 30500;%[0 24000 30000 36000];
    elseif plt == 2
        cd /Volumes/mahadevanlab/mathieudever/PSOM_backup/PSOM_outputs/Papa_winter
        timesteps = 44000;%[0 32000 44000 56000];
    end
    counter = 0;
    for tt = timesteps
        counter = counter +1;
        
        rho = ncread(['full_',num2str(tt,'%05d'),'.cdf'],'rho');
        rho(1,:,:) = NaN; rho(end,:,:) = NaN;
        u = ncread(['full_',num2str(tt,'%05d'),'.cdf'],'u');
        u(1,:,:) = NaN; u(end,:,:) = NaN;
        v = ncread(['full_',num2str(tt,'%05d'),'.cdf'],'v');
        v(1,:,:) = NaN; v(end,:,:) = NaN;
        w = ncread(['full_',num2str(tt,'%05d'),'.cdf'],'w');
        w(1,:,:) = NaN; w(end,:,:) = NaN;
        % vor = ncread(['full_',num2str(ii,'%05d'),'.cdf'],'vor');
        % vor(1,:,:) = NaN; vor(end,:,:) = NaN;
        
        x = ncread(['full_',num2str(tt,'%05d'),'.cdf'],'xc')*1000;
        y = ncread(['full_',num2str(tt,'%05d'),'.cdf'],'yc')*1000;
        Z2 = ncread(['full_',num2str(tt,'%05d'),'.cdf'],'zc');
        z = squeeze(Z2(1,1,:));
        [X,Y] = meshgrid(x,y);
        [X2,Y2,Z2] = meshgrid(x,y,squeeze(Z2(1,1,:)));
        X2 = permute(X2,[2 1 3]);
        Y2 = permute(Y2,[2 1 3]);
        Z2 = permute(Z2,[2 1 3]);
        
        % Computes Strain
        ux = (u(3:end,:,:)-u(1:end-2,:,:))./(X2(3:end,:,:)-X2(1:end-2,:,:));
        uy = (u(:,3:end,:)-u(:,1:end-2,:))./(Y2(:,3:end,:)-Y2(:,1:end-2,:));
        vx = (v(3:end,:,:)-v(1:end-2,:,:))./(X2(3:end,:,:)-X2(1:end-2,:,:));
        vy = (v(:,3:end,:)-v(:,1:end-2,:))./(Y2(:,3:end,:)-Y2(:,1:end-2,:));
        
        strain = sqrt((ux(:,2:end-1,:)-vy(2:end-1,:,:)).^2 + (vx(:,2:end-1,:)+uy(2:end-1,:,:)).^2);
        vorticity = (vx(:,2:end-1,:)-uy(2:end-1,:,:))./f;
        
        M2x = -g./nanmean(rho(:)).*...
            (rho(3:end,:,:)-rho(1:end-2,:,:))./(X2(3:end,:,:)-X2(1:end-2,:,:));
        M2y = -g./nanmean(rho(:)).*...
            (rho(:,3:end,:)-rho(:,1:end-2,:))./(Y2(:,3:end,:)-Y2(:,1:end-2,:));
        M2 = sqrt(M2x(:,2:end-1,:).^2+M2y(2:end-1,:,:).^2);
        
        %         M2x = -g./nanmean(rho(:)).*...
        %             (rho(12:end,:,:)-rho(1:end-11,:,:))./(X2(12:end,:,:)-X2(1:end-11,:,:));
        %         M2y = -g./nanmean(rho(:)).*...
        %             (rho(:,12:end,:)-rho(:,1:end-11,:))./(Y2(:,12:end,:)-Y2(:,1:end-11,:));
        %         M2 = sqrt(M2x(:,6:end-6,:).^2+M2y(6:end-6,:,:).^2);
        %
        N2 = -9.81./nanmean(rho(:)).*...
            (rho(:,:,3:end)-rho(:,:,1:end-2))./(Z2(:,:,3:end)-Z2(:,:,1:end-2));
        
        [X2,Y2,Z2] = meshgrid(x,y,squeeze(Z2(1,1,:)));
        M2 = permute(M2,[2 1 3]);
        w = permute(w,[2 1 3]);
        vorticity = permute(vorticity,[2 1 3]);
        rho = permute(rho,[2 1 3]);
        
        for ii = 2:size(rho,1)-1
            for jj = 2:size(rho,2)-1
                % Get density profile
                temp = squeeze(rho(ii,jj,:));
                % Interpolate density on 1-m grid
                thedepth = 1:.1:1000;
                rho0 = interp1(-z(isnan(temp)==0),temp(isnan(temp)==0),thedepth);
                MLD(ii,jj) = thedepth(find(rho0>rho0(thedepth==10)+0.03,1,'first'));
            end; clear ii
        end; clear jj
        
        
        %% PLOT M2
        theminus = 2;
        if plt == 1
            subplot(3,3,1)
        elseif plt == 2
            subplot(3,3,3)
        end
        slice(X2(istart:ifinish,2:end-1,1:end-theminus)/1000,Y2(istart:ifinish,2:end-1,1:end-theminus)/1000,Z2(istart:ifinish,2:end-1,1:end-theminus),M2(istart:ifinish,:,1:end-theminus),x(3)/1000,y(istart)/1000,z(end-theminus))
        shading interp; cmocean('amp');
        axis tight;
        zlim([-300 z(end-theminus)+2]);
        set(gca,'fontsize',28)
        xlabel('x [km]','rotation',13);ylabel('y [km]','rotation',-18);zlabel('Z [km]')
        hold on
        plot3(x(1)/1000*ones(size(MLD(istart:ifinish,:),1),1),y(istart:ifinish)/1000,-MLD(istart:ifinish,2),'k','linewidth',3);
        plot3(x(2:end-1)/1000,y(istart)/1000*ones(size(MLD(istart:ifinish,2:end),2),1),-MLD(istart,2:end)','k','linewidth',3);
        if plt == 1
            caxis([0 0.8e-7])
        elseif plt == 2
            caxis([0 3e-7])
        end
        
        POS = get(gca,'position');
        H = colorbar;
        set(gca,'position',POS);
        POS = get(H,'position');
        set(H,'position',[POS(1) POS(2)+0.2*POS(4) POS(3) .55*POS(4)])
        clear POS H
        %contourslice(X2(istart:ifinish,:,1:end-theminus)/1000,Y2(istart:ifinish,:,1:end-theminus)/1000,Z2(istart:ifinish,:,1:end-theminus),rho(istart:ifinish,:,1:end-theminus),x(3)/1000,y(istart)/1000,z(end-theminus))
        
        %% PLOT VORTICITY
        if plt == 1
            subplot(3,3,4)
        elseif plt == 2
            subplot(3,3,6)
        end
        
        slice(X2(istart-1:ifinish,2:end-1,1:end-theminus)/1000,Y2(istart-1:ifinish,2:end-1,1:end-theminus)/1000,Z2(istart-1:ifinish,2:end-1,1:end-theminus),vorticity(istart-1:ifinish,:,1:end-theminus),x(3)/1000,y(istart-1)/1000,z(end-theminus))
        shading interp; cmocean('curl');
        axis tight;
        zlim([-300 z(end-theminus)+2]);
        set(gca,'fontsize',28)
        xlabel('x [km]','rotation',13);ylabel('y [km]','rotation',-18);zlabel('Z [km]')
        hold on     
        if plt == 1
            caxis([-.1 .1])
        elseif plt == 2
            caxis([-1 1])
        end
        POS = get(gca,'position');
        H = colorbar;
        set(gca,'position',POS);
        POS = get(H,'position');
        set(H,'position',[POS(1) POS(2)+0.2*POS(4) POS(3) .55*POS(4)])
        clear POS H
        %contourslice(X2(istart:ifinish,:,1:end-theminus)/1000,Y2(istart:ifinish,:,1:end-theminus)/1000,Z2(istart:ifinish,:,1:end-theminus),rho(istart:ifinish,:,1:end-theminus),x(3)/1000,y(istart)/1000,z(end-theminus))
        % PLOT MLD 
        XLIM = get(gca,'xlim');
        YLIM = get(gca,'ylim');
        plot3(XLIM(1)*ones(size(MLD(istart:ifinish,:),1),1),y(istart:ifinish)/1000,-MLD(istart:ifinish,2),'k','linewidth',3);
        plot3(x(2:end-1)/1000,YLIM(1)*ones(size(MLD(istart:ifinish,2:end),2),1),-MLD(istart,2:end)','k','linewidth',3);   
        
        %% PLOT W
        
        if plt == 1
            subplot(3,3,7)
        elseif plt == 2
            subplot(3,3,9)
        end
        
        slice(X2(istart-1:ifinish,:,1:end-theminus)/1000,Y2(istart-1:ifinish,:,1:end-theminus)/1000,Z2(istart-1:ifinish,:,1:end-theminus),w(istart-1:ifinish,:,1:end-theminus)*86400,x(3)/1000,y(istart-1)/1000,z(end-theminus))
        shading interp; cmocean('balance');
        axis tight;zlim([-300 z(end-theminus)+2]);
        set(gca,'fontsize',28)
        xlabel('x [km]','rotation',13);ylabel('y [km]','rotation',-18);zlabel('Z [km]')
        hold on

        if plt == 1
            caxis([-1 1])
        elseif plt == 2
            caxis([-30 30])
        end
        POS = get(gca,'position');
        H = colorbar;
        set(gca,'position',POS);
        POS = get(H,'position');
        set(H,'position',[POS(1) POS(2)+0.2*POS(4) POS(3) .55*POS(4)])
        clear POS H
        %contourslice(X2(istart:ifinish,:,1:end-theminus)/1000,Y2(istart:ifinish,:,1:end-theminus)/1000,Z2(istart:ifinish,:,1:end-theminus),rho(istart:ifinish,:,1:end-theminus),x(3)/1000,y(istart)/1000,z(end-theminus))
        % PLOT MLD 
        XLIM = get(gca,'xlim');
        YLIM = get(gca,'ylim');
        plot3(XLIM(1)*ones(size(MLD(istart:ifinish,:),1),1),y(istart:ifinish)/1000,-MLD(istart:ifinish,2),'k','linewidth',3);
        plot3(x(2:end-1)/1000,YLIM(1)*ones(size(MLD(istart:ifinish,2:end),2),1),-MLD(istart,2:end)','k','linewidth',3);
        
        %% PLOT HISTOGRAMS
        
        subplot(3,3,2)
        myh(plt) = histogram(M2(istart:ifinish,:,16:end),0:1e-8:5e-7,'Normalization','probability');
        Nmodel = histcounts(M2(istart:ifinish,:,16:end),0:1e-8:5e-7,'Normalization','probability');
        hold on
        xlabel('M2 [s^-^2]'); ylabel('PDF')
        set(gca,'fontsize',28,'yscale','log','ylim',[1e-6 1e0],'ytick',[1e-6 1e-4 1e-2 1e0])
        if plt == 1
            % Display Glider data
            % 2*Rossby radius
            r = 84;
            ind = find(glider.M2time(:,2)==7 & glider.M2dist'>=r/2*1000-1000 & glider.M2dist'<=r/2*1000+1000 & glider.M2depth'<300);
            Ndata = histcounts(glider.M2(ind),0:1e-8:5e-7,'Normalization','probability');
            theh(plt) = plot(1e-8/2:1e-8:5e-7-1e-8/2, Ndata,'color',[0 0.4470 0.7410],'linewidth',3);
            R1 = corrcoef(Nmodel,Ndata);
            
        elseif plt == 2
            POS = get(gca,'position');
            set(gca,'position',[POS(1)+0.01 POS(2)+0.22*POS(4) POS(3) POS(4)-0.38*POS(4)])
            clear POS
            set(myh(plt),'facecolor',[0.8500 0.3250 0.0980]);
            
            % Display Glider data
            % 2*Rossby radius
            r = 30;
            ind = find(glider.M2time(:,2)==2 & glider.M2dist'>=r/2*1000-1000 & glider.M2dist'<=r/2*1000+1000 & glider.M2depth'<300);
            Ndata = histcounts(glider.M2(ind),0:1e-8:5e-7,'Normalization','probability');
            theh(plt) = plot(1e-8/2:1e-8:5e-7-1e-8/2, Ndata,'color',[0.8500 0.3250 0.0980],'linewidth',3);
            R2 = corrcoef(Nmodel,Ndata);

            %XLIM = get(gca,'xlim');
            %theh = plot(-2:-1,[0 0],'k','linewidth',3);
            %set(gca,'xlim',XLIM); clear XLIM
            HH = legend([myh theh],'Papa\_meso','Papa\_submeso',['Glider (Jul; r = ',num2str(R1(1,2),'%1.2f'),')'], ['Glider (Feb; r = ',num2str(R2(1,2),'%1.2f'),')']);
            clear theh
        end
        
        subplot(3,3,5)
        blop = vorticity(istart:ifinish,:,16:end);
        %skewness(blop(:))
        histogram(blop(:),-2:.05:2,'Normalization','probability')
        hold on
        xlabel('\zeta /f'); ylabel('PDF')
        set(gca,'fontsize',28,'yscale','log','ylim',[1e-6 1e0],'ytick',[1e-6 1e-4 1e-2 1e0])
        if plt == 2
            %text(1.22,1e5,['\sigma = ',num2str(nanstd(blop(:)),'%.2f')],'fontsize',20,'color',[0.8500    0.3250    0.0980])
            POS = get(gca,'position');
            set(gca,'position',[POS(1)+0.01 POS(2)+0.22*POS(4) POS(3) POS(4)-0.38*POS(4)])
            clear POS
        elseif plt == 1
            %text(-2,1e5,['\sigma = ',num2str(nanstd(blop(:)),'%.2f')],'fontsize',20,'color',[0    0.4470    0.7410])
        end
        
        subplot(3,3,8)
        blop = w(istart:ifinish,:,16:end)*86400;
        %skewness(blop(:))        
        histogram(blop(:),[-10e-4:2e-5:10e-4]*86400,'Normalization','probability')
        hold on
        xlabel('w [m/day]'); ylabel('PDF')
        set(gca,'fontsize',28,'yscale','log','ylim',[1e-6 1e0],'ytick',[1e-6 1e-4 1e-2 1e0])
        if plt == 2
            %text(65,1e5,['\sigma = ',num2str(nanstd(blop(:)),'%.2f'),' m/day'],'fontsize',20,'color',[0.8500    0.3250    0.0980],'HorizontalAlignment','center')
            POS = get(gca,'position');
            set(gca,'position',[POS(1)+0.01 POS(2)+0.22*POS(4) POS(3) POS(4)-0.38*POS(4)])
            clear POS
        elseif plt == 1
            %text(-65,1e5,['\sigma = ',num2str(nanstd(blop(:)),'%.2f'),' m/day'],'fontsize',20,'color',[0    0.4470    0.7410],'HorizontalAlignment','center')
        end
        
    end; clear ii counter
end

text(-.85,6.1,'Papa\_summer','fontsize',30,'HorizontalAlignment','center','Units','normalized','fontweight','b')
text(1.82,6.1,'Papa\_winter','fontsize',30,'HorizontalAlignment','center','Units','normalized','fontweight','b')

cd(thefolder)
set(gcf,'color','w')
%export_fig -r200 Fig5_dynamic_signature_w_valid.png
%close