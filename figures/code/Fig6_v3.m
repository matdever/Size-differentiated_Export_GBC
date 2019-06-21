clear
close all

istart = 200;
ifinish = 400;

%fullfigure

%%
thefolder = pwd;
for plt = 2:-1:1
    if plt == 1
        cd /Volumes/mahadevanlab/mathieudever/PSOM_backup/PSOM_outputs/Papa_summer/highrez
        timesteps = 24000:1000:36000;
    elseif plt == 2
        cd /Volumes/mahadevanlab/mathieudever/PSOM_backup/PSOM_outputs/Papa_winter/highrez
        timesteps = 32000:1000:56000;
    end
    counter = 0;
    thew = [];
    for tt = timesteps
        counter = counter +1;
        
        
        w = ncread(['full_',num2str(tt,'%05d'),'.cdf'],'w');
        w(1,:,:) = NaN; w(end,:,:) = NaN;
        
        x = ncread(['full_',num2str(tt,'%05d'),'.cdf'],'xc')*1000;
        y = ncread(['full_',num2str(tt,'%05d'),'.cdf'],'yc')*1000;
        Z2 = ncread(['full_',num2str(tt,'%05d'),'.cdf'],'zc');
        z = squeeze(Z2(1,1,:));
        [X,Y] = meshgrid(x,y);
        [X2,Y2,Z2] = meshgrid(x,y,squeeze(Z2(1,1,:)));
        X2 = permute(X2,[2 1 3]);
        Y2 = permute(Y2,[2 1 3]);
        Z2 = permute(Z2,[2 1 3]);
        
        w = permute(w,[2 1 3]);
        
        %% PLOT HISTOGRAMS
        
        blop = w(istart:ifinish,:,16:end)*86400;
        thew = cat(1,thew,blop(:));
        if tt == 36000 && plt == 1
            %histogram(thew,-10:.5:10,'Normalization','probability')
            yyaxis right
            [N,EDGES] = histcounts(thew,-10:.1:10,'Normalization','probability');
            h2 = plot((EDGES(2:end)+EDGES(1:end-1))/2,N*100,'color',[0.8500    0.3250    0.0980],'linewidth',3);
            set(gca,'ycolor',[0.8500    0.3250    0.0980])
            ylabel('PDF [%]')
        elseif tt == 56000 && plt == 2
            fullfigure
%            histogram(thew,-10:.5:10,'Normalization','probability')
            [N,EDGES] = histcounts(thew,-10:.1:10,'Normalization','probability');
            N(101) = (N(100)+N(102))/2;
            h1 = plot((EDGES(2:end)+EDGES(1:end-1))/2,N*100,'color',[ 0    0.4470    0.7410],'linewidth',3);
            set(gca,'ycolor',[ 0    0.4470    0.7410])
            hold on
            xlabel('w [m/day]'); ylabel('PDF [%]')
            set(gca,'fontsize',28)%,'yscale','linear','ylim',[1e-6 1e0],'ytick',[1e-6 1e-4 1e-2 1e0])
        end
    end
    
end; clear ii counter

yyaxis left
set(gca,'ycolor',[ 0    0.4470    0.7410])

line([-5 -5],[0 0.7],'linewidth',1.5,'color','k')
line([-1 -1],[0 0.7],'linewidth',1.5,'color','k')
line([-.05 -.05],[0 0.7],'linewidth',1.5,'color','k')
line([-.025 -.025],[0 0.7],'linewidth',1.5,'color','k')

HH = legend([h1 h2],'Summer','Winter','location','northeast');
set(HH,'fontsize',30)

axes('position',get(gca,'position'),...
    'color','none',...
    'xaxislocation','top',...
    'xlim',get(gca,'xlim'),...
    'xtick',[-5 -1 -.05 -0.025],...
    'ytick',[],...
    'tickdir','out',...
    'fontsize',28)

cd(thefolder)
set(gcf,'color','w')
%export_fig Fig6_v3.png