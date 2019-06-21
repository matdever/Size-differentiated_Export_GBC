clear

%% Import (x,y,z) for each size class in winter and summer
[x0025s,y0025s,z0025s] = importxyz('/Volumes/mdever/particleposition/doy25/wtotal0025mday_summer.csv');
[x1s,y1s,z1s] = importxyz('/Volumes/mdever/particleposition/doy25/wtotal1mday_summer.csv');
[x5s,y5s,z5s] = importxyz('/Volumes/mdever/particleposition/doy25/wtotal5mday_summer.csv');

[x0025,y0025,z0025] = importxyz('/Volumes/mdever/particleposition/doy25/wtotal0025mday_winter.csv');
[x1,y1,z1] = importxyz('/Volumes/mdever/particleposition/doy25/wtotal1mday_winter.csv');
[x5,y5,z5] = importxyz('/Volumes/mdever/particleposition/doy25/wtotal5mday_winter.csv');

color0025 = [215,48,39]/256;
color1 = [39,100,25]/256;
color5 = [49,54,149]/256;

%% Bin the data into 1x1 km and compute the median vertical displacement
% with respect to the gravitational model
thedoy = 25;
countx = 1;
rez = 1;
for x = 0:rez:111
    county = 1;
    for y = 100:rez:200
        
        %% 0025 m/day
        ind = find(x0025>=x & x0025<=x+rez & y0025>=y & y0025<=y+rez);
        if isempty(ind)
            medi0025(countx,county) = NaN;
        else
            medi0025(countx,county) = median(z0025(ind)-(80+0.025*thedoy));
        end
        
        %% 1 m/day
        ind = find(x1>=x & x1<=x+rez & y1>=y & y1<=y+rez);
        if isempty(ind)
            medi1(countx,county) = NaN;
        else
            medi1(countx,county) = median(z1(ind)-(80+1*thedoy));
        end
        
        %% 5 m/day
        ind = find(x5>=x & x5<=x+rez & y5>=y & y5<=y+rez);
        if isempty(ind)
            medi5(countx,county) = NaN;
        else
            medi5(countx,county) = median(z5(ind)-(80+5*thedoy));
        end
        county = county + 1;
        
    end; clear y
    countx = countx + 1;
end; clear x

%%

figure
% 0.025 m/day subplot
subplot(3,2,1)
pcolor(0:rez:111,100:rez:200,medi0025'); shading interp
cmocean('balance'); caxis([-50 50])
set(gca,'fontsize',14,'xticklabels',{},'tickdir','out')
H = colorbar; set(H,'location','northoutside','TickLength',0.04)
hold on
contour(0:rez:111,100:rez:200,medi0025',[0 0],'k')
axis tight; ylabel('y [km]')
set(gca,'position',[0.1 .66 .33 .25])
set(H,'position',[.12 .92 .3 .015])
text(109,188,'0.025 m/day','BackgroundColor','w','EdgeColor','k','fontsize',14,'HorizontalAlignment','right','color',color0025)

subplot(3,2,3)
pcolor(0:rez:111,100:rez:200,medi1'); shading interp
cmocean('balance'); caxis([-50 50])
set(gca,'fontsize',14,'xticklabels',{},'tickdir','out')
hold on
contour(0:rez:111,100:rez:200,medi1',[0 0],'k')
axis tight; ylabel('y [km]')
set(gca,'position',[0.1 .365 .33 .25])
text(109,188,'1 m/day','BackgroundColor','w','EdgeColor','k','fontsize',14,'HorizontalAlignment','right','color',color1);

subplot(3,2,5)
pcolor(0:rez:111,100:rez:200,medi5'); shading interp
cmocean('balance'); caxis([-50 50])
set(gca,'fontsize',14,'tickdir','out')
hold on
contour(0:rez:111,100:rez:200,medi5',[0 0],'k')
axis tight; ylabel('y [km]')
set(gca,'position',[0.1 .07 .33 .25])
text(109,188,'5 m/day','BackgroundColor','w','EdgeColor','k','fontsize',14,'HorizontalAlignment','right','color',color5);

% Compute PDFs
Z = 0:.5:300;
[N0025s,~] = histcounts(z0025s,Z,'Normalization','probability');
[N1s,~] = histcounts(z1s,Z,'Normalization','probability');
[N5s,~] = histcounts(z5s,Z,'Normalization','probability');
[N0025,~] = histcounts(z0025,Z,'Normalization','probability');
[N1,~] = histcounts(z1,Z,'Normalization','probability');
[N5,~] = histcounts(z5,Z,'Normalization','probability');
Z = (Z(1:end-1)+Z(2:end))/2;

% plot PDFs
subplot(3,2,[2 4 6])
h1 = patch([N0025s 0 0 N0025s(1)],-[Z Z(end) Z(1) Z(1)],'k','edgecolor','k','facealpha',0.1,'edgealpha',.1);
patch([N0025s 0 0 N0025s(1)],[Z Z(end) Z(1) Z(1)],color0025,'edgecolor',color0025,'facealpha',0.1,'edgealpha',.1);
hold on; ylim([0 300])
patch([N1s 0 0 N1s(1)],[Z Z(end) Z(1) Z(1)],color1,'edgecolor',color1,'facealpha',0.1,'edgealpha',.1)
patch([N5s 0 0 N5s(1)],[Z Z(end) Z(1) Z(1)],color5,'edgecolor',color5,'facealpha',0.1,'edgealpha',.1)
h2 = plot(N0025,-Z,'color','k','LineWidth',2);
plot(N0025,Z,'color',color0025,'LineWidth',2);
plot(N1,Z,'color',color1,'LineWidth',2);
plot(N5,Z,'color',color5,'LineWidth',2);
xlabel('Number of particles');ylabel('Depth [m]')
set(gca,'XAxisLocation','top','FontSize',14,'YDir','reverse','XTickLabel',{})
grid on; box on
h3 = line(get(gca,'xlim'),[80 80],'color','k','linestyle','--','linewidth',2);
plot(N0025,Z,'color',color0025,'LineWidth',1.5)
plot(N1,Z,'color',color1,'LineWidth',1.5)
plot(N5,Z,'color',color5,'LineWidth',1.5)
set(gca,'Position',[0.6 0.06 0.3 0.85])

% Mark the reference depth for vertical displacement
scatter(max(get(gca,'xlim')),80+0.025*25,100,'v','filled','markerfacecolor',color0025)
scatter(max(get(gca,'xlim')),80+1*25,100,'v','filled','markerfacecolor',color1)
scatter(max(get(gca,'xlim')),80+5*25,100,'v','filled','markerfacecolor',color5)
legend([h1 h2 h3],'Summer','Winter','seeded region','location','southeast')

set(gcf,'color','w')
export_fig -r300 Fig6.png

%% INLINE FUNCTION

function [x,y,z] = importxyz(filename)

%% Initialize variables.
delimiter = ',';
startRow = 1;
endRow = inf;

%% Format for each line of text:
formatSpec = '%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');


%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
x = dataArray{:, 1}/1000;
y = dataArray{:, 2}/1000;
z = -dataArray{:, 3};
end


