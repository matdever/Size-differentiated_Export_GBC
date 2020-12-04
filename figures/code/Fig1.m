clear

%% This code plot the model domain with the wind profile and the het flux 
% meridoinal anomaly in 3 panels
%%
% ii = 0;
% x = ncread(['/Volumes/mahadevanlab/mathieudever/PSOM_backup/PSOM_outputs/Papa_summer/full_',num2str(ii,'%05d'),'.cdf'],'xc')*1000;
% y = ncread(['/Volumes/mahadevanlab/mathieudever/PSOM_backup/PSOM_outputs/Papa_summer/full_',num2str(ii,'%05d'),'.cdf'],'yc')*1000;
% Z = ncread(['/Volumes/mahadevanlab/mathieudever/PSOM_backup/PSOM_outputs/Papa_summer/full_',num2str(ii,'%05d'),'.cdf'],'zc');
% z = squeeze(Z(1,1,:));
% rho = ncread(['/Volumes/mahadevanlab/mathieudever/PSOM_backup/PSOM_outputs/Papa_summer/full_',num2str(ii,'%05d'),'.cdf'],'rho');
% rho(1,:,:) = NaN; rho(end,:,:) = NaN;
% clear ii
% 
% % Import the scale factor for restauration time scale
% filename = '/Volumes/mahadevanlab/mathieudever/PSOM_backup/PSOM_outputs/scalefactor.out';
% delimiter = ','; startRow = 2; formatSpec = '%*s%s%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
% for col=1:length(dataArray)-1
%     raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
% end
% numericData = NaN(size(dataArray{1},1),size(dataArray,2));
% rawData = dataArray{1};
% for row=1:size(rawData, 1)
%     % Create a regular expression to detect and remove non-numeric prefixes and
%     % suffixes.
%     regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
%     try
%         result = regexp(rawData(row), regexstr, 'names');
%         numbers = result.numbers;
%         
%         % Detected commas in non-thousand locations.
%         invalidThousandsSeparator = false;
%         if numbers.contains(',')
%             thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
%             if isempty(regexp(numbers, thousandsRegExp, 'once'))
%                 numbers = NaN;
%                 invalidThousandsSeparator = true;
%             end
%         end
%         % Convert numeric text to numbers.
%         if ~invalidThousandsSeparator
%             numbers = textscan(char(strrep(numbers, ',', '')), '%f');
%             numericData(row, 1) = numbers{1};
%             raw{row, 1} = numbers{1};
%         end
%     catch
%         raw{row, 1} = rawData{row};
%     end
% end
% scalefactor = cell2mat(raw(:, 1));
% clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp;
% 
% save Fig1_data.mat
%%
clear
load dataset/Fig1_data.mat

figure
pcolor(x/1000,y(2:end-1)/1000,repmat(scalefactor,1,length(x)));
shading flat
cmocean('amp')
set(gca,'fontsize',16,'yticklabel',{})
xlabel('x [km]'); %ylabel('y[km]')
hold on
contour(x/1000,y/1000,rho(:,:,end)',1025:.01:1026,'k')
caxis([0 1])
line(get(gca,'xlim'),[y(2)/1000 y(2)/1000],'linewidth',6,'color','k')
line(get(gca,'xlim'),[y(end-1)/1000 y(end-1)/1000],'linewidth',6,'color','k')
%text(max(x)/1000,175,'~','fontsize',60,'HorizontalAlignment','center')
%text(min(x)/1000,175,'~','fontsize',60,'HorizontalAlignment','center')
text(5,15,'(b)','fontsize',18,'FontWeight','b','color','w')
set(gca,'DataAspectRatio',[1 2 1])

axes('position',[0.76 0.1129 0.15 0.8121])
set(gca,'fontsize',16,'xaxislocation','top','yaxislocation','right');%,'ycolor','w')
hold on
box on; grid on
plot(-1/24*y/1000+1/24*y(304)/1000,y/1000,'k'); axis tight
set(gca,'xtick',-5:5:5)
text(-6,15,'(c)','fontsize',18,'FontWeight','b')

axes('position',[0.125 0.1129 0.15 0.8121])
set(gca,'fontsize',16,'xaxislocation','top');%,'ycolor','w')
hold on
box on; grid on
edge = 0.03;
ywindmin = abs(atanh(2e-2 - 1)/edge/pi) ;
ywindmax = y(end)/1000-ywindmin;
WS(1:305) = 0.5*(tanh(0.03*(y(1:305)/1000-ywindmin)*pi)+1);
WS(306:610) = -0.5*(tanh(0.03*(y(306:610)/1000-ywindmax)*pi)-1);
plot(WS,y/1000,'k'); axis tight
set(gca,'xtick',0:.5:1,'xdir','reverse','xlim',[0 1])
text(0.95,15,'(a)','fontsize',18,'FontWeight','b')
ylabel('y [km]')


set(gcf,'color','w')
%export_fig('Fig1_model_domain.png','-r200')
