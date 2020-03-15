clear
% Set of marameters
wmax = 5; wmin = 0.025;
slope = 3; B0 = 1;

% Compute total biomass based on the paramters
Btot = B0/(wmax.^((3-slope)/2))*(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));
B0 = B0/Btot;
B0025_3 = B0*(0.025/wmax).^((3-slope)/2);
B1_3 = B0*(1/wmax).^((3-slope)/2);
B5_3 = B0*(5/wmax).^((3-slope)/2);

% For Junge slope = 2
slope = 2;
% Compute biomass of 1 largest particle necessary to keep total biomass constant
B0 = Btot*(wmax.^((3-slope)/2))/(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));
B0025_2 = B0*(0.025/wmax).^((3-slope)/2);
B1_2 = B0*(1/wmax).^((3-slope)/2);
B5_2 = B0*(5/wmax).^((3-slope)/2);

% For Junge slope = 4
slope = 4;
% Compute biomass of 1 largest particle necessary to keep total biomass constant
B0 = Btot*(wmax.^((3-slope)/2))/(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));
B0025_4 = B0*(0.025/wmax).^((3-slope)/2);
B1_4 = B0*(1/wmax).^((3-slope)/2);
B5_4 = B0*(5/wmax).^((3-slope)/2);

% For Junge slope = 5
slope = 5;
% Compute biomass of 1 largest particle necessary to keep total biomass constant
B0 = Btot*(wmax.^((3-slope)/2))/(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));
B0025_5 = B0*(0.025/wmax).^((3-slope)/2);
B1_5 = B0*(1/wmax).^((3-slope)/2);
B5_5 = B0*(5/wmax).^((3-slope)/2);

% For Junge slope = 6
slope = 6;
% Compute biomass of 1 largest particle necessary to keep total biomass constant
B0 = Btot*(wmax.^((3-slope)/2))/(2/(5-slope)*wmax.^((5-slope)/2) - 2/(5-slope)*wmin.^((5-slope)/2));
B0025_6 = B0*(0.025/wmax).^((3-slope)/2);
B1_6 = B0*(1/wmax).^((3-slope)/2);
B5_6 = B0*(5/wmax).^((3-slope)/2);
%%

% Set colors for each class
color0025 = [215,48,39]/256;
color1 = [39,100,25]/256;
color5 = [49,54,149]/256;

figure(2)
thew = [0.025 1 5];
h1 = loglog(thew,[B0025_3 B1_3 B5_3]./sum([B0025_3 B1_3 B5_3])*100,'-k','linewidth',2);
hold on
h2 = loglog(thew,[B0025_2 B1_2 B5_2]./sum([B0025_2 B1_2 B5_2])*100,':k','linewidth',2);
h3 = loglog(thew,[B0025_4 B1_4 B5_4]./sum([B0025_4 B1_4 B5_4])*100,'--k','linewidth',2);


scatter(thew,[B0025_3 B1_3 B5_3]./sum([B0025_3 B1_3 B5_3])*100,100,cat(1,color0025,color1,color5),'s','filled');
scatter(thew,[B0025_2 B1_2 B5_2]./sum([B0025_2 B1_2 B5_2])*100,100,cat(1,color0025,color1,color5),'s','filled');
scatter(thew,[B0025_4 B1_4 B5_4]./sum([B0025_4 B1_4 B5_4])*100,100,cat(1,color0025,color1,color5),'s','filled');


grid on
set(gca,'fontsize',16,'ytick',[0 10 100],'ylim',[0 100],'yticklabel',{'0','10','100'})
xlabel('Sinking velocity [m/day]')
ylabel('Relative Biomass content [%]')
legend([h2 h1 h3],'\xi = 2','\xi = 3','\xi = 4','location','north')
%s1 = gca;
%h = findobj(gcf,'type','axes') % Find the axes object in the GUI
%s = copyobj(h,gcf) % Copy axes object h into figure f1
%set(s,'xaxislocation','top','color','none','XLabel',' ','YLabel',' ')


% get the figure smaller to accomodate the 2nd x-axis later
POS = get(gca,'position');
set(gca,'position',[POS(1) POS(2) POS(3) POS(4)-0.05])

% create 2 x-axis
axes('position',get(gca,'position'),...
    'xlim',get(gca,'xlim'),...
    'xscale','log',...
    'ylim',get(gca,'ylim'),...
    'ytick',get(gca,'ytick'),...
    'yticklabel',{},...
    'yscale','log',...
    'fontsize',16,...
    'color','none',...
    'XAxisLocation','top')
hold on
% Parameters for Stokes sinking velocity equation
g = 9.81;                   % m/s2
rho_water = 1026;           % kg/m3
rho_part = 1.01*rho_water;   % kg/m3
viscosity = 1.48e-3;        % Pa.s

% If I want to set ticks at specific radii
% myxtick = [0:10:80];
% wseq = (2.*(rho_part-rho_water).*g.*(myxtick/1e6).^2)./(9.*viscosity)*86400;
% set(gca,'xtick',wseq,'xticklabel',myxtick)

% If I want to compute the radius at the velocities already ticked on the
% 1st x-axis
Xtick = get(gca,'xtick');
Xtick = [0.02 0.05 0.1 0.5 1 5 10];
for ii = 1:length(Xtick)
    
    % equivalent stokes' radius
    xticklab{ii} = num2str(sqrt((9*viscosity*Xtick(ii)/86400)./(2*g*(rho_part-rho_water)))*1e6,'%.0f');
    
end

set(gca,'xtick',Xtick,'xticklabel',xticklab)
xlabel('Equivalent Stokes radius [\mum]')
set(gcf,'color','w')
%export_fig -r300 Fig4a.png
%%


%% This code plots the distribution of Particle, given a spectral slope and the number of LARGEST particle
clear
% Set colors for each class
color0025 = [215,48,39]/256;
color1 = [39,100,25]/256;
color5 = [49,54,149]/256;
counter = 1;
for spectral_slope = [2 3 4]
    
    % compute for the sinking velocity actually modeled
    ws = [0.025 1 5];
    N0 = 1971717;
    N = (ws/max(ws)).^-(spectral_slope/2)*N0;
    B = (ws/max(ws)).^((3-spectral_slope)/2);
    
    % Create the theoretical curve
    thexaxis = ws;
    theyaxis = (ws/max(ws)).^-(spectral_slope/2);
    theyaxis./sum(theyaxis)*100
    %
    
    figure(1)
    % theoretical curve
    if spectral_slope>3
        style = '--';
    elseif spectral_slope == 3
        style = '-';
    elseif spectral_slope <3
        style = ':';
    end

    % Plot N/N0
    h(counter) = plot(ws,theyaxis./sum(theyaxis)*100,'linestyle',style,'color','k','linewidth',2);
    hold on
    scatter(ws,theyaxis./sum(theyaxis)*100,100,cat(1,color0025,color1,color5),'s','filled')
    
    set(gca,'yscale','log','xscale','log','fontsize',16)
    grid on; axis tight;    
    
    xlabel('Sinking velocity [m/day]')
    %ylabel('Number of particles')
    xlim([0.01 10])
    counter = counter +1;
end
legend(h,'\xi = 2','\xi = 3','\xi = 4','location','northeast')
ylabel('Relative number of particles [%]')

% get the figure smaller to accomodate the 2nd x-axis later
POS = get(gca,'position');
set(gca,'position',[POS(1) POS(2) POS(3) POS(4)-0.05],'ytick',[0.01 0.1 1 10 100],'Yticklabel',{'0.01','0.1','1','10','100'}, 'ylim',[0.001 100])

% create 2 x-axis
axes('position',get(gca,'position'),...
    'xlim',get(gca,'xlim'),...
    'xscale','log',...
    'ylim',get(gca,'ylim'),...
    'ytick',get(gca,'ytick'),...
    'yticklabel',{},...
    'yscale','log',...
    'fontsize',16,...
    'color','none',...
    'XAxisLocation','top')
hold on
% Parameters for Stokes sinking velocity equation
g = 9.81;                   % m/s2
rho_water = 1026;           % kg/m3
rho_part = 1.01*rho_water;   % kg/m3
viscosity = 1.48e-3;        % Pa.s

% If I want to set ticks at specific radii
% myxtick = [0:10:80];
% wseq = (2.*(rho_part-rho_water).*g.*(myxtick/1e6).^2)./(9.*viscosity)*86400;
% set(gca,'xtick',wseq,'xticklabel',myxtick)

% If I want to compute the radius at the velocities already ticked on the
% 1st x-axis
Xtick = get(gca,'xtick');
Xtick = [0.02 0.05 0.1 0.5 1 5 10];
for ii = 1:length(Xtick)
    
    % equivalent stokes' radius
    xticklab{ii} = num2str(sqrt((9*viscosity*Xtick(ii)/86400)./(2*g*(rho_part-rho_water)))*1e6,'%.0f');
    
end

set(gca,'xtick',Xtick,'xticklabel',xticklab)
xlabel('Equivalent Stokes radius [\mum]')
set(gcf,'color','w')
%export_fig -r300 Fig4b.png