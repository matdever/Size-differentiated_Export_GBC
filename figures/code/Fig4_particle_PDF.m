clear

%% This code plots the distribution of Particle, given a spectral slope and the number of LARGEST particle

counter = 1;
for spectral_slope = [2 3 4]
    
    % compute for the sinking velocity actually modeled
    ws = [0.025 0.5 1 5];
    N0 = 1971717;
    N = (ws/max(ws)).^-(spectral_slope/2)*N0;
    B = (ws/max(ws)).^((3-spectral_slope)/2);
    
    % Create the theoretical curve
    thexaxis = [min(ws):0.001:max(ws)];
    theyaxis = (thexaxis/thexaxis(4976)).^-(spectral_slope/2);
    theyaxisB = (thexaxis/thexaxis(4976)).^((3-spectral_slope)/2);
    
    %%
    
    figure(1)
    % theoretical curve
    if spectral_slope>3
        color = 'r';
    elseif spectral_slope == 3
        color = 'k';
    elseif spectral_slope <3
        color = 'b';
    end

    % Plot N/N0
    h(1) = plot(thexaxis,theyaxis,'--','color',color,'linewidth',2);
    hold on
    scatter(ws,N./N0,100,'s','markeredgecolor',color,'markerfacecolor',color)
    
    % Plot B/B0
    h(2) = plot(thexaxis,theyaxisB,'-','color',color,'linewidth',2);
    g(counter) = h(2);
    scatter(ws,B,100,'s','markeredgecolor',color,'markerfacecolor',color)
    
    set(gca,'yscale','log','xscale','log','fontsize',16)
    grid on; axis tight;    
    
    xlabel('Sinking velocity [m/day]')
    %ylabel('Number of particles')
    xlim([0.01 10])
    counter = counter +1;
end
H = legend([h(1) h(2)],'N/N_0','B/B_0','location','east');
set(H,'Position',[0.7554 0.56 0.1321 0.1250])


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
%text(3,3e11,'\zeta = -3.5','fontsize',18)

counter = 1;
for spectral_slope = [4 3 2]
    
    % compute for the sinking velocity actually modeled
    ws = [0.025 0.5 1 5];
    N0 = 1971717;
    N = (ws/max(ws)).^-(spectral_slope/2)*N0;
    B = (ws/max(ws)).^((3-spectral_slope)/2);
    
    % Create the theoretical curve
    thexaxis = [min(ws):0.001:max(ws)];
    theyaxis = (thexaxis/thexaxis(4976)).^-(spectral_slope/2);
    theyaxisB = (thexaxis/thexaxis(4976)).^((3-spectral_slope)/2);
    
    %%
    
    figure(1)
    % theoretical curve
    if spectral_slope>3
        color = 'r';
    elseif spectral_slope == 3
        color = 'k';
    elseif spectral_slope <3
        color = 'b';
    end

    % Plot N/N0
    h(1) = plot(thexaxis,theyaxis,'--','color',color,'linewidth',2);
    hold on
    scatter(ws,N./N0,100,'s','markeredgecolor',color,'markerfacecolor',color)
    
    % Plot B/B0
    h(2) = plot(thexaxis,theyaxisB,'-','color',color,'linewidth',2);
    g(counter) = h(2);
    scatter(ws,B,100,'s','markeredgecolor',color,'markerfacecolor',color)
        
    counter = counter +1;
end
H = legend(g,'\xi = 4','\xi = 3','\xi = 2','location','northeast');
set(H,'Color',[1 1 1])


set(gcf,'color','w')
export_fig('Fig4_particle_PDF.png','-r200')
