clear
%% This code extracts all the necessary data from the SQLite database to
% produce Fig7 of the GBC manuscript.
% This figure plots the histogram of the relative biomass vs vertical
% velocity at DOY 25 in Summer/Winter for xi = 2 and xi = 4

% Add the local path with the importw function
addpath('/vortexfs1/home/mdever/Fig7')  

%% Import w for each size class in winter and summer

counter = 1;
for tt = 1:0.125:25
    disp(['Day ',num2str(tt),' out of 25 (',num2str(counter/193*100),'%)'])
    if mod(59.5+tt,1)==0
        summernum = num2str(59.5+tt,'%2.1f');
    else
        summernum = num2str(59.5+tt);
    end
    [wt0025s(:,counter)] = importw(['/vortexfs1/share/mahadevanlab/mathieudever/particle_backup/particlespeed/Papa_summer/w0025mday_remin_summer_doy',summernum,'.csv']);
    [wt1s(:,counter)] = importw(['/vortexfs1/share/mahadevanlab/mathieudever/particle_backup/particlespeed/Papa_summer/w1mday_remin_summer_doy',summernum,'.csv']);
    [wt5s(:,counter)] = importw(['/vortexfs1/share/mahadevanlab/mathieudever/particle_backup/particlespeed/Papa_summer/w5mday_remin_summer_doy',summernum,'.csv']);
    
    if mod(39+tt,1)==0
        winternum = num2str(39+tt,'%2.1f');
    else
        winternum = num2str(39+tt);
    end    
    [wt0025w(:,counter)] = importw(['/vortexfs1/share/mahadevanlab/mathieudever/particle_backup/particlespeed/Papa_winter/w0025mday_remin_winter_doy',winternum,'.csv']);
    [wt1w(:,counter)] = importw(['/vortexfs1/share/mahadevanlab/mathieudever/particle_backup/particlespeed/Papa_winter/w1mday_remin_winter_doy',winternum,'.csv']);
    [wt5w(:,counter)] = importw(['/vortexfs1/share/mahadevanlab/mathieudever/particle_backup/particlespeed/Papa_winter/w5mday_remin_winter_doy',winternum,'.csv']);
    counter = counter +1;
end; clear tt counter summernum winternum

%% Histogram on total time series
% Compute histograms in summer
[N0025s,E] = histcounts(wt0025s(:)*86400,-10.05:.1:10.05);
[N1s,~] = histcounts(wt1s(:)*86400,-10.05:.1:10.05);
[N5s,~] = histcounts(wt5s(:)*86400,-10.05:.1:10.05);
E = (E(1:end-1)+E(2:end))/2;

N0025s(N0025s==0) = 1e-2; % to lok better in log plots
N1s(N1s==0) = 1e-2;
N5s(N5s==0) = 1e-2;

% Compute histograms in winter
[N0025w,Ew] = histcounts(wt0025w(:)*86400,-155:10:155);
[N1w,~] = histcounts(wt1w(:)*86400,-155:10:155);
[N5w,~] = histcounts(wt5w(:)*86400,-155:10:155);
Ew = (Ew(1:end-1)+Ew(2:end))/2;

N0025w(N0025w==0) = 1e-2;
N1w(N1w==0) = 1e-2;
N5w(N5w==0) = 1e-2;

%% Histogram on each time steps

for ii = 1:size(wt0025s,2)
    % Compute histograms in summer
    [TS.N0025s(:,ii),~] = histcounts(wt0025s(:,ii)*86400,-10.05:.1:10.05);
    [TS.N1s(:,ii),~] = histcounts(wt1s(:,ii)*86400,-10.05:.1:10.05);
    [TS.N5s(:,ii),~] = histcounts(wt5s(:,ii)*86400,-10.05:.1:10.05);
    %E = (E(1:end-1)+E(2:end))/2;
    
    TS.N0025s(TS.N0025s(:,ii)==0,ii) = 1e-2; % to look better in log plots
    TS.N1s(TS.N1s(:,ii)==0,ii) = 1e-2;
    TS.N5s(TS.N5s(:,ii)==0,ii) = 1e-2;
    
    % Compute histograms in winter
    [TS.N0025w(:,ii),~] = histcounts(wt0025w(:,ii)*86400,-155:10:155);
    [TS.N1w(:,ii),~] = histcounts(wt1w(:,ii)*86400,-155:10:155);
    [TS.N5w(:,ii),~] = histcounts(wt5w(:,ii)*86400,-155:10:155);
    %Ew = (Ew(1:end-1)+Ew(2:end))/2;
    
    TS.N0025w(TS.N0025w(:,ii)==0,ii) = 1e-2;
    TS.N1w(TS.N1w(:,ii)==0,ii) = 1e-2;
    TS.N5w(TS.N5w(:,ii)==0,ii) = 1e-2;
end; clear tt

%% Histogram on average of time series
% Compute histograms in summer
N0025s_mean = nanmean(TS.N0025s,2);
N1s_mean = nanmean(TS.N1s,2);
N5s_mean = nanmean(TS.N5s,2);

% Compute histograms in winter
N0025w_mean = nanmean(TS.N0025w,2);
N1w_mean = nanmean(TS.N1w,2);
N5w_mean = nanmean(TS.N5w,2);
%%
% Delete the ws because huge matrix and we can't save it
clear wt*
% Save the rest
save /vortexfs1/home/mdever/Fig7/Fig8_data.mat

%%
