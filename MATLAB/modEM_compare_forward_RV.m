% This code is used to plot different data sets (observed or predicted) so that they can be
% compared
% ── USER INPUT ────────────────────────────────────────────────────────────────
% Set N_datasets to 1, 2, 3, or 4 — only the first N entries below are used.
N_datasets =2;

% Define up to 4 datasets. Unused entries are ignored.
ds(1).file  = 'fwd.Data.dat';
ds(1).label = 'Forward Data';
ds(1).color = [0.00 0.45 0.70];  % blue

ds(2).file  = 'Datafile_raw.data';
ds(2).label = 'Real Data';
ds(2).color = [0.47 0.67 0.19];  % green 

ds(3).file = '0_5mgrid/final_model_with_river/fwd.Finaldata05_with_river_controlled.dat';
ds(3).label = '05m grid river';
ds(3).color = [0.85 0.33 0.10];  % orange

ds(4).file  = '2mgrid/final_model_with_river/fwd.Finaldata_with_river_controlled.dat';
ds(4).label = '2mgrid river';
ds(4).color = [0.49 0.18 0.56];  % purple

ds(5).file  = 'fwd.Data_padding25.dat';
ds(5).label = 'X = 6000m';
ds(5).color = [0.30 0.75 0.93];  % choose any RGB color

station = 20;
ms      = 20;  % marker size
% ─────────────────────────────────────────────────────────────────────────────

% Load only the datasets we need
for n = 1:N_datasets
    d = load_data_modem(ds(n).file);
    ds(n).T   = d.T;
    ds(n).rho = d.rho;
    ds(n).pha = d.pha;
end

rho_titles = {'\rho_{xx}', '\rho_{xy}', '\rho_{yx}', '\rho_{yy}'};
pha_titles = {'\phi_{xx}', '\phi_{xy}', '\phi_{yx}', '\phi_{yy}'};

figure('Position', [100 100 1400 600])

for k = 1:4

    % ── Apparent resistivity ─────────────────────────────────────────────────
    subplot(2, 4, k)
    for n = 1:N_datasets
        scatter(ds(n).T, ds(n).rho(:,k,station), ms, ds(n).color, ...
            'filled', 'DisplayName', ds(n).label); hold on
    end
    set(gca, 'XScale', 'log', 'YScale', 'log')
    xlabel('Period (s)'); ylabel('\rho (\Omegam)')
    title(rho_titles{k}); grid on
    if k == 1, legend('Location', 'best'); end

    % ── Phase ────────────────────────────────────────────────────────────────
    subplot(2, 4, k+4)
    for n = 1:N_datasets
        scatter(ds(n).T, ds(n).pha(:,k,station), ms, ds(n).color, ...
            'filled', 'DisplayName', ds(n).label); hold on
    end
    set(gca, 'XScale', 'log')
    xlabel('Period (s)'); ylabel('Phase (°)')
    title(pha_titles{k}); grid on

end

sgtitle(['Station ' num2str(station)], 'FontSize', 13)