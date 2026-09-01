clear all
data = load_data_modem('fwd.Data.dat');
T      = data.T;
rho    = data.rho;
pha    = data.pha;
rhoerr = data.rhoerr;
phaerr = data.phaerr;
station = 17;

% ---------------------------------------------------------------
% Optional second dataset, plotted gray and BEHIND the main data
% (e.g. observed vs. forward, or a comparison inversion run)
% Leave data2_file = '' to skip plotting a second dataset.
% ---------------------------------------------------------------
data2_file  = '';                 % e.g. 'obs.Data.dat'
plot_second = ~isempty(data2_file);
if plot_second
    data2   = load_data_modem(data2_file);
    T2      = data2.T;
    rho2    = data2.rho;
    pha2    = data2.pha;
    rhoerr2 = data2.rhoerr;
    phaerr2 = data2.phaerr;
    freq2   = 1 ./ T2;
end

% ---------------------------------------------------------------
% Toggle: add extra breathing room around the rho y-axis limits?
% (Phase is always fixed to [-180, 180], see below.)
% ---------------------------------------------------------------
use_yaxis_padding = true;   % set to false to use tight rho axis limits

% ---------------------------------------------------------------
% Toggle: compressed_view = true  -> ONE figure, 2x2 grid:
%           top-left = rho, diagonal components (xx & yy overlaid)
%           top-right = rho, off-diagonal components (xy & yx overlaid)
%           bottom-left = phase, diagonal components
%           bottom-right = phase, off-diagonal components
%         compressed_view = false -> original behavior: ONE figure,
%           2x4 grid, each of the 4 components in its own column.
% ---------------------------------------------------------------
compressed_view = false;

% ---------------------------------------------------------------
% Toggle: plot the rho y-axis on a LINEAR scale instead of log.
% (Phase is unaffected -- it's always linear, fixed to [-180, 180].)
% ---------------------------------------------------------------
plot_linear = true;   % true = linear rho y-axis; false = log rho y-axis (original)

% Convert period (s) to frequency (Hz) to match Python convention
freq = 1 ./ T;

% Style definitions (matching Python COLORS / MARKERS)
comps   = {'xx', 'xy', 'yx', 'yy'};
colors  = {[0.851 0.310 0.118], ...   % #d94f1e
           [0.102 0.427 0.788], ...   % #1a6dc9
           [0.118 0.561 0.306], ...   % #1e8f4e
           [0.439 0.188 0.627]};      % #7030a0
markers = {'o', 's', '^', 'd'};       % MATLAB has no diamond "D" string, use 'd'

% Style for the background/secondary dataset
gray_color = [0.75 0.75 0.75];
gray_alpha = 0.5;   % marker transparency (best-effort, see try/catch below)

% Plot range: lower bound fixed at 1 kHz, upper bound from data
freq_min = 1e3;
freq_max = max(freq) * 1.1;
if plot_second
    freq_max = max(freq_max, max(freq2) * 1.1);
end

% Y-axis padding factor for rho (only used if use_yaxis_padding = true)
rho_pad_factor = 0.5;   % extra "decades" fraction added top/bottom (log scale)

% Fixed phase axis limits (always applied)
pha_ylim = [-180 180];

% Rho y-axis scale string, derived from plot_linear
if plot_linear
    rho_yscale = 'linear';
else
    rho_yscale = 'log';
end

% When plot_linear = true, use these FIXED y-axis ranges instead of
% computing them from the data (diagonal xx/yy vs off-diagonal xy/yx
% typically sit on very different scales).
rho_ylim_lin_diag    = [0 0.2];   % xx, yy
rho_ylim_lin_offdiag = [14 31];   % xy, yx

if compressed_view
    % =================================================================
    % COMPRESSED VIEW: one figure, 2x2 grid. Each subplot overlays the
    % two components of its group (diagonal or off-diagonal) together.
    % =================================================================
    groups(1).cols = [1 4];   % xx, yy
    groups(1).name = 'Diagonal (xx, yy)';
    groups(2).cols = [2 3];   % xy, yx
    groups(2).name = 'Off-diagonal (xy, yx)';

    fig = figure('Color', 'w', 'Position', [100 100 900 700]);

    for g = 1:2
        cols = groups(g).cols;

        % ---------------- Rho subplot ----------------
        ax_rho = subplot(2, 2, g);
        hold(ax_rho, 'on');
        rho_handles = gobjects(1, numel(cols));
        rho_labels  = cell(1, numel(cols));
        all_rho_vals = [];

        for j = 1:numel(cols)
            col  = cols(j);
            colr = colors{col};
            mk   = markers{col};

            rho_v = rho(:, col, station);
            rho_e = rhoerr(:, col, station);
            all_rho_vals = [all_rho_vals; rho_v(:) - rho_e(:); rho_v(:) + rho_e(:)]; %#ok<AGROW>

            if plot_second
                rho_v2 = rho2(:, col, station+1);
                rho_e2 = rhoerr2(:, col, station+1);
                all_rho_vals = [all_rho_vals; rho_v2(:) - rho_e2(:); rho_v2(:) + rho_e2(:)]; %#ok<AGROW>

                h2 = errorbar(ax_rho, freq2, rho_v2, rho_e2, mk, ...
                    'MarkerSize', 6.5, 'MarkerFaceColor', 'none', ...
                    'MarkerEdgeColor', gray_color, 'Color', gray_color, ...
                    'LineWidth', 1.0, 'LineStyle', 'none', 'CapSize', 3);
                drawnow;
                try
                    h2.MarkerHandle.EdgeColorData(4) = uint8(255*gray_alpha);
                    h2.Bar.EdgeColorData(4)          = uint8(255*gray_alpha);
                    h2.Line.ColorData(4)             = uint8(255*gray_alpha);
                catch
                end
            end

            h1 = errorbar(ax_rho, freq, rho_v, rho_e, mk, ...
                'MarkerSize', 6.5, 'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', colr, 'Color', colr, ...
                'LineWidth', 1.2, 'LineStyle', 'none', 'CapSize', 3);
            rho_handles(j) = h1;
            rho_labels{j}  = ['\rho_{' comps{col} '}'];
        end

        set(ax_rho, 'XScale', 'log', 'YScale', rho_yscale);
        title(ax_rho, groups(g).name, 'Color', [0.07 0.07 0.07]);
        if g == 1
            ylabel(ax_rho, '\rho_a [\Omega\cdotm]', 'Color', [0.07 0.07 0.07]);
        end
        legend(ax_rho, rho_handles, rho_labels, 'Location', 'best', 'Box', 'off');
        style_axes(ax_rho, freq_min, freq_max);
        if plot_linear
            if g == 1
                ylim(ax_rho, rho_ylim_lin_diag);
            else
                ylim(ax_rho, rho_ylim_lin_offdiag);
            end
        else
            if use_yaxis_padding
                set_log_ylim(ax_rho, all_rho_vals, rho_pad_factor);
            else
                set_log_ylim(ax_rho, all_rho_vals, 0);
            end
        end

        % ---------------- Phase subplot ----------------
        ax_phi = subplot(2, 2, g + 2);
        hold(ax_phi, 'on');
        pha_handles = gobjects(1, numel(cols));
        pha_labels  = cell(1, numel(cols));

        for j = 1:numel(cols)
            col  = cols(j);
            colr = colors{col};
            mk   = markers{col};

            pha_v = wrap_phase(pha(:, col, station));
            pha_e = phaerr(:, col, station);

            if plot_second
                pha_v2 = wrap_phase(pha2(:, col, station));
                pha_e2 = phaerr2(:, col, station);

                h2p = errorbar(ax_phi, freq2, pha_v2, pha_e2, mk, ...
                    'MarkerSize', 6.5, 'MarkerFaceColor', 'none', ...
                    'MarkerEdgeColor', gray_color, 'Color', gray_color, ...
                    'LineWidth', 1.0, 'LineStyle', 'none', 'CapSize', 3);
                drawnow;
                try
                    h2p.MarkerHandle.EdgeColorData(4) = uint8(255*gray_alpha);
                    h2p.Bar.EdgeColorData(4)          = uint8(255*gray_alpha);
                    h2p.Line.ColorData(4)             = uint8(255*gray_alpha);
                catch
                end
            end

            h1p = errorbar(ax_phi, freq, pha_v, pha_e, mk, ...
                'MarkerSize', 6.5, 'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', colr, 'Color', colr, ...
                'LineWidth', 1.2, 'LineStyle', 'none', 'CapSize', 3);
            pha_handles(j) = h1p;
            pha_labels{j}  = ['\phi_{' comps{col} '}'];
        end

        set(ax_phi, 'XScale', 'log');
        title(ax_phi, groups(g).name, 'Color', [0.07 0.07 0.07]);
        if g == 1
            ylabel(ax_phi, 'Phase [\circ]', 'Color', [0.07 0.07 0.07]);
        end
        xlabel(ax_phi, 'Frequency [Hz]', 'Color', [0.07 0.07 0.07]);
        legend(ax_phi, pha_handles, pha_labels, 'Location', 'best', 'Box', 'off');
        style_axes(ax_phi, freq_min, freq_max);
        ylim(ax_phi, pha_ylim);
        set(ax_phi, 'YTick', -180:60:180);
    end

    sgtitle(fig, ['Station ' num2str(station)], 'FontSize', 13);

    out_svg = sprintf('station_%d_resistivity_phase_compressed.svg', station);
    saveas(fig, out_svg, 'svg');

else
    % =================================================================
    % ORIGINAL VIEW: one figure, 2x4 grid, one component per column.
    % All 4 rho subplots share a common y-axis.
    % =================================================================

    % Pre-pass: gather rho values across all 4 components (and the
    % second dataset, if present) for a common rho y-axis.
    all_rho_vals = [];
    for col = 1:4
        rho_v = rho(:, col, station);
        rho_e = rhoerr(:, col, station);
        all_rho_vals = [all_rho_vals; rho_v(:) - rho_e(:); rho_v(:) + rho_e(:)]; %#ok<AGROW>
        if plot_second
            rho_v2 = rho2(:, col, station+1);
            rho_e2 = rhoerr2(:, col, station+1);
            all_rho_vals = [all_rho_vals; rho_v2(:) - rho_e2(:); rho_v2(:) + rho_e2(:)]; %#ok<AGROW>
        end
    end

    fig = figure('Color', 'w', 'Position', [100 100 1400 560]);
    ax_rho_list = gobjects(1, 4);

    for col = 1:4
        c    = comps{col};
        colr = colors{col};
        mk   = markers{col};

        % ---------------- Rho subplot (top row) ----------------
        ax_rho = subplot(2, 4, col);
        ax_rho_list(col) = ax_rho;
        hold(ax_rho, 'on');

        rho_v  = rho(:, col, station);
        rho_e  = rhoerr(:, col, station);

        if plot_second
            rho_v2 = rho2(:, col, station+1);
            rho_e2 = rhoerr2(:, col, station+1);

            h2 = errorbar(ax_rho, freq2, rho_v2, rho_e2, mk, ...
                'MarkerSize', 6.5, 'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', gray_color, 'Color', gray_color, ...
                'LineWidth', 1.0, 'LineStyle', 'none', 'CapSize', 3);
            drawnow;
            try
                h2.MarkerHandle.EdgeColorData(4) = uint8(255*gray_alpha);
                h2.Bar.EdgeColorData(4)          = uint8(255*gray_alpha);
                h2.Line.ColorData(4)             = uint8(255*gray_alpha);
            catch
            end
        end

        h1 = errorbar(ax_rho, freq, rho_v, rho_e, mk, ...
            'MarkerSize', 6.5, 'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', colr, 'Color', colr, ...
            'LineWidth', 1.2, 'LineStyle', 'none', 'CapSize', 3);

        set(ax_rho, 'XScale', 'log', 'YScale', rho_yscale);
        title(ax_rho, ['\rho_{' c '}'], 'Color', [0.07 0.07 0.07]);
        if col == 1
            ylabel(ax_rho, '\rho_a [\Omega\cdotm]', 'Color', [0.07 0.07 0.07]);
        end
        if col == 4 && plot_second
            legend(ax_rho, [h1 h2], {'forward computed data', 'observed data'}, 'Location', 'best', 'Box', 'off');
        end
        style_axes(ax_rho, freq_min, freq_max);
        % Common y-axis applied to all rho subplots AFTER this loop.

        % ---------------- Phase subplot (bottom row) ----------------
        ax_phi = subplot(2, 4, col + 4);
        hold(ax_phi, 'on');

        pha_v = wrap_phase(pha(:, col, station));
        pha_e = phaerr(:, col, station);

        if plot_second
            pha_v2 = wrap_phase(pha2(:, col, station));
            pha_e2 = phaerr2(:, col, station);

            h2p = errorbar(ax_phi, freq2, pha_v2, pha_e2, mk, ...
                'MarkerSize', 6.5, 'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', gray_color, 'Color', gray_color, ...
                'LineWidth', 1.0, 'LineStyle', 'none', 'CapSize', 3);
            drawnow;
            try
                h2p.MarkerHandle.EdgeColorData(4) = uint8(255*gray_alpha);
                h2p.Bar.EdgeColorData(4)          = uint8(255*gray_alpha);
                h2p.Line.ColorData(4)             = uint8(255*gray_alpha);
            catch
            end
        end

        errorbar(ax_phi, freq, pha_v, pha_e, mk, ...
            'MarkerSize', 6.5, 'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', colr, 'Color', colr, ...
            'LineWidth', 1.2, 'LineStyle', 'none', 'CapSize', 3);

        set(ax_phi, 'XScale', 'log');
        title(ax_phi, ['\phi_{' c '}'], 'Color', [0.07 0.07 0.07]);
        if col == 1
            ylabel(ax_phi, 'Phase [\circ]', 'Color', [0.07 0.07 0.07]);
        end
        xlabel(ax_phi, 'Frequency [Hz]', 'Color', [0.07 0.07 0.07]);
        style_axes(ax_phi, freq_min, freq_max);
        ylim(ax_phi, pha_ylim);
        set(ax_phi, 'YTick', -180:60:180);
    end

    % Apply rho y-axis limits to all 4 rho subplots.
    % - Linear mode: fixed ranges, different for diagonal (xx, yy)
    %   vs off-diagonal (xy, yx) components.
    % - Log mode: one shared range computed from all 4 components,
    %   as before.
    if plot_linear
        for col = 1:4
            if col == 1 || col == 4
                ylim(ax_rho_list(col), rho_ylim_lin_diag);
            else
                ylim(ax_rho_list(col), rho_ylim_lin_offdiag);
            end
        end
    else
        if use_yaxis_padding
            rho_ylim = get_log_ylim(all_rho_vals, rho_pad_factor);
        else
            rho_ylim = get_log_ylim(all_rho_vals, 0);
        end
        for col = 1:4
            ylim(ax_rho_list(col), rho_ylim);
        end
    end

    out_svg = sprintf('station_%d_resistivity_phase_components.svg', station);
    saveas(fig, out_svg, 'svg');
end


% Helper function: apply Python-like axes styling
function style_axes(ax, xmin, xmax)
    set(ax, 'Color', 'w', ...
        'XColor', [0.13 0.13 0.13], ...
        'YColor', [0.13 0.13 0.13], ...
        'GridColor', [0.6 0.6 0.6], ...
        'GridLineStyle', ':', ...
        'GridAlpha', 0.35, ...
        'MinorGridLineStyle', ':', ...
        'MinorGridAlpha', 0.18, ...
        'FontName', 'Helvetica', ...
        'Box', 'on');
    grid(ax, 'on');
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';
    xlim(ax, [xmin xmax]);
end

% Helper function: pad y-limits on a log-scale axis and SET them
% directly on the given axis (used by the compressed_view branch,
% where each subplot has its own independent rho range).
function set_log_ylim(ax, vals, pad_factor)
    ylim(ax, get_log_ylim(vals, pad_factor));
end

% Helper function: compute padded y-limits on a log-scale axis (adds a
% fraction of a decade above/below the data range so points aren't
% flush with the edges). Returns [lo hi] so the same limits can be
% reused across several axes (used by the non-compressed branch to
% give all 4 rho subplots one shared range).
function lims = get_log_ylim(vals, pad_factor)
    vals = vals(isfinite(vals) & vals > 0);
    if isempty(vals)
        lims = [1 10];
        return
    end
    lo = min(vals);
    hi = max(vals);
    if lo == hi
        lo = lo / 2;
        hi = hi * 2;
    end
    log_lo = log10(lo);
    log_hi = log10(hi);
    span   = log_hi - log_lo;
    pad    = max(span * pad_factor, 0.1);  % minimum pad in case span is tiny
    lims = [10^(log_lo - pad), 10^(log_hi + pad)];
end

% Helper function: wrap phase values (in degrees) into the range
% (-180, 180], so all points are shown on a single, fixed phase axis
% regardless of any 360-degree ambiguity in the raw data.
function wrapped = wrap_phase(pha_deg)
    wrapped = mod(pha_deg + 180, 360) - 180;
end

% %% ---------------- Tipper plot ----------------
T   = data.T;
Tzx = data.tip(:,1,station);
Tzy = data.tip(:,2,station);

if plot_second
    Tzx2 = data2.tip(:,1,station);
    Tzy2 = data2.tip(:,2,station);
end

figure()

% Real parts
subplot(2,2,1)
hold on
if plot_second
    scatter(T2, real(Tzx2), 20, gray_color, 'filled', 'MarkerFaceAlpha', gray_alpha)
end
scatter(T, real(Tzx), 20, 'b', 'filled')
set(gca, 'XScale','log')
xlabel('Period (s)'); ylabel('Re(T_{zx})')
title('Re(T_{zx})'); grid on

subplot(2,2,2)
hold on
if plot_second
    scatter(T2, real(Tzy2), 20, gray_color, 'filled', 'MarkerFaceAlpha', gray_alpha)
end
scatter(T, real(Tzy), 20, 'r', 'filled')
set(gca, 'XScale','log')
xlabel('Period (s)'); ylabel('Re(T_{zy})')
title('Re(T_{zy})'); grid on

% Imaginary parts
subplot(2,2,3)
hold on
if plot_second
    scatter(T2, imag(Tzx2), 20, gray_color, 'filled', 'MarkerFaceAlpha', gray_alpha)
end
scatter(T, imag(Tzx), 20, 'b', 'filled')
set(gca, 'XScale','log')
xlabel('Period (s)'); ylabel('Im(T_{zx})')
title('Im(T_{zx})'); grid on

subplot(2,2,4)
hold on
if plot_second
    scatter(T2, imag(Tzy2), 20, gray_color, 'filled', 'MarkerFaceAlpha', gray_alpha)
end
scatter(T, imag(Tzy), 20, 'r', 'filled')
set(gca, 'XScale','log')
xlabel('Period (s)'); ylabel('Im(T_{zy})')
title('Im(T_{zy})'); grid on

sgtitle(['Station ' num2str(station)], 'FontSize', 13)

%%
%% Tipper induction arrows (MT_lab-consistent style)
% Mirrors MT_lab/mtmap.m's update_map 'IV' case:
%   - uses arrow() instead of quiver()
%   - isotropic scaling (single scale_iv factor on both components)
%   - real and imaginary overlaid in one set of axes (color-coded),
%     same as the ivpv_real / ivpv_imag toggle in the GUI

Tzx = data.tip(:,1,station);
Tzy = data.tip(:,2,station);

% Wise convention
ReX = real(Tzx);
ReY = real(Tzy);
ImX = imag(Tzx);
ImY = imag(Tzy);

x = log10(T);
good = ~isnan(x) & ~isnan(ReX);
x   = x(good);
ReX = ReX(good);
ReY = ReY(good);
ImX = ImX(good);
ImY = ImY(good);

% --- isotropic scaling, equivalent to MT_lab's GUI 'scale_iv' parameter ---
% (MT_lab scales both arrow components by the same factor so arrow
% direction/angle is preserved; no separate horizontal stretch)
scale_iv = 1;   % set this like the scale_iv edit box in mtmap GUI
ReX = ReX*scale_iv;  ReY = ReY*scale_iv;
ImX = ImX*scale_iv;  ImY = ImY*scale_iv;

% which components to draw, equivalent to ivpv_real / ivpv_imag checkboxes
plot_real = true;
plot_imag = true;

% colors (colorblind friendly) - role matches real_col / imag_col in mtmap
realColor = [0.00 0.45 0.74];   % blue
imagColor = [0.85 0.33 0.10];   % orange
lwdth     = 1.2;

figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[5 5 16 8]);

ax = axes; hold(ax,'on'); box(ax,'on');
yline(ax,0,'k-','LineWidth',0.8)

% tail points: one row per period, y = 0 (station axis position)
Xtail = [x(:) zeros(numel(x),1)];

% Real-part arrows (tail -> tail+[ReX ReY]), vectorized over all periods,
% exactly the same call pattern as update_map's 'IV' case:
%   arrow(X, X+Y, 'length',5,'facecolor',col,'edgecolor',col,'Linewidth',lwdth)
if plot_real
    arrow(Xtail, Xtail+[ReY ReX], 'length',5, ...
        'facecolor',realColor,'edgecolor',realColor,'Linewidth',lwdth);
end

% Imaginary-part arrows, overlaid in the same axes
if plot_imag
    arrow(Xtail, Xtail+[ImY ImX], 'length',5, ...
        'facecolor',imagColor,'edgecolor',imagColor,'Linewidth',lwdth);
end

ylabel('T_y','FontSize',11)
xlabel('log_{10}(Period [s])','FontSize',11)
title('Tipper Induction Arrows (Wise Convention)', ...
    'FontWeight','bold','FontSize',12)

% simple manual legend, since arrow() doesn't register with legend()
text(min(x)+0.1, 0.42, 'Real',      'Color',realColor,'FontWeight','bold','FontSize',11)
text(min(x)+0.1, 0.34, 'Imaginary', 'Color',imagColor,'FontWeight','bold','FontSize',11)

grid on
ax.GridAlpha = 0.15;
ax.FontSize  = 10;
ax.LineWidth = 1;
ylim([-0.5 0.5])
xlim([min(x)-0.2 max(x)+0.2])

