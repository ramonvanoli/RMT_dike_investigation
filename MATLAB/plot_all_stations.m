dobs = load_data_modem('Datafile.data'); % loads observed data only

% =========================================================
% PLOT OPTIONS - set true/false to enable/disable each plot
% =========================================================
plot_impedance   = true;    % plot real/imag of Z components
plot_rho_phase   = false;   % plot apparent resistivity and phases
plot_tipper      = true;    % plot real/imag of tipper components
show_errorbars   = false;   % errorbars get messy with many stations overlaid

show_stations_as_legend = false;  % true = station legend, false = colorbar
% =========================================================


% =========================================================
% PAPER FIGURE STYLE
% =========================================================
panel_bg = [0.92 0.92 0.95];   % ggplot2-style panel gray
grid_col = [1.00 1.00 1.00];   % white gridlines

% Font sizes
tick_fontsize      = 15;
axislabel_fontsize = 19;   % increased
title_fontsize     = 17;
colorbar_fontsize  = 19;   % increased

% Line and marker sizes
line_width   = 1.1;
marker_size  = 3.5;

% =========================================================


ns = dobs.ns;
cmap = parula(ns);   % one color per station

fig = figure(11);
clf;
fig.Color = 'w';

comp_names = {'ZXX','ZXY','ZYX','ZYY'};
comp_idx   = [1,2,3,4];

% =========================================================
% COUNT ROWS
% =========================================================
n_rows = 0;

if plot_impedance
    n_rows = n_rows + 2;
end

if plot_rho_phase
    n_rows = n_rows + 2;
end

if plot_tipper
    n_rows = n_rows + 1;
end

if n_rows == 0
    error('No plots enabled.');
end


% =========================================================
% TILED LAYOUT
% =========================================================
tl = tiledlayout(n_rows, 4, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

row = 0;


% =========================================================
% IMPEDANCE
% =========================================================
if plot_impedance

    % -----------------------------------------------------
    % Real Z
    % -----------------------------------------------------
    row = row + 1;
    ax_row = gobjects(1,4);

    for ic = 1:4

        ax = nexttile((row-1)*4 + ic);
        hold(ax,'on');

        for is = 1:ns

            f_obs = 1 ./ dobs.T;
            yv = real(dobs.Z(:,comp_idx(ic),is));

            if show_errorbars

                errorbar(ax, f_obs, yv, ...
                    real(dobs.Zerr(:,comp_idx(ic),is)), ...
                    'o-', ...
                    'Color', cmap(is,:), ...
                    'MarkerFaceColor', cmap(is,:), ...
                    'MarkerSize', marker_size, ...
                    'CapSize', 2, ...
                    'LineWidth', line_width);

            else

                plot(ax, f_obs, yv, '-o', ...
                    'Color', cmap(is,:), ...
                    'MarkerFaceColor', cmap(is,:), ...
                    'MarkerSize', marker_size, ...
                    'LineWidth', line_width);

            end

        end

        title(ax, comp_names{ic}, ...
            'FontWeight','bold');

        ylabel(ax, 'Re(Z)');
        xlabel(ax, 'Frequency (Hz)');

        set(ax,'XScale','log');

        style_ggplot(ax, panel_bg, grid_col, ...
            tick_fontsize, axislabel_fontsize, title_fontsize);

        ax_row(ic) = ax;

    end

    % Match diagonal components
    match_ylim(ax_row(1), ax_row(4));

    % Match off-diagonal components
    match_ylim(ax_row(2), ax_row(3));


    % -----------------------------------------------------
    % Imaginary Z
    % -----------------------------------------------------
    row = row + 1;
    ax_row = gobjects(1,4);

    for ic = 1:4

        ax = nexttile((row-1)*4 + ic);
        hold(ax,'on');

        for is = 1:ns

            f_obs = 1 ./ dobs.T;
            yv = imag(dobs.Z(:,comp_idx(ic),is));

            if show_errorbars

                errorbar(ax, f_obs, yv, ...
                    imag(dobs.Zerr(:,comp_idx(ic),is)), ...
                    'o-', ...
                    'Color', cmap(is,:), ...
                    'MarkerFaceColor', cmap(is,:), ...
                    'MarkerSize', marker_size, ...
                    'CapSize', 2, ...
                    'LineWidth', line_width);

            else

                plot(ax, f_obs, yv, '-o', ...
                    'Color', cmap(is,:), ...
                    'MarkerFaceColor', cmap(is,:), ...
                    'MarkerSize', marker_size, ...
                    'LineWidth', line_width);

            end

        end

        title(ax, comp_names{ic}, ...
            'FontWeight','bold');

        ylabel(ax, 'Im(Z)');
        xlabel(ax, 'Frequency (Hz)');

        set(ax,'XScale','log');

        style_ggplot(ax, panel_bg, grid_col, ...
            tick_fontsize, axislabel_fontsize, title_fontsize);

        ax_row(ic) = ax;

    end

    % Match diagonal components
    match_ylim(ax_row(1), ax_row(4));

    % Match off-diagonal components
    match_ylim(ax_row(2), ax_row(3));

end


% =========================================================
% APPARENT RESISTIVITY + PHASE
% =========================================================
if plot_rho_phase

    % -----------------------------------------------------
    % Apparent resistivity
    % -----------------------------------------------------
    row = row + 1;
    ax_row = gobjects(1,4);

    for ic = 1:4

        ax = nexttile((row-1)*4 + ic);
        hold(ax,'on');

        for is = 1:ns

            f_obs = 1 ./ dobs.T;
            yv = dobs.rho(:,comp_idx(ic),is);

            if show_errorbars

                errorbar(ax, f_obs, yv, ...
                    dobs.rhoerr(:,comp_idx(ic),is), ...
                    'o-', ...
                    'Color', cmap(is,:), ...
                    'MarkerFaceColor', cmap(is,:), ...
                    'MarkerSize', marker_size, ...
                    'CapSize', 2, ...
                    'LineWidth', line_width);

            else

                loglog(ax, f_obs, yv, '-o', ...
                    'Color', cmap(is,:), ...
                    'MarkerFaceColor', cmap(is,:), ...
                    'MarkerSize', marker_size, ...
                    'LineWidth', line_width);

            end

        end

        title(ax, comp_names{ic}, ...
            'FontWeight','bold');

        ylabel(ax, '\rho_a (\Omega m)');
        xlabel(ax, 'Frequency (Hz)');

        set(ax,'XScale','log','YScale','log');

        style_ggplot(ax, panel_bg, grid_col, ...
            tick_fontsize, axislabel_fontsize, title_fontsize);

        ax_row(ic) = ax;

    end

    % Match diagonal components
    match_ylim(ax_row(1), ax_row(4));

    % Match off-diagonal components
    match_ylim(ax_row(2), ax_row(3));


    % -----------------------------------------------------
    % Phase
    % -----------------------------------------------------
    row = row + 1;
    ax_row = gobjects(1,4);

    for ic = 1:4

        ax = nexttile((row-1)*4 + ic);
        hold(ax,'on');

        for is = 1:ns

            f_obs = 1 ./ dobs.T;
            yv = dobs.pha(:,comp_idx(ic),is);

            if show_errorbars

                errorbar(ax, f_obs, yv, ...
                    dobs.phaerr(:,comp_idx(ic),is), ...
                    'o-', ...
                    'Color', cmap(is,:), ...
                    'MarkerFaceColor', cmap(is,:), ...
                    'MarkerSize', marker_size, ...
                    'CapSize', 2, ...
                    'LineWidth', line_width);

            else

                plot(ax, f_obs, yv, '-o', ...
                    'Color', cmap(is,:), ...
                    'MarkerFaceColor', cmap(is,:), ...
                    'MarkerSize', marker_size, ...
                    'LineWidth', line_width);

            end

        end

        title(ax, comp_names{ic}, ...
            'FontWeight','bold');

        ylabel(ax, 'Phase (°)');
        xlabel(ax, 'Frequency (Hz)');

        set(ax,'XScale','log');

        style_ggplot(ax, panel_bg, grid_col, ...
            tick_fontsize, axislabel_fontsize, title_fontsize);

        ax_row(ic) = ax;

    end

    % Match diagonal components
    match_ylim(ax_row(1), ax_row(4));

    % Match off-diagonal components
    match_ylim(ax_row(2), ax_row(3));

end


% =========================================================
% TIPPER
% Re(Tx) Im(Tx) Re(Ty) Im(Ty)
% =========================================================
if plot_tipper

    row = row + 1;

    tip_specs = {
        1, 'real', 'Tx', 'Re(T)';
        1, 'imag', 'Tx', 'Im(T)';
        2, 'real', 'Ty', 'Re(T)';
        2, 'imag', 'Ty', 'Im(T)';
        };

    for k = 1:4

        comp = tip_specs{k,1};
        fn   = tip_specs{k,2};
        ttl  = tip_specs{k,3};
        ylab = tip_specs{k,4};

        ax = nexttile((row-1)*4 + k);
        hold(ax,'on');

        for is = 1:ns

            f_obs = 1 ./ dobs.T;
            yv = feval(fn, dobs.tip(:,comp,is));

            if show_errorbars

                errorbar(ax, f_obs, yv, ...
                    feval(fn, dobs.tiperr(:,comp,is)), ...
                    'o-', ...
                    'Color', cmap(is,:), ...
                    'MarkerFaceColor', cmap(is,:), ...
                    'MarkerSize', marker_size, ...
                    'CapSize', 2, ...
                    'LineWidth', line_width);

            else

                plot(ax, f_obs, yv, '-o', ...
                    'Color', cmap(is,:), ...
                    'MarkerFaceColor', cmap(is,:), ...
                    'MarkerSize', marker_size, ...
                    'LineWidth', line_width);

            end

        end

        title(ax, ttl, ...
            'FontWeight','bold');

        ylabel(ax, ylab);
        xlabel(ax, 'Frequency (Hz)');

        set(ax,'XScale','log');

        style_ggplot(ax, panel_bg, grid_col, ...
            tick_fontsize, axislabel_fontsize, title_fontsize);

    end

end


% =========================================================
% STATION LEGEND OR COLORBAR
% =========================================================
if show_stations_as_legend

    % -----------------------------------------------------
    % Station legend inside upper-left subplot
    % -----------------------------------------------------

    ax_legend = nexttile(1);
    hold(ax_legend,'on');

    % Create invisible handles for the legend
    legend_handles = gobjects(ns,1);

    for is = 1:ns

        legend_handles(is) = plot(ax_legend, NaN, NaN, '-o', ...
            'Color', cmap(is,:), ...
            'MarkerFaceColor', cmap(is,:), ...
            'MarkerSize', marker_size, ...
            'LineWidth', line_width);

    end

    % Generate station labels
    legend_labels = arrayfun(@(x) ...
        sprintf('Station %d', x), ...
        1:ns, ...
        'UniformOutput', false);

    % Create legend inside upper-left corner
    lgd = legend(ax_legend, ...
        legend_handles, ...
        legend_labels, ...
        'Location', 'northwest', ...
        'FontName', 'Helvetica', ...
        'FontSize', tick_fontsize, ...
        'Box', 'off');

else

    % -----------------------------------------------------
    % Colorbar
    % -----------------------------------------------------

    colormap(cmap);

    cb = colorbar('eastoutside');

    clim([1 ns]);

    cb.Layout.Tile = 'east';

    cb.Label.String = 'Station index';
    cb.Label.FontSize = colorbar_fontsize;
    cb.Label.FontWeight = 'normal';
    cb.Label.FontName = 'Helvetica';

    cb.FontSize = tick_fontsize;
    cb.FontName = 'Helvetica';

    cb.LineWidth = 1.0;

end


% ============================================================
% OPTIONAL OVERALL TITLE
% ============================================================
% title(tl, ...
%     ['All stations (N = ', num2str(ns), ') - observed data only'], ...
%     'FontSize', 16, ...
%     'FontWeight', 'bold');


%% ============================================================
function match_ylim(ax_a, ax_b)
% Forces two axes to share the same y-axis limits, using the
% union of their individual auto ranges.

    yl_a = ylim(ax_a);
    yl_b = ylim(ax_b);

    yl = [ ...
        min(yl_a(1), yl_b(1)), ...
        max(yl_a(2), yl_b(2)) ...
        ];

    ylim(ax_a, yl);
    ylim(ax_b, yl);

end


%% ============================================================
function style_ggplot(ax, panel_bg, grid_col, ...
    tick_fontsize, axislabel_fontsize, title_fontsize)
% Applies a ggplot2-inspired theme optimized for paper figures.

    % ---------------------------------------------------------
    % General axes appearance
    % ---------------------------------------------------------
    ax.Color = panel_bg;
    ax.Box = 'off';
    ax.TickDir = 'out';

    ax.FontName = 'Helvetica';
    ax.FontSize = tick_fontsize;

    ax.XColor = [0.35 0.35 0.35];
    ax.YColor = [0.35 0.35 0.35];

    ax.LineWidth = 1.0;


    % ---------------------------------------------------------
    % Grid
    % ---------------------------------------------------------
    grid(ax, 'on');

    ax.GridColor = grid_col;
    ax.GridAlpha = 1;

    ax.MinorGridColor = grid_col;
    ax.MinorGridAlpha = 0.6;

    ax.GridLineStyle = '-';
    ax.Layer = 'bottom';

    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';


    % ---------------------------------------------------------
    % Component title
    % ---------------------------------------------------------
    title(ax, ax.Title.String, ...
        'FontWeight', 'bold', ...
        'FontSize', title_fontsize, ...
        'Color', [0.2 0.2 0.2], ...
        'FontName', 'Helvetica');


    % ---------------------------------------------------------
    % X-axis label
    % ---------------------------------------------------------
    ax.XLabel.FontSize = axislabel_fontsize;
    ax.XLabel.FontWeight = 'normal';
    ax.XLabel.FontName = 'Helvetica';
    ax.XLabel.Color = [0.15 0.15 0.15];


    % ---------------------------------------------------------
    % Y-axis label
    % ---------------------------------------------------------
    ax.YLabel.FontSize = axislabel_fontsize;
    ax.YLabel.FontWeight = 'normal';
    ax.YLabel.FontName = 'Helvetica';
    ax.YLabel.Color = [0.15 0.15 0.15];

end