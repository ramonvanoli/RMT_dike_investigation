% This script visualizes the inversion results and data fit of ModEM
% inversions and is based on mtcode functions
m = load_model_modem('DIKE3_NLCG_030.rho'); %loads model either .res or .model
dobs = load_data_modem('Datafile.data'); %loads observed data
dpred = load_data_modem('DIKE3_NLCG_030.dat'); %loads predicted data
[m,dobs] = link_model_data(m,dobs); %geo-references model
plot_diagonal_section(m,dobs); %plot diagonal slice
%plot_model_3D(m,dobs)
% Plot misfit
invlog = load_logfile_modem_RV('DIKE3_NLCG.log');
plot_misfit_convergence(invlog)
% compare_data(dobs,dpred00,dpred26)
s = detailed_statistics(dobs, dpred);
% =========================================================
% PLOT OPTIONS - set true/false to enable/disable each plot
% =========================================================
plot_impedance   = true;   % plot real/imag of Z components
plot_rho_phase   = false;   % plot apparent resistivity and phases
plot_tipper      = true;   % plot real/imag of tipper components
% =========================================================
is = 1;

% ggplot2-inspired palette
col_obs  = [0.20 0.20 0.20];   % dark gray markers for observed
col_pred = [0.00 0.45 0.70];   % muted blue for predicted line
panel_bg = [0.92 0.92 0.95];   % classic ggplot2 panel gray
grid_col = [1.00 1.00 1.00];   % white gridlines

while 1

    fig = figure(10); clf;
    fig.Color = 'w';

    comp_names = {'ZXX','ZXY','ZYX','ZYY'};
    comp_idx   = [1,2,3,4];

    % Count rows
    n_rows = 0;
    if plot_impedance,  n_rows = n_rows + 2; end
    if plot_rho_phase,  n_rows = n_rows + 2; end
    if plot_tipper,     n_rows = n_rows + 1; end

    if n_rows == 0
        disp('No plots enabled.');
        break;
    end

    tl = tiledlayout(n_rows, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    f_obs  = 1 ./ dobs.T;
    f_pred = 1 ./ dpred.T;

    row = 0;
    legend_added = false;

    % ==============================
    % IMPEDANCE
    % ==============================
    if plot_impedance

        % Real Z
        row = row + 1;
        ax_row = gobjects(1,4);
        for ic = 1:4
            ax = nexttile((row-1)*4+ic);
            hold(ax,'on');

            errorbar(ax, f_obs, ...
                real(dobs.Z(:,comp_idx(ic),is)), ...
                real(dobs.Zerr(:,comp_idx(ic),is)), ...
                'o', 'Color', col_obs, 'MarkerFaceColor', col_obs, ...
                'MarkerSize', 4, 'CapSize', 3, 'LineWidth', 0.75);

            semilogx(ax, f_pred, ...
                real(dpred.Z(:,comp_idx(ic),is)), ...
                '-', 'Color', col_pred, 'LineWidth', 1.8);

            title(ax, comp_names{ic}, 'FontWeight','normal')
            ylabel(ax, 'Re(Z)')
            xlabel(ax, 'Frequency (Hz)')
            set(ax,'XScale','log')
            style_ggplot(ax, panel_bg, grid_col)

            if row == 1 && ic == 4
                legend(ax, 'observed','predicted','Location','northeast','Box','off')
                legend_added = true;
            end

            ax_row(ic) = ax;
        end
        match_ylim(ax_row(1), ax_row(4));  % ZXX / ZYY (diagonal)
        match_ylim(ax_row(2), ax_row(3));  % ZXY / ZYX (off-diagonal)

        % Imaginary Z
        row = row + 1;
        ax_row = gobjects(1,4);
        for ic = 1:4
            ax = nexttile((row-1)*4+ic);
            hold(ax,'on');

            errorbar(ax, f_obs, ...
                imag(dobs.Z(:,comp_idx(ic),is)), ...
                imag(dobs.Zerr(:,comp_idx(ic),is)), ...
                'o', 'Color', col_obs, 'MarkerFaceColor', col_obs, ...
                'MarkerSize', 4, 'CapSize', 3, 'LineWidth', 0.75);

            semilogx(ax, f_pred, ...
                imag(dpred.Z(:,comp_idx(ic),is)), ...
                '-', 'Color', col_pred, 'LineWidth', 1.8);

            title(ax, comp_names{ic}, 'FontWeight','normal')
            ylabel(ax, 'Im(Z)')
            xlabel(ax, 'Frequency (Hz)')
            set(ax,'XScale','log')
            style_ggplot(ax, panel_bg, grid_col)

            ax_row(ic) = ax;
        end
        match_ylim(ax_row(1), ax_row(4));  % ZXX / ZYY (diagonal)
        match_ylim(ax_row(2), ax_row(3));  % ZXY / ZYX (off-diagonal)

    end

    % ==============================
    % APPARENT RESISTIVITY + PHASE
    % ==============================
    if plot_rho_phase

        % Apparent resistivity
        row = row + 1;
        ax_row = gobjects(1,4);
        for ic = 1:4
            ax = nexttile((row-1)*4+ic);
            hold(ax,'on');

            errorbar(ax, f_obs, ...
                dobs.rho(:,comp_idx(ic),is), ...
                dobs.rhoerr(:,comp_idx(ic),is), ...
                'o', 'Color', col_obs, 'MarkerFaceColor', col_obs, ...
                'MarkerSize', 4, 'CapSize', 3, 'LineWidth', 0.75);

            loglog(ax, f_pred, ...
                dpred.rho(:,comp_idx(ic),is), ...
                '-', 'Color', col_pred, 'LineWidth', 1.8);

            title(ax, comp_names{ic}, 'FontWeight','normal')
            ylabel(ax, '\rho_a (\Omega m)')
            xlabel(ax, 'Frequency (Hz)')
            set(ax,'XScale','log','YScale','log')
            style_ggplot(ax, panel_bg, grid_col)

            if row == 1 && ic == 4 && ~legend_added
                legend(ax, 'observed','predicted','Location','northeast','Box','off')
                legend_added = true;
            end

            ax_row(ic) = ax;
        end
        match_ylim(ax_row(1), ax_row(4));  % ZXX / ZYY (diagonal)
        match_ylim(ax_row(2), ax_row(3));  % ZXY / ZYX (off-diagonal)

        % Phase
        row = row + 1;
        ax_row = gobjects(1,4);
        for ic = 1:4
            ax = nexttile((row-1)*4+ic);
            hold(ax,'on');

            errorbar(ax, f_obs, ...
                dobs.pha(:,comp_idx(ic),is), ...
                dobs.phaerr(:,comp_idx(ic),is), ...
                'o', 'Color', col_obs, 'MarkerFaceColor', col_obs, ...
                'MarkerSize', 4, 'CapSize', 3, 'LineWidth', 0.75);

            semilogx(ax, f_pred, ...
                dpred.pha(:,comp_idx(ic),is), ...
                '-', 'Color', col_pred, 'LineWidth', 1.8);

            title(ax, comp_names{ic}, 'FontWeight','normal')
            ylabel(ax, 'Phase (°)')
            xlabel(ax, 'Frequency (Hz)')
            set(ax,'XScale','log')
            style_ggplot(ax, panel_bg, grid_col)

            ax_row(ic) = ax;
        end
        match_ylim(ax_row(1), ax_row(4));  % ZXX / ZYY (diagonal)
        match_ylim(ax_row(2), ax_row(3));  % ZXY / ZYX (off-diagonal)

    end

    % ==============================
    % TIPPER
    % Re(Tx) Im(Tx) Re(Ty) Im(Ty)
    % ==============================
    if plot_tipper

        row = row + 1;
        tip_specs = {
            1, 'real', 'Tx', 'Re(T)';
            1, 'imag', 'Tx', 'Im(T)';
            2, 'real', 'Ty', 'Re(T)';
            2, 'imag', 'Ty', 'Im(T)';
            };

        for k = 1:4
            comp   = tip_specs{k,1};
            fn     = tip_specs{k,2};
            ttl    = tip_specs{k,3};
            ylab   = tip_specs{k,4};

            ax = nexttile((row-1)*4+k);
            hold(ax,'on');

            errorbar(ax, f_obs, ...
                feval(fn, dobs.tip(:,comp,is)), ...
                feval(fn, dobs.tiperr(:,comp,is)), ...
                'o', 'Color', col_obs, 'MarkerFaceColor', col_obs, ...
                'MarkerSize', 4, 'CapSize', 3, 'LineWidth', 0.75);

            semilogx(ax, f_pred, ...
                feval(fn, dpred.tip(:,comp,is)), ...
                '-', 'Color', col_pred, 'LineWidth', 1.8);

            title(ax, ttl, 'FontWeight','normal')
            ylabel(ax, ylab)
            xlabel(ax, 'Frequency (Hz)')
            set(ax,'XScale','log')
            style_ggplot(ax, panel_bg, grid_col)

            if row == 1 && k == 4 && ~legend_added
                legend(ax, 'observed','predicted','Location','northeast','Box','off')
                legend_added = true;
            end
        end

    end

    title(tl, ['Site: ', dobs.site{is}, ...
        ' (', num2str(is), '/', num2str(dobs.ns), ') RMS = ', ...
        num2str(s.rms_site(is), '%.2f')], ...
        'FontSize', 13, 'FontWeight', 'bold')

    choice = menu(['Site: ' dobs.site{is}], ...
        'Next', 'Previous', 'Select from list', 'Quit');

    if choice == 1
        is = min(is+1, dobs.ns);
    elseif choice == 2
        is = max(is-1, 1);
    elseif choice == 3
        is = listdlg('PromptString','Select station:', ...
            'ListString', dobs.site, ...
            'SelectionMode', 'single', ...
            'ListSize', [300 300]);
        if isempty(is)
            is = 1;
        end
    else
        break
    end

end


%% ============================================================
function match_ylim(ax_a, ax_b)
% Forces two axes to share the same y-axis limits, using the union
% of their individual (auto) ranges. Call this AFTER both axes have
% been fully plotted (so autoscaling has already picked sensible
% limits) and BEFORE any manual ylim override you might want later.

    yl_a = ylim(ax_a);
    yl_b = ylim(ax_b);

    yl = [min(yl_a(1), yl_b(1)), max(yl_a(2), yl_b(2))];

    ylim(ax_a, yl);
    ylim(ax_b, yl);
end


%% ============================================================
function style_ggplot(ax, panel_bg, grid_col)
% Applies a ggplot2-inspired theme to a MATLAB axes.
 
    ax.Color     = panel_bg;
    ax.Box       = 'off';
    ax.TickDir   = 'out';
    ax.FontName  = 'Helvetica';
    ax.FontSize  = 12;          % bigger tick labels
    ax.XColor    = [0.35 0.35 0.35];
    ax.YColor    = [0.35 0.35 0.35];
    ax.LineWidth = 1.0;
 
    grid(ax, 'on');
    ax.GridColor      = grid_col;
    ax.GridAlpha      = 1;
    ax.MinorGridColor = grid_col;
    ax.MinorGridAlpha = 0.6;
    ax.GridLineStyle  = '-';
    ax.Layer          = 'bottom';   % gridlines render behind data points
 
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';
 
    % Title styling: this is the component label (ZXX, Re(Z), etc.)
    % keep it bold so it's clearly distinct from the axis labels below
    title(ax, ax.Title.String, 'FontWeight','bold', 'FontSize', 14, ...
        'Color',[0.2 0.2 0.2]);
 
    % Axis labels: bigger, but NOT bold, so the component title stands out
    ax.XLabel.FontSize   = 13;
    ax.XLabel.FontWeight = 'normal';
    ax.XLabel.Color      = [0.15 0.15 0.15];
 
    ax.YLabel.FontSize   = 13;
    ax.YLabel.FontWeight = 'normal';
    ax.YLabel.Color      = [0.15 0.15 0.15];
end


%%
%plot_misfit_rms_map(dobs,dpred,s)
%%
%%
function s_freq = plot_misfit_by_frequency(varargin)
% PLOT_MISFIT_BY_FREQUENCY
% Aggregate misfit (obs vs pred) across ALL sites and components,
% and show how well the fit is at each frequency.
%
% Supports one or more surveys/inversions plotted together for comparison.
%
% Usage:
%
%   plot_misfit_by_frequency(dobs, dpred)
%
%   plot_misfit_by_frequency(dobs, dpred, dobs2, dpred2)
%
%   plot_misfit_by_frequency(..., ...
%       'labels', {'Run A','Run B'})
%
%   plot_misfit_by_frequency(..., 'metric', 'nrms')
%   plot_misfit_by_frequency(..., 'metric', 'abs')
%   plot_misfit_by_frequency(..., 'metric', 'rms')
%
%   plot_misfit_by_frequency(..., 'quantity', 'rho_a')
%
%
% -------------------------------------------------------------------------
% METRICS
% -------------------------------------------------------------------------
%
% 'abs'
%     Mean absolute raw residual:
%
%         mean(abs(obs - pred))
%
%
% 'nrms'
%     Mean absolute normalized residual:
%
%         mean(abs((obs - pred) ./ error))
%
%     This is useful for showing the average deviation in units of the
%     assigned data error.
%
%
% 'rms'
%     Root-mean-square normalized residual:
%
%         sqrt(mean(((obs - pred) ./ error).^2))
%
%     This is the metric that is directly comparable to the usual
%     normalized RMS reported by an inversion.
%
%
% -------------------------------------------------------------------------
% QUANTITIES
% -------------------------------------------------------------------------
%
% 'Z'
%     Misfit calculated separately for the real and imaginary parts of
%     the impedance tensor.
%
% 'rho_a'
%     Misfit calculated on apparent resistivity derived from |Z|:
%
%         rho_a = rho_const * T * |Z|^2
%
%     Phase is ignored.
%
%
% -------------------------------------------------------------------------
% OPTIONS
% -------------------------------------------------------------------------
%
% 'labels'
%     Cell array of names, one per (dobs,dpred) pair.
%
%     Default:
%         {'Survey 1','Survey 2',...}
%
% 'metric'
%     'nrms' (default), 'abs', or 'rms'
%
% 'quantity'
%     'Z' (default) or 'rho_a'
%
% 'rho_const'
%     Scalar prefactor used in
%
%         rho_a = rho_const * T * |Z|^2
%
%     Default: 0.2
%
% 'include_tipper'
%     true/false (default true).
%
%     Ignored for quantity = 'rho_a'.
%
% 'freq_tol'
%     Relative tolerance in log10 frequency used to group frequencies
%     that are nominally identical but differ slightly due to floating
%     point precision.
%
%     Default: 1e-6
%
% 'save_path'
%     If given, exports the figure to this path.
%
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
%
% Returns s_freq, a struct array with fields:
%
%     label
%     Z
%     T
%
% Each quantity contains:
%
%     freq
%     n
%     mean_abs
%     rms
%     std_abs
%     std_res
%
% For 'abs' and 'nrms':
%
%     mean_abs = mean(abs(residual))
%     std_abs  = std(abs(residual))
%
% For 'rms':
%
%     rms      = sqrt(mean(residual.^2))
%     std_res  = std(residual)
%
% where the residual for RMS is the NORMALIZED residual:
%
%     residual = (obs - pred) ./ error
%
% -------------------------------------------------------------------------


%% ========================================================================
% 1) Split inputs into dataset pairs vs Name-Value options
% ========================================================================

dsets = {};
i = 1;

while i <= numel(varargin) && isstruct(varargin{i})

    if i+1 > numel(varargin) || ~isstruct(varargin{i+1})
        error('plot_misfit_by_frequency:pairing', ...
            'Each dobs must be followed by a matching dpred struct.');
    end

    dsets{end+1} = struct( ... %#ok<AGROW>
        'dobs', varargin{i}, ...
        'dpred', varargin{i+1});

    i = i + 2;
end

opts_in = varargin(i:end);

if isempty(dsets)
    error('plot_misfit_by_frequency:noData', ...
        'Provide at least one (dobs, dpred) pair.');
end

nsurv = numel(dsets);


%% ========================================================================
% 2) Input parameters
% ========================================================================

p = inputParser;

addParameter(p, 'metric', 'nrms');
addParameter(p, 'quantity', 'Z');
addParameter(p, 'rho_const', 0.2);
addParameter(p, 'include_tipper', true);
addParameter(p, 'freq_tol', 1e-6);
addParameter(p, 'labels', {});
addParameter(p, 'save_path', '');

parse(p, opts_in{:});

metric         = lower(p.Results.metric);
quantity       = p.Results.quantity;
rho_const      = p.Results.rho_const;
include_tipper = p.Results.include_tipper;
freq_tol       = p.Results.freq_tol;
labels         = p.Results.labels;
save_path      = p.Results.save_path;


%% ========================================================================
% 3) Validate input parameters
% ========================================================================

if ~ismember(lower(quantity), {'z', 'rho_a'})
    error('plot_misfit_by_frequency:quantity', ...
        'Unknown quantity "%s". Use ''Z'' or ''rho_a''.', quantity);
end

if ~ismember(metric, {'nrms', 'abs', 'rms'})
    error('plot_misfit_by_frequency:metric', ...
        'Unknown metric "%s". Use ''nrms'', ''abs'', or ''rms''.', metric);
end

is_rho_a = strcmpi(quantity, 'rho_a');

if is_rho_a && include_tipper

    warning('plot_misfit_by_frequency:tipperIgnored', ...
        ['Tipper has no apparent-resistivity analog; include_tipper is ' ...
         'being forced to false because quantity=''rho_a''.']);

    include_tipper = false;
end


%% ========================================================================
% 4) Labels
% ========================================================================

if isempty(labels)

    labels = arrayfun(@(k) sprintf('Survey %d', k), ...
        1:nsurv, ...
        'UniformOutput', false);

elseif numel(labels) ~= nsurv

    error('plot_misfit_by_frequency:labels', ...
        ['Number of labels (%d) must match number of ' ...
         '(dobs,dpred) pairs (%d).'], ...
        numel(labels), ...
        nsurv);

end


%% ========================================================================
% 5) ggplot2-style panel styling
% ========================================================================

panel_bg = [0.92 0.92 0.95];
grid_col = [1.00 1.00 1.00];


% One hue per survey

base_palette = [ ...
    0.12 0.47 0.71;   % blue
    0.89 0.10 0.11;   % red
    0.20 0.63 0.17;   % green
    1.00 0.50 0.00;   % orange
    0.42 0.24 0.60;   % purple
    0.65 0.34 0.16];  % brown


if nsurv > size(base_palette,1)

    base_palette = [ ...
        base_palette; ...
        lines(nsurv - size(base_palette,1))];

end

surv_cols = base_palette(1:nsurv, :);


%% ========================================================================
% 6) Compute per-frequency statistics for each survey
% ========================================================================

s_freq = struct('label', {}, 'Z', {}, 'T', {});


for k = 1:nsurv

    dobs  = dsets{k}.dobs;
    dpred = dsets{k}.dpred;


    % ---------------------------------------------------------------------
    % Frequency
    % ---------------------------------------------------------------------

    f = 1 ./ dobs.T(:);

    ncomp = size(dobs.Z, 2);
    nsite = size(dobs.Z, 3);

    freqArr = repmat(f, [1, ncomp, nsite]);


    % ---------------------------------------------------------------------
    % Apparent-resistivity misfit
    % ---------------------------------------------------------------------

    if is_rho_a

        Tper = repmat(dobs.T(:), [1, ncomp, nsite]);


        % Apparent resistivity

        rho_obs = ...
            rho_const .* Tper .* abs(dobs.Z).^2;

        rho_pred = ...
            rho_const .* Tper .* abs(dpred.Z).^2;


        % ---------------------------------------------------------------
        % Propagate Z error onto rho_a
        %
        % rho_a = rho_const * T * |Z|^2
        %
        % First-order error propagation:
        %
        % sigma_rho_a ~= 2 * rho_const * T * |Z| * sigma_Z
        % ---------------------------------------------------------------

        rho_err = ...
            2 .* rho_const .* Tper .* ...
            abs(dobs.Z) .* dobs.Zerr;


        switch metric

            case 'abs'

                % Raw absolute residual.
                %
                % The absolute value itself is applied later during
                % frequency binning.

                res = rho_obs - rho_pred;


            case 'nrms'

                % Normalized residual

                res = ...
                    (rho_obs - rho_pred) ./ rho_err;

                res(rho_err <= 0) = NaN;


            case 'rms'

                % IMPORTANT:
                %
                % RMS is now calculated from the NORMALIZED residual.
                % This makes the frequency-dependent RMS comparable to
                % the overall normalized RMS of an inversion.

                res = ...
                    (rho_obs - rho_pred) ./ rho_err;

                res(rho_err <= 0) = NaN;

        end


        res_Z  = res(:);
        freq_Z = freqArr(:);


    % ---------------------------------------------------------------------
    % Impedance misfit
    % ---------------------------------------------------------------------

    else

        % Real part

        [res_Zr, freq_Zr] = local_residuals( ...
            real(dobs.Z), ...
            real(dpred.Z), ...
            real(dobs.Zerr), ...
            freqArr, ...
            metric);


        % Imaginary part

        [res_Zi, freq_Zi] = local_residuals( ...
            imag(dobs.Z), ...
            imag(dpred.Z), ...
            imag(dobs.Zerr), ...
            freqArr, ...
            metric);


        res_Z = [ ...
            res_Zr(:); ...
            res_Zi(:)];


        freq_Z = [ ...
            freq_Zr(:); ...
            freq_Zi(:)];

    end


    % ---------------------------------------------------------------------
    % Bin impedance / apparent-resistivity residuals by frequency
    % ---------------------------------------------------------------------

    s_freq(k).label = labels{k};

    s_freq(k).Z = local_bin_by_frequency( ...
        freq_Z, ...
        res_Z, ...
        freq_tol);


    % ---------------------------------------------------------------------
    % Tipper
    % ---------------------------------------------------------------------

    have_tipper = ...
        include_tipper && ...
        isfield(dobs, 'tip');


    if have_tipper

        ncomp_t = size(dobs.tip, 2);

        freqArr_t = repmat( ...
            f, ...
            [1, ncomp_t, nsite]);


        % Real part

        [res_Tr, freq_Tr] = local_residuals( ...
            real(dobs.tip), ...
            real(dpred.tip), ...
            real(dobs.tiperr), ...
            freqArr_t, ...
            metric);


        % Imaginary part

        [res_Ti, freq_Ti] = local_residuals( ...
            imag(dobs.tip), ...
            imag(dpred.tip), ...
            imag(dobs.tiperr), ...
            freqArr_t, ...
            metric);


        res_T = [ ...
            res_Tr(:); ...
            res_Ti(:)];


        freq_T = [ ...
            freq_Tr(:); ...
            freq_Ti(:)];


        s_freq(k).T = local_bin_by_frequency( ...
            freq_T, ...
            res_T, ...
            freq_tol);

    end

end


%% ========================================================================
% 7) Plot
% ========================================================================

have_tipper = ...
    include_tipper && ...
    isfield(dsets{1}.dobs, 'tip');


fig = figure(20);
clf;

fig.Color = 'w';

fig.Units = 'centimeters';

fig.Position = [2 2 18 12];


ax1 = axes(fig);

hold(ax1, 'on');


quant_name = 'Z';

if is_rho_a
    quant_name = '\rho_a';
end


for k = 1:nsurv

    col_dark = surv_cols(k,:);

    col_light = ...
        col_dark + ...
        (1 - col_dark) * 0.55;


    % ---------------------------------------------------------------------
    % Impedance / apparent resistivity
    % ---------------------------------------------------------------------

    plot_one_metric( ...
        ax1, ...
        s_freq(k).Z, ...
        col_dark, ...
        '-', ...
        sprintf('%s \\cdot %s', ...
            s_freq(k).label, ...
            quant_name), ...
        metric);


    % ---------------------------------------------------------------------
    % Tipper
    % ---------------------------------------------------------------------

    if have_tipper

        plot_one_metric( ...
            ax1, ...
            s_freq(k).T, ...
            col_light, ...
            '--', ...
            sprintf('%s \\cdot T', ...
                s_freq(k).label), ...
            metric);

    end

end


%% ========================================================================
% Axis labels
% ========================================================================

ylabel(ax1, ...
    local_metric_label(metric, quantity), ...
    'Interpreter', 'tex');


xlabel(ax1, 'Frequency (Hz)');

set(ax1, 'XScale', 'log');


%% ========================================================================
% Target line
% ========================================================================

% A normalized RMS of 1 corresponds approximately to residuals that are
% on the scale of the assigned data errors.
%
% Only add this line for RMS because RMS is now a normalized quantity.

if strcmp(metric, 'nrms')

    yline(ax1, ...
        1, ...
        '--', ...
        'target = 1', ...
        'Color', [0.35 0.35 0.35], ...
        'LineWidth', 1, ...
        'FontSize', 9, ...
        'LabelHorizontalAlignment', 'left', ...
        'HandleVisibility', 'off');

end


%% ========================================================================
% Legend
% ========================================================================

lg = legend(ax1, ...
    'Location', 'northoutside', ...
    'Orientation', 'horizontal', ...
    'Box', 'off', ...
    'NumColumns', min(nsurv*2, 4), ...
    'FontSize', 14);


lg.ItemTokenSize = [18 18];


%% ========================================================================
% Style
% ========================================================================

style_ggplot_local( ...
    ax1, ...
    panel_bg, ...
    grid_col);


ax1.LineWidth = 0.9;


%% ========================================================================
% Save
% ========================================================================

if ~isempty(save_path)

    exportgraphics( ...
        fig, ...
        save_path, ...
        'Resolution', 300);

end


end % main function



%% =========================================================================
function [res, freq] = local_residuals(obs, pred, err, freqArr, metric)

% LOCAL_RESIDUALS
%
% Calculates the residuals used for the requested metric.
%
% 'abs'
%     Raw residual:
%
%         obs - pred
%
% 'nrms'
%     Normalized residual:
%
%         (obs - pred) ./ err
%
% 'rms'
%     Normalized residual:
%
%         (obs - pred) ./ err
%
% The distinction between nrms and rms is made later during the frequency
% aggregation:
%
%     nrms -> mean(abs(residual))
%
%     rms  -> sqrt(mean(residual.^2))


switch metric

    case 'abs'

        % Raw residual

        res = obs - pred;


    case {'nrms', 'rms'}

        % Normalized residual

        res = ...
            (obs - pred) ./ err;

        % Invalid / zero errors

        res(err <= 0) = NaN;


    otherwise

        error( ...
            'Unknown metric "%s". Use ''nrms'', ''abs'', or ''rms''.', ...
            metric);

end


freq = freqArr;

end



%% =========================================================================
function s = local_bin_by_frequency(freq, res, freq_tol)

% LOCAL_BIN_BY_FREQUENCY
%
% Groups (frequency,residual) pairs by approximately identical frequency
% and calculates the statistical quantities used for plotting.
%
%
% For 'abs' and 'nrms':
%
%     mean_abs = mean(abs(res))
%     std_abs  = std(abs(res))
%
%
% For 'rms':
%
%     rms      = sqrt(mean(res.^2))
%     std_res  = std(res)
%
%
% The function calculates all quantities because the metric is selected
% later by plot_one_metric().


%% ------------------------------------------------------------------------
% Remove invalid values
% -------------------------------------------------------------------------

valid = ...
    ~isnan(freq) & ...
    ~isnan(res) & ...
    isfinite(freq) & ...
    isfinite(res);


freq = freq(valid);

res = res(valid);


%% ------------------------------------------------------------------------
% Group frequencies in log space
% -------------------------------------------------------------------------

logf = log10(freq);


[ufreq, ~, ic] = uniquetol( ...
    logf, ...
    freq_tol);


ufreq = 10.^ufreq;

nbins = numel(ufreq);


%% ------------------------------------------------------------------------
% Preallocate
% -------------------------------------------------------------------------

n = zeros(nbins,1);

mean_abs = zeros(nbins,1);

rms_val = zeros(nbins,1);

std_abs = zeros(nbins,1);

std_res = zeros(nbins,1);


%% ------------------------------------------------------------------------
% Calculate statistics for each frequency
% -------------------------------------------------------------------------

for kk = 1:nbins

    r = res(ic == kk);


    n(kk) = numel(r);


    % Mean absolute residual

    mean_abs(kk) = ...
        mean(abs(r));


    % RMS residual

    rms_val(kk) = ...
        sqrt(mean(r.^2));


    % Standard deviation of absolute residuals

    std_abs(kk) = ...
        std(abs(r));


    % Standard deviation of signed residuals

    std_res(kk) = ...
        std(r);

end


%% ------------------------------------------------------------------------
% Sort from high frequency to low frequency
% -------------------------------------------------------------------------

[ufreq, sortIdx] = ...
    sort(ufreq, 'descend');


s.freq = ufreq;

s.n = n(sortIdx);

s.mean_abs = mean_abs(sortIdx);

s.rms = rms_val(sortIdx);

s.std_abs = std_abs(sortIdx);

s.std_res = std_res(sortIdx);

end



%% =========================================================================
function plot_one_metric(ax, s, col, lstyle, name, metric)

% PLOT_ONE_METRIC
%
% Selects the quantity to plot and the corresponding standard deviation.
%
% For 'abs':
%
%     mean(abs(raw residual))
%
% For 'nrms':
%
%     mean(abs(normalized residual))
%
% For 'rms':
%
%     sqrt(mean(normalized residual.^2))


switch lower(metric)

    case {'abs', 'nrms'}

        % Mean absolute residual

        y = s.mean_abs;


        % Standard deviation of absolute residuals

        sigma = s.std_abs;


    case 'rms'

        % RMS of NORMALIZED residuals

        y = s.rms;


        % Standard deviation of signed normalized residuals

        sigma = s.std_res;


    otherwise

        error( ...
            'Unknown metric "%s".', ...
            metric);

end


%% ------------------------------------------------------------------------
% Shaded +/- one standard deviation
% -------------------------------------------------------------------------

lo = max(y - sigma, 0);

hi = y + sigma;


fill(ax, ...
    [s.freq; flipud(s.freq)], ...
    [lo; flipud(hi)], ...
    col, ...
    'FaceAlpha', 0.22, ...
    'EdgeColor', col, ...
    'EdgeAlpha', 0.35, ...
    'LineWidth', 0.5, ...
    'HandleVisibility', 'off');


%% ------------------------------------------------------------------------
% Central line
% -------------------------------------------------------------------------

plot(ax, ...
    s.freq, ...
    y, ...
    ['o' lstyle], ...
    'Color', col, ...
    'MarkerFaceColor', col, ...
    'MarkerEdgeColor', 'w', ...
    'MarkerSize', 5, ...
    'LineWidth', 2, ...
    'DisplayName', name);

end



%% =========================================================================
function lbl = local_metric_label(metric, quantity)

is_rho_a = strcmpi(quantity, 'rho_a');


switch lower(metric)

    case 'nrms'

        if is_rho_a

            lbl = 'Normalized residual ';
            ylabel(lbl, 'Interpreter', 'latex');

        else

            lbl = 'Normalized residual ';
            ylabel(lbl, 'Interpreter', 'latex');

        end


    case 'abs'

        if is_rho_a

            lbl = ...
                'Mean |\rho_{a,obs} - \rho_{a,pred}| (\Omega\cdotm)';

        else

            lbl = ...
                'Mean |obs - pred|';

        end


    case 'rms'

        if is_rho_a

            lbl = ...
                'RMS normalized \rho_a residual (\sigma units)';

        else

            lbl = ...
                'RMS normalized residual (\sigma units)';

        end


    otherwise

        lbl = 'Misfit';

end

end



%% =========================================================================
function style_ggplot_local(ax, panel_bg, grid_col)

ax.Color = panel_bg;

ax.Box = 'off';

ax.TickDir = 'out';

ax.FontName = 'Helvetica';

ax.FontSize = 12;


ax.XColor = [0.35 0.35 0.35];

ax.YColor = [0.35 0.35 0.35];

ax.LineWidth = 0.75;


grid(ax, 'on');


ax.GridColor = grid_col;

ax.GridAlpha = 1;


ax.MinorGridColor = grid_col;

ax.MinorGridAlpha = 0.6;


ax.GridLineStyle = '-';

ax.Layer = 'bottom';


ax.XMinorGrid = 'on';

ax.YMinorGrid = 'on';

end



%% =========================================================================
% Example data loading and plotting
% =========================================================================

% Uncomment and adapt these lines if the data are not already loaded.

dobs2 = load_data_modem( ...
    '\\wsl.localhost\Ubuntu\home\ramonvanoli\MasterThesis\Inversion\Dike1\2mgrid\28ohmm\errfl_10_5\FINAL\Datafile.data');

dpred2 = load_data_modem( ...
    "\\wsl.localhost\Ubuntu\home\ramonvanoli\MasterThesis\Inversion\Dike1\2mgrid\28ohmm\errfl_10_5\FINAL\DIKE1_NLCG_021.dat");


% -------------------------------------------------------------------------
% Choose metric:
%
% 'abs'  = mean absolute raw residual
%
% 'nrms' = mean absolute normalized residual
%
% 'rms'  = RMS normalized residual
%
% -------------------------------------------------------------------------

plot_misfit_by_frequency( ...
    dobs, dpred, ...
    dobs2, dpred2, ...
    'labels', {'DIKE 3', 'DIKE 1'}, ...
    'include_tipper', true, ...
    'metric', 'nrms');


% -------------------------------------------------------------------------
% Export figure
% -------------------------------------------------------------------------
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, ...
    'C:\Users\Baar\OneDrive - Vanoli AG\HS_25\Master Thesis\Figures\Inversion\misfit_by_freq_nrms.png', ...
    'Resolution', 600);
