%% 
%last updated 3/1/26
% Multi-file processing with per-file full plotting and clean Excel output
clear

% === File Selection ===
[fileNames, pathName] = uigetfile('*.csv', 'Select CSV File(s)', 'MultiSelect', 'on');
if isequal(fileNames, 0)
    disp('User canceled file selection.');
    return;
end
if ischar(fileNames)
    fileNames = {fileNames};
end

numFiles = numel(fileNames);

% === Initialize Storage (outside loop) ===
P = nan(numFiles,1);                % average power per file

% Mode / interval stats (clean, preallocated column vectors)
pct_modeA_all     = nan(numFiles,1);
pct_modeB_all     = nan(numFiles,1);
mean_ms_all       = nan(numFiles,1);
std_ms_all        = nan(numFiles,1);
median_ms_all     = nan(numFiles,1);
sample_count_all  = nan(numFiles,1);
sample_count_all_green = nan(numFiles,1);
square_count_all   = nan(numFiles,1);
triangle_count_all = nan(numFiles,1);

% Peak voltage statistics (per file)
peakV_mean_all   = nan(numFiles,1);
peakV_count_all  = nan(numFiles,1);

peakV_std_all    = nan(numFiles,1);
peakV_median_all = nan(numFiles,1);

avg_maxV_all = nan(numFiles,1);
std_maxV_all = nan(numFiles,1);

avg_minV_all = nan(numFiles,1);
std_minV_all = nan(numFiles,1);

overall_meanV_all = nan(numFiles,1);
overall_stdV_all  = nan(numFiles,1);

% figure handles (stable numbering)
fig_wave = 3;        % voltage/current figure
fig_summary = 20;    % per-file stacked bar + histogram

header_row = 14;     % NumHeaderLines for readmatrix

%% === MAIN LOOP: process each file and produce all plots ===
for i = 1:numFiles
    fprintf('\nProcessing file %d of %d: %s\n', i, numFiles, fileNames{i});
    fullpath = fullfile(pathName, fileNames{i});
    data = readmatrix(fullpath, 'NumHeaderLines', header_row);

    % Column mapping (unchanged)
    T_idx = 1;   % TIME is first column
    LowI_idx = 2;   % CH2
    HighI_idx = 3;  % CH3
    V_idx = 6;      % MATH1

    % Per-file assignments (unchanged)
    T = data(:, T_idx);
    lowI = data(:, LowI_idx);
    highI = data(:, HighI_idx);
    V = data(:, V_idx);

    T = fill_invalid_with_nearest(T);
    V = fill_invalid_with_nearest(V);
    lowI = fill_invalid_with_nearest(lowI);
    highI = fill_invalid_with_nearest(highI);

    % Current Fixing (unchanged)
    max_signal = max(lowI);
    highI(abs(highI) < max_signal) = 0;
    I = highI + lowI;

   
    % Sampling frequency
    fs = 1 / mean(diff(T)); % Hz

    % Subplot 1: Voltage vs Time
    subplot(2,1,1);
    plot(T, V, 'k', 'LineWidth', 1.1);
    xlabel('Time [s]', 'FontSize', 12);
    ylabel('Voltage [V]', 'FontSize', 12);
    title('Voltage vs Time', 'FontSize', 13);
    grid on;

    % Subplot 2: STFT Intensity vs Time
    subplot(2,1,2);

    window = 512;
    noverlap = round(window * 0.75);
    nfft = 2048;
    [~, ~, t_stft, ps] = spectrogram(V, window, noverlap, nfft, fs, 'yaxis');

    intensity = mean(abs(ps).^2, 1);
    % === MOVING AVERAGE FOR STABILITY (NEW) ==

    ma_window = 5;  % number of STFT points to smooth over (adjust if needed)
    intensity_ma = movmean(intensity, ma_window);

    [peak_vals, peak_locs] = findpeaks(intensity_ma);
    peak_times = t_stft(peak_locs);

    %[peak_vals, peak_locs] = findpeaks(intensity);
    %peak_times = t_stft(peak_locs);

    % Parameters
    pct_within = 0.3;
    within_multiplier = 1 - pct_within;
    drop_pct = pct_within;
    min_required_peaks = 7;
    
    %min_peak_level = 15;
    drop_check_length = 3;      % Number of consecutive peaks required for drop validation
    Min_AmperageSpike = 0.175;
    TimeTolerance = 0.000125; %(100 microseconds)
    AmperageSearchWindow = 0.0002750;   % ± window (microseconds) for counting current peaks
    Mode_AmperageLimit = 1.75*Min_AmperageSpike;
    CurrentSpacingVariable = 1e-7;   % 



    [min_peak_values, min_peak_locs] = findpeaks(-intensity_ma);

    % Convert back to actual minima values
    min_vals = -min_peak_values;
    
    % Mean level of the minima
    min_peak_level = mean(min_vals) * 0.30;
    
    disp(min_peak_level)


    plot(t_stft, intensity, 'r', ...
    'LineWidth', 1.2, ...
    'Color', [1 0 0 0.25]); hold on;   % raw (transparent)

    xlabel('Time [s]'); ylabel('Relative Intensity (a.u.)');
    grid off;
    
    current_run_idx = [];
    prev_peak = [];
    green_marker_times = [];
    purple_marker_times = [];
    purple_marker_amps  = [];
    drop_check_start_idx = 0;   % Index where the drop checking sequence started
    
    
    for k = 1:numel(peak_vals)
    pk = peak_vals(k);
    pk_time = peak_times(k);
    
    if isempty(current_run_idx)
        % Start a new run
        current_run_idx = k;
        prev_peak = pk;
        drop_check_start_idx = 0; % Reset drop check start
        continue;
    end
    
    % --- RUN CONTINUATION CRITERIA (Original Logic) ---
    % Continues if the peak increases (pk >= prev_peak) & Above 50 OR stays within the multiplier
    is_valid_peak = pk >= min_peak_level;

    continues = is_valid_peak && ...
           ((pk >= prev_peak) || (pk >= prev_peak * within_multiplier));
    
    if continues
        % A run-continuation criteria was met, so any prior drop-check is invalid
        drop_check_start_idx = 0; 
        
        current_run_idx(end+1) = k;
        prev_peak = pk;
        continue;
    end
    
    % --- PERSISTENT DROP VALIDATION (New Logic) ---
    
    % If we reach here, 'continues' is false. This peak (pk) failed to continue the run.
    
    % The critical peak to reference for the drop is the *last peak successfully added to the run*.
    % Use peak_vals(current_run_idx(end)) as the reference for drop_pct
    
    % Safety check: If current_run_idx is empty, something is wrong, but typically it holds at least one peak.
    if isempty(current_run_idx)
        prev_run_peak_val = 0; % Should not happen based on loop structure
    else
        prev_run_peak_val = peak_vals(current_run_idx(end));
    end
    
    % Check if the current peak (pk) is sufficiently low compared to the *last run peak*
    is_low_enough = (pk < prev_run_peak_val * drop_pct);

    if is_low_enough
        if drop_check_start_idx == 0
            % This is the first peak in the potential drop sequence
            drop_check_start_idx = k;
            
        elseif (k - drop_check_start_idx + 1) >= drop_check_length
            % This is the N-th (5th) consecutive peak below the drop threshold! Run END.

            % 1. Evaluate and record the SUCCESSFUL run (up to the peak BEFORE the drop sequence)
            successful_run_indices = current_run_idx;
            
            if numel(successful_run_indices) >= min_required_peaks
                run_times = peak_times(successful_run_indices);
                end_idx = successful_run_indices(end); % Last peak BEFORE the drop sequence
                
                % Record the green window
                extended_end_idx = min(end_idx + 2, numel(peak_times));
                extended_end_time = peak_times(extended_end_idx);
                
                % Plot marker at the end of the successful run
                plot(peak_times(end_idx), 0, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
                green_marker_times(end+1,1) = peak_times(end_idx);
            end
            
            % 2. Restart a new run from the CURRENT peak (pk)
            current_run_idx = k;
            prev_peak = pk;
            drop_check_start_idx = 0;
            continue;
        end
    else
        % The current peak (pk) failed continuation, AND it was NOT low enough 
        % to continue the drop check sequence.
        
        % This single peak breaks the persistent drop sequence.
        drop_check_start_idx = 0;
        
        % The previous run ENDED at current_run_idx(end).
        % A new run MUST START at the current peak (k).
        
        % Check the run just ended for validity
        if numel(current_run_idx) >= min_required_peaks
            run_times = peak_times(current_run_idx);
            end_idx = current_run_idx(end);
            
            
            % Plot marker at the end of the successful run
            plot(peak_times(end_idx), 0, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
            green_marker_times(end+1,1) = peak_times(end_idx);
        end

        % Start a new run from the current peak (pk)
        current_run_idx = k;
        prev_peak = pk;
        continue;
    end
    
end

% --- EVALUATE LAST RUN ---
if ~isempty(current_run_idx) && numel(current_run_idx) >= min_required_peaks && drop_check_start_idx == 0
    % Only evaluate the last run if it didn't end prematurely and is not currently stuck in a drop check
    run_times = peak_times(current_run_idx);
    end_idx = current_run_idx(end);
    extended_end_idx = min(end_idx + 2, numel(peak_times));
    extended_end_time = peak_times(extended_end_idx);
    last_pk_time = peak_times(current_run_idx(end));
    plot(last_pk_time, 0, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
    green_marker_times(end+1,1) = last_pk_time;
end


% ============================================================
% BLUE MARKERS  CURRENT-VALIDATED EVENTS
% ============================================================

blue_marker_times = [];   % SEPARATE storage (do NOT mix with green)

for g = 1:numel(green_marker_times)
    t0 = green_marker_times(g);

    % Time tolerance window
    t_low  = t0 - TimeTolerance;
    t_high = t0 + TimeTolerance;

    % Indices in the tolerance window
    idx = find(T >= t_low & T <= t_high);

    if isempty(idx)
        continue;
    end

    % Find samples exceeding the amperage threshold
    valid_idx = idx(abs(I(idx)) >= Min_AmperageSpike);

    if isempty(valid_idx)
        continue;
    end

    % Choose the spike closest in time to the green marker
    [~, k] = min(abs(T(valid_idx) - t0));
    best_idx = valid_idx(k);

    % Store the aligned time
    blue_t = T(best_idx);
    blue_marker_times(end+1,1) = blue_t;

    % Plot BLUE marker at the ACTUAL current spike time
    plot(blue_t, 0, 'bo', ...
        'MarkerFaceColor','b', ...
        'MarkerSize', 6);
end

% ============================================================
% PEAK VOLTAGE ANALYSIS AROUND CONFIRMED (BLUE) MARKERS
% ============================================================

peak_voltage_vals = nan(numel(blue_marker_times),1);

for b = 1:numel(blue_marker_times)

    t0 = blue_marker_times(b);

    % Define search window ± TimeTolerance
    t_low  = t0 - TimeTolerance;
    t_high = t0 + TimeTolerance;

    idx = find(T >= t_low & T <= t_high);

    if isempty(idx)
        continue;
    end

    % Maximum absolute voltage in window
    peak_voltage_vals(b) = max(abs(V(idx)));

end

% Remove any NaNs (in case of edge effects)
peak_voltage_vals = peak_voltage_vals(~isnan(peak_voltage_vals));

% ============================================================
% VOLTAGE STATISTICS
% ============================================================

if ~isempty(peak_voltage_vals)

    % Absolute voltage statistics
    peakV_mean_all(i)   = mean(peak_voltage_vals);
    peakV_std_all(i)    = std(peak_voltage_vals);
    peakV_median_all(i) = median(peak_voltage_vals);
    peakV_count_all(i)  = numel(peak_voltage_vals);

    % Signed voltage peaks (true polarity)
    signed_peaks = nan(numel(blue_marker_times),1);

    for b = 1:numel(blue_marker_times)
        t0 = blue_marker_times(b);

        t_low  = t0 - TimeTolerance;
        t_high = t0 + TimeTolerance;

        idx = find(T >= t_low & T <= t_high);

        if isempty(idx)
            continue;
        end

        % signed peak voltage
        [~,k] = max(abs(V(idx)));
        signed_peaks(b) = V(idx(k));
    end

    signed_peaks = signed_peaks(~isnan(signed_peaks));

    % Positive peaks (maxima)
    max_peaks = signed_peaks(signed_peaks > 0);

    % Negative peaks (minima)
    min_peaks = signed_peaks(signed_peaks < 0);

    % Average max voltage
    avg_maxV_all(i) = mean(max_peaks);
    std_maxV_all(i) = std(max_peaks);

    % Average min voltage
    avg_minV_all(i) = mean(min_peaks);
    std_minV_all(i) = std(min_peaks);

    % Overall average (signed)
    overall_meanV_all(i) = mean(signed_peaks);
    overall_stdV_all(i)  = std(signed_peaks);

else

    peakV_mean_all(i)   = NaN;
    peakV_std_all(i)    = NaN;
    peakV_median_all(i) = NaN;
    peakV_count_all(i)  = 0;

end

% ============================================================
% HISTOGRAM OF PEAK VOLTAGES
% ============================================================

fig_peakV = 150;   % choose unused figure number
figure(fig_peakV); clf(fig_peakV);
set(fig_peakV,'Position',[300 250 800 500]);

if ~isempty(peak_voltage_vals)

    histogram(peak_voltage_vals,20);
    xlabel('Peak |Voltage| within \pm TimeTolerance (V)');
    ylabel('Count');
    title('Distribution of Confirmed-Peak Voltages');
    grid off;

    % ---- Statistics ----
    Np = numel(peak_voltage_vals);
    mean_peakV = mean(peak_voltage_vals);
    std_peakV  = std(peak_voltage_vals);
    median_peakV = median(peak_voltage_vals);

    txt = sprintf(['N = %d samples\n' ...
                   'Mean:   %.2f V\n' ...
                   'Std:    %.2f V\n' ...
                   'Median: %.2f V'], ...
                   Np, mean_peakV, std_peakV, median_peakV);

    xlim_vals = xlim;
    ylim_vals = ylim;

    text(xlim_vals(1) + 0.65*(xlim_vals(2)-xlim_vals(1)), ...
         ylim_vals(1) + 0.80*(ylim_vals(2)-ylim_vals(1)), ...
         txt, ...
         'FontSize',12, ...
         'BackgroundColor',[1 1 1 .7], ...
         'EdgeColor',[0 0 0], ...
         'Margin',5);

else
    text(0.5,0.5,'No confirmed peaks available for voltage analysis', ...
        'HorizontalAlignment','center','FontSize',12);
    axis off;
end










% ============================================================
% CYAN MARKER CLASSIFICATION (CURRENT PEAK COUNT)
% ============================================================

square_marker_times   = [];
triangle_marker_times = [];

for b = 1:numel(blue_marker_times)
    t0 = blue_marker_times(b);

    % Define amperage search window
    t_low  = t0;
    t_high = t0 + AmperageSearchWindow;

    idx = (T >= t_low & T <= t_high);

    if ~any(idx)
        continue;
    end

    % Use absolute current
    I_seg = abs(I(idx));

    % --- Find raw current peaks ---
[amp_pks, locs] = findpeaks(I_seg);
peak_times_I = T(idx);
peak_times_I = peak_times_I(locs);

I_window = I(idx);                 % signed current segment
signed_amp_pks = I_window(locs);   % signed peak amplitudes

% --- Cluster peaks within 1 microsecond ---
keep_mask = false(size(amp_pks));
used = false(size(amp_pks));

    for p = 1:numel(amp_pks)
        if used(p), continue; end
    
        % Find peaks within spacing window
        dt = abs(peak_times_I - peak_times_I(p));
        cluster_idx = find(dt <= CurrentSpacingVariable);
    
        % Select largest magnitude peak in cluster
        [~, imax] = max(amp_pks(cluster_idx));
        keep_idx = cluster_idx(imax);
    
        keep_mask(keep_idx) = true;
        used(cluster_idx) = true;
    end
    
    % --- Retained peaks only ---
    amp_pks = amp_pks(keep_mask);
    peak_times_I = peak_times_I(keep_mask);
    
    % --- Count peaks above Mode_AmperageLimit ---
    valid_peaks = amp_pks(amp_pks >= Mode_AmperageLimit);
    n_peaks = numel(valid_peaks);
    
    strong_idx = amp_pks >= Mode_AmperageLimit;
    purple_marker_times = [purple_marker_times; peak_times_I(strong_idx)];
    purple_marker_amps  = [purple_marker_amps;  signed_amp_pks(strong_idx)];


    % Classification
    if n_peaks >= 0 && n_peaks <= 4 %0 covers for case marker gets placed after a single current spike
        square_marker_times(end+1,1) = t0;
    elseif n_peaks >= 5
        triangle_marker_times(end+1,1) = t0;
    end
end




green_count = numel(green_marker_times);
blue_count  = numel(blue_marker_times);
% If no peaks, skip overlays & mark NaNs later
if isempty(peak_times)
    hold off;
    pct_modeA_all(i) = NaN;
    pct_modeB_all(i) = NaN;
    mean_ms_all(i) = NaN;
    std_ms_all(i) = NaN;
    median_ms_all(i) = NaN;
    sample_count_all(i) = 0;
    sample_count_all_green(i) = 0;
    continue; 
end

    if numel(blue_marker_times) >= 2
        dt_sec = diff(blue_marker_times);
        dt_ms = dt_sec * 1000;

        mean_ms_all(i) = mean(dt_ms);
        std_ms_all(i)  = std(dt_ms);
        median_ms_all(i) = median(dt_ms);
        sample_count_all(i) = numel(dt_ms);
        sample_count_all_green(i) = numel(green_marker_times) - 1;
    else
        mean_ms_all(i) = NaN;
        std_ms_all(i)  = NaN;
        median_ms_all(i) = NaN;
        sample_count_all(i) = 0;
        sample_count_all_green(i) = 0;
        dt_ms = [];
    end



    square_count   = numel(square_marker_times);
    triangle_count = numel(triangle_marker_times);
    classified_total = square_count + triangle_count;

    square_count_all(i)   = square_count;
    triangle_count_all(i) = triangle_count;

    % Sanity check
    if classified_total ~= blue_count
        warning('Current classification mismatch: %d classified vs %d blue markers', ...
            classified_total, blue_count);
    end

    %% ============================================================
    %                CREATE NEW SUMMARY FIGURE (per-file)
    % ============================================================
    figure(fig_summary);
    clf(fig_summary); % clear previous contents each file
    set(fig_summary, 'Position', [200 200 950 700]);

    % 1) TOP PLOT  STACKED BAR CHART FOR MODES A & B
    subplot(3,1,1);
    
    bar_vals = [green_count, blue_count];
    b = bar(bar_vals);
    
    b.FaceColor = 'flat';
    b.CData(1,:) = [0 1 0];   % green
    b.CData(2,:) = [0 0 1];   % blue
    
    set(gca, ...
        'XTick', 1:2, ...
        'XTickLabel', {'Green markers', 'Blue markers'}, ...
        'FontSize', 12);
    
    ylabel('Sample Count');
    title(['Event Sample Sizes  ' fileNames{i}], 'Interpreter','none');
    grid off;

     
    % ============================================================
    % 2) MODE A vs MODE B MARKERS (TRIANGLE vs SQUARE)
    % ============================================================
    subplot(3,1,2);
    
    mode_vals = [triangle_count, square_count];   % Mode A, Mode B
    b2 = bar(mode_vals);
    
    b2.FaceColor = 'flat';
    b2.CData(1,:) = [0 0.8 0.8];   % cyan triangle (Mode A)
    b2.CData(2,:) = [0.6 0.9 0.9]; % lighter cyan square (Mode B)
    
    set(gca, ...
        'XTick', 1:2, ...
        'XTickLabel', {'Mode A (Triangle)', 'Mode B (Square)'}, ...
        'FontSize', 12);
    
    ylabel('Sample Count');
    title('Current-Validated Mode Classification');
    grid off;
    

    % ==================================================
    % 3) BOTTOM PLOT  HISTOGRAM + TEXT STATISTICS (ORIGINAL BLOCK)
    % ==================================================
    subplot(3,1,3);

    % map variables to the original names used in your pasted block
    avg_dt_ms = mean_ms_all(i);
    std_dt_ms = std_ms_all(i);
    median_dt_ms = median_ms_all(i);

    if ~isempty(dt_ms)
        histogram(dt_ms, 20);
        xlabel('Interval Between Blue Markers (ms)');
        ylabel('Count');
        title('Distribution of Mode-A Marker Intervals');
        grid off;

        % ----- Prepare statistics text -----
        N = numel(dt_ms);
        txt = sprintf(['N = %d samples\n' ...
                       'Mean:   %.2f ms\n' ...
                       'Std:    %.2f ms\n' ...
                       'Median: %.2f ms'], ...
                       N, avg_dt_ms, std_dt_ms, median_dt_ms);

        % Position textbox in the TOP RIGHT of histogram
        xlim_vals = xlim;
        ylim_vals = ylim;

        text(xlim_vals(1) + 0.65*(xlim_vals(2)-xlim_vals(1)), ...
             ylim_vals(1) + 0.80*(ylim_vals(2)-ylim_vals(1)), ...
             txt, 'FontSize',12, 'BackgroundColor',[1 1 1 .7], ...
             'EdgeColor',[0 0 0], 'Margin',5);
    else
        text(0.5, 0.5, 'Not enough markers for interval statistics', ...
             'HorizontalAlignment','center','FontSize',12);
        axis off;
    end

    drawnow;

%% === Results Table & Exports ===
% Mode stats Excel (one row per file)
results_table = table( ...
    fileNames(:), ...
    mean_ms_all(:), ...
    std_ms_all(:), ...
    median_ms_all(:), ...
    sample_count_all(:), ...
    sample_count_all_green(:), ...
    square_count_all(:), ...   
    triangle_count_all(:), ... 
    'VariableNames', { ...
        'Filename', ...
        'MeanInterval_ms', ...
        'StdInterval_ms', ...
        'MedianInterval_ms', ...
        'Blue Marker Count', ...
        'Green Marker Count', ...
        'Mode B Count', ...
        'Mode A Count', ...
    } ...
);


output_xlsx = fullfile(pathName, 'mode_statistics.xlsx');
writetable(results_table, output_xlsx);
fprintf('Mode A/B interval statistics saved to: %s\n', output_xlsx);


%% === Peak Voltage Statistics Export (Separate Excel File) ===

peakV_table = table( ...
    fileNames(:), ...
    peakV_mean_all(:), ...
    peakV_std_all(:), ...
    peakV_median_all(:), ...
    avg_maxV_all(:), ...
    std_maxV_all(:), ...
    avg_minV_all(:), ...
    std_minV_all(:), ...
    overall_meanV_all(:), ...
    overall_stdV_all(:), ...
    peakV_count_all(:), ...
    'VariableNames', { ...
        'Filename', ...
        'MeanPeakVoltage_V', ...
        'StdPeakVoltage_V', ...
        'MedianPeakVoltage_V', ...
        'AvgMaxVoltage_V', ...
        'StdMaxVoltage_V', ...
        'AvgMinVoltage_V', ...
        'StdMinVoltage_V', ...
        'OverallMeanVoltage_V', ...
        'OverallStdVoltage_V', ...
        'SampleSize' ...
    } ...
);

output_peakV_xlsx = fullfile(pathName, 'confirmed_peak_voltage_statistics.xlsx');
writetable(peakV_table, output_peakV_xlsx);

fprintf('Confirmed peak voltage statistics saved to: %s\n', output_peakV_xlsx);


%% === Export Raw + STFT Data to Excel (Optional) ===

do_export = strcmpi(strtrim(input('Export signal data to Excel? (true/false): ', 's')), 'true');

if do_export

    % --- User-defined export time range (milliseconds) ---
    t_start_ms = 15;    % <-- set export window start here
    t_end_ms   = 35;   % <-- set export window end here

    % Convert to seconds for masking — no rounding applied
    t_start_s = t_start_ms / 1000;
    t_end_s   = t_end_ms   / 1000;

    % --- Trim raw signals to time window, convert to ms ---
    raw_mask = (T >= t_start_s) & (T <= t_end_s);
    T_export = T(raw_mask) * 1000;   % full double precision, ms
    V_export = V(raw_mask);
    I_export = I(raw_mask);

    % --- Trim STFT to time window, convert to ms ---
    stft_mask            = (t_stft >= t_start_s) & (t_stft <= t_end_s);
    t_stft_export        = t_stft(stft_mask)'    * 1000;   % ms, column vector
    intensity_raw_export = intensity(stft_mask)';           % raw STFT, column vector
    intensity_ma_export  = intensity_ma(stft_mask)';        % moving average, column vector

    % --- Filter blue markers to time window, convert to ms ---
    blue_in_window = blue_marker_times( ...
        blue_marker_times >= t_start_s & blue_marker_times <= t_end_s);
    blue_export = blue_in_window * 1000;   % ms, full double precision

    % --- Filter green markers to time window, convert to ms ---
    green_in_window = green_marker_times( ...
        green_marker_times >= t_start_s & green_marker_times <= t_end_s);
    green_export = green_in_window * 1000;   % ms, full double precision

    % --- NaN-pad all columns to equal length ---
    n_raw   = numel(T_export);
    n_stft  = numel(t_stft_export);
    n_blue  = numel(blue_export);
    n_green = numel(green_export);
    n_rows  = max([n_raw, n_stft, n_blue, n_green]);

    T_export             = [T_export;             nan(n_rows - n_raw,   1)];
    V_export             = [V_export;             nan(n_rows - n_raw,   1)];
    I_export             = [I_export;             nan(n_rows - n_raw,   1)];
    t_stft_export        = [t_stft_export;        nan(n_rows - n_stft,  1)];
    intensity_raw_export = [intensity_raw_export; nan(n_rows - n_stft,  1)];
    intensity_ma_export  = [intensity_ma_export;  nan(n_rows - n_stft,  1)];
    blue_export          = [blue_export;          nan(n_rows - n_blue,  1)];
    green_export         = [green_export;         nan(n_rows - n_green, 1)];

    % --- Build table ---
    export_table = table( ...
        T_export, ...
        V_export, ...
        I_export, ...
        t_stft_export, ...
        intensity_raw_export, ...
        intensity_ma_export, ...
        blue_export, ...
        green_export, ...
        'VariableNames', { ...
            'RawSignal_Time_ms', ...
            'RawSignal_Voltage_V', ...
            'RawSignal_Current_A', ...
            'STFT_Time_ms', ...
            'STFT_Intensity_Raw', ...
            'STFT_Intensity_MovingAvg', ...
            'BlueMarker_Time_ms', ...
            'GreenMarker_Time_ms' ...
        });

    % --- Name file after source CSV + time range ---
    [~, base_name, ~]  = fileparts(fileNames{i});
    range_tag          = sprintf('%gto%gms', t_start_ms, t_end_ms);
    output_signal_xlsx = fullfile(pathName, [base_name '_signal_' range_tag '.xlsx']);

    writetable(export_table, output_signal_xlsx);
    fprintf('Signal export (%g\x2013%g ms) saved to: %s\n', t_start_ms, t_end_ms, output_signal_xlsx);

else
    fprintf('Signal export skipped.\n');
end




%% === NEW FIGURE: STFT Intensity + Current (with green markers) ===
fig_new = 99;   % Choose any unused figure number
figure(fig_new); clf(fig_new);
set(fig_new, 'Position', [120 120 900 500]);


subplot(3,1,1);
plot(T,V,'k','LineWidth', 1.2); hold on;
xlabel('Time [s]');
ylabel('Voltage [V]');
grid off;



% --- Subplot 1: STFT Intensity vs Time with GREEN MARKERS ---
subplot(3,1,2);
subplot(3,1,2);

% Raw STFT intensity (semi-transparent)
plot(t_stft, intensity, 'r', ...
    'LineWidth', 1.2, ...
    'Color', [1 0 0 0.25]); 
hold on;

% Moving-average STFT intensity (used for detection)
plot(t_stft, intensity_ma, 'black', 'LineWidth', 1.8);


% Re-plot GREEN MARKERS at detected run ends
if ~isempty(green_marker_times)
    plot(green_marker_times, ...
         zeros(size(green_marker_times)), ...
         'go', ...
         'MarkerFaceColor','g', ...
         'MarkerSize',6);
end


xlabel('Time [s]');
ylabel('STFT Intensity (Voltage) (a.u.)');
grid off;


% --- Subplot 2: Current vs Time ---
subplot(3,1,3);
plot(T, I, 'blue', 'LineWidth', 1.2); hold on;

% % Plot BLUE markers at y = 0 (time-only event markers)
if ~isempty(blue_marker_times)
    plot(blue_marker_times, ...
         zeros(size(blue_marker_times)), ...
         'bo', ...
         'MarkerFaceColor','cyan', ...
         'MarkerSize',6);
end

% Square markers: 1–2 current peaks
if ~isempty(square_marker_times)
    plot(square_marker_times, ...
         zeros(size(square_marker_times)), ...
         'cs', ...
         'MarkerFaceColor','c', ...
         'MarkerSize',7);
end

% Triangle markers: 3+ current peaks
if ~isempty(triangle_marker_times)
    plot(triangle_marker_times, ...
         zeros(size(triangle_marker_times)), ...
         'c^', ...
         'MarkerFaceColor','c', ...
         'MarkerSize',7);
end

% --- PURPLE DIAMONDS: strong current spikes (true amplitude) ---
% if ~isempty(purple_marker_times)
%     plot(purple_marker_times, ...
%         purple_marker_amps, ...
%          'd', ...
%          'MarkerEdgeColor',[0.45 0 0.75], ...
%          'MarkerFaceColor',[0.6 0 0.9], ...
%          'MarkerSize',7);
% end



xlabel('Time [s]');
ylabel('Current [A]');
title('Current vs Time');
grid off;

linkaxes(findall(gcf,'Type','axes'),'x');

end


%% ============================================================
% FFT Amplitude vs Frequency
% ============================================================

fig_fft = 201;
figure(fig_fft); clf(fig_fft);
set(fig_fft,'Position',[350 250 800 450]);

% --- FFT computation ---
N = length(V);
Y = fft(V);

P2 = abs(Y/N);                 % two-sided spectrum
P1 = P2(1:floor(N/2)+1);       % single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);   % energy correction

% Frequency axis
f = fs*(0:floor(N/2))/N;

% Plot
plot(f, P1, 'k', 'LineWidth', 1.2); hold on
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('FFT Amplitude Spectrum');
grid off

% --- Find maximum amplitude (ignore DC component) ---
[fft_max, idx] = max(P1(2:end));
idx = idx + 1;

f_max = f(idx);

% Marker
plot(f_max, fft_max, 'ro', ...
    'MarkerFaceColor','r', ...
    'MarkerSize',8);

% Floating textbox
txt = sprintf('Max FFT Amplitude: %.4f\nFrequency: %.2f Hz', fft_max, f_max);

xlim_vals = xlim;
ylim_vals = ylim;

text(xlim_vals(1) + 0.65*(xlim_vals(2)-xlim_vals(1)), ...
     ylim_vals(1) + 0.85*(ylim_vals(2)-ylim_vals(1)), ...
     txt, ...
     'FontSize',12, ...
     'BackgroundColor',[1 1 1 0.8], ...
     'EdgeColor',[0 0 0], ...
     'Margin',5);









%% === Helper Functions ===
function filled = fill_invalid_with_nearest(vec)
    invalid = isnan(vec) | isinf(vec);
    valid_idx = find(~invalid);
    if isempty(valid_idx)
        error('No valid data');
    end

    filled = vec;
    filled(invalid) = interp1(valid_idx, vec(valid_idx), find(invalid), 'nearest', 'extrap');
end
