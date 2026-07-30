% Usage:
% - Run the script and select one or more CSV files when prompted.
% - Input files must contain time, current (low and high channels), and voltage data with a 14-line header.
%
% Outputs:
% - Plots of voltage, current, charge, and energy vs. time for each file.
% - FFT analysis plot of current frequency content.
% - Command window display of average power for each file and summary statistics.


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

% === Initialize Storage ===
file_list = {};
energy_list_G = [];
power_list_G = [];
energy_list_I = [];
power_list_I = [];
Zeta_diel_list = [];
C_cell_list = [];
power_list = [];


% === Extract all timestamps ===
timestamp_str = cellfun(@(f) regexp(f, '_3_(\d{17})\.3', 'tokens'), fileNames, 'UniformOutput', false);
timestamp_str = cellfun(@(c) c{1}{1}, timestamp_str, 'UniformOutput', false);
t_all = datetime(timestamp_str, 'InputFormat', 'yyyyMMddHHmmssSSS');
t_all = t_all - t_all(1);


% === Setup Plots ===
color_map = lines(100); % Enough distinct colors
color_idx = 1;

% === Loop Over Files ===
for i = 1:length(fileNames)
    fullpath = fullfile(pathName, fileNames{i});
        header_row = 14;  % Assuming 13 header rows to skip
        data = readmatrix(fullpath, 'NumHeaderLines', header_row);
        
        % Now manually define column indices based on your oscilloscope channel mapping:
        T_idx = 1;   % TIME is first column
        LowI_idx = 2;   % CH2
        HighI_idx = 3;  % CH3
        V_idx = 6;      % MATH1
        
        T = data(:, T_idx);
        lowI = data(:, LowI_idx);
        highI = data(:, HighI_idx);
        V = data(:, V_idx);

        T = fill_invalid_with_nearest(T);
        V = fill_invalid_with_nearest(V);
        lowI = fill_invalid_with_nearest(lowI);
        highI = fill_invalid_with_nearest(highI);

        %Current Fixing
        max_signal = max(lowI);   % Maximum current range of lowI signal
        highI(abs(highI) < max_signal) = 0;
        I = highI + lowI;

       
  


 

%% Plotting   
    % Plot entire V(t) with peaks
    figure(3);
    clf
    subplot(2,1,1)    
    plot(T, V, 'k'); hold on;
    title(['Voltage over Time with Maxima - ' fileNames{i}], 'Interpreter', 'none');
    xlabel('Time [s]'); ylabel('Voltage [V]'); grid on;

    % Plot entire I(t) with peaks
    subplot(2,1,2)
    plot(T, I, 'k'); hold on;
    title(['Current over Time with Maxima - ' fileNames{i}], 'Interpreter', 'none');
    xlabel('Time [s]'); ylabel('Current [A]'); grid on; 

%% === Results Calcs ===

% --- Avg Power Calculation
period = T(end)-T(1);
P_inst = V .* I;
E_avg = trapz(T, P_inst);  % Energy in joules
P_avg = E_avg / period;
P(i) = P_avg;
StdDev = std(P);

E_cumu = cumtrapz(T, P_inst);

end

%% === Results Table ===
hold off;

fprintf('\n-------------------------------------------------------------------------------------------------------------------------------------------------');
fprintf('\n=== Energy, Power Results ===\n');
fprintf('%-30s %-50s \n', 'Plasma Power [W]', 'File Name')

for i = 1:length(fileNames)
    fprintf('%3.4f                        ',P(i))
    fprintf('%-50s\n', fileNames{i});
end

fprintf('\n-----------------------------------------------------------------------------------------------\n');
fprintf('%-30s %3.4f  \n','Average Plasma Power [W]',mean(P))
fprintf('%-30s %3.4f  \n','Standard Deviation [W]', StdDev)

% Prepare output data for CSV
output_table = table(fileNames(:), P(:), 'VariableNames', {'Filename', 'AveragePower_W'});

% Specify output CSV filename
output_csv = fullfile(pathName, 'power_results.csv');

% Write table to CSV
writetable(output_table, output_csv);
fprintf('Power results saved to: %s\n', output_csv);


%% Charge Passed Plot

figure(9); clf; grid on;
plot(T,E_cumu)
ylabel('Energy (J)')

%% Power Over Time Plot
figure(7); clf; hold on;
plot(minutes(t_all), P, '.', 'MarkerSize', 20);
ylabel('Average Power [W]','FontSize',12);
xlabel('Time [min]','FontSize',12);
title('Evolution of Plasma Characterics','FontSize',15);
legend('Power','Location','best')
box on
set(gca, 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k');  % thicker black border


%%
%fftanalysis(I,T,1e6, 50e6,true);
%fftanalysis(V,T,1e6, 50e6,true);


%% === Helper Functions ===
function filled = fill_invalid_with_nearest(vec)
    invalid = isnan(vec) | isinf(vec);
    valid_idx = find(~invalid);
    if isempty(valid_idx)
        error('No valid data');
    end
    
    % Fill invalids by nearest valid using linear interpolation + nearest extrapolation
    filled = vec;
    filled(invalid) = interp1(valid_idx, vec(valid_idx), find(invalid), 'nearest', 'extrap');
end

function [discharge_sum, main_freq_kHz] = fftanalysis(I,T,freq_min, freq_max, qplot)
    % FFT Analysis
        Y = fft(I);
        fs = 1 / mean(diff(T));
        f = (0:length(Y)-1) * fs / length(Y);
        magnitude = real(Y).^2 + imag(Y).^2;

    % === Peak Detection ===
        [peak_vals, locs] = findpeaks(magnitude);
        peak_freqs = f(locs);
        [~, idx] = maxk(peak_vals, 4);
        top_peaks = [peak_freqs(idx)', peak_vals(idx)];
        top_peaks_sorted = sortrows(top_peaks, 1);
        ratios = top_peaks_sorted(:,1) / top_peaks_sorted(1,1);
        % === Sum of Magnitudes from 1 MHz to 50 MHz ===
        
        % Select only frequencies within range
        in_range_idx = (f >= freq_min) & (f <= freq_max);
        mag_in_range = magnitude(in_range_idx);
        
        % Sum of raw FFT magnitude (power-like quantity)
        total_magnitude_sum = sum(mag_in_range);
        discharge_sum = total_magnitude_sum;

        main_freq_kHz = top_peaks_sorted(1,1);


        if qplot == true
                    fprintf('Total FFT magnitude from 1 MHz to 50 MHz: %.4e\n', total_magnitude_sum);
                    figure(6)
                    clf
                    %semilogx(f, magnitude);
                    loglog(f, magnitude)
                    hold on;
                    plot(top_peaks_sorted(1:4,1), top_peaks_sorted(1:4,2), '*', 'MarkerSize', 12);
                    xlabel('Frequency (Hz, log scale)');
                    ylabel('Magnitude');
                    title(sprintf('FFT of Current (Main Freq: %.0f kHz)', main_freq_kHz / 1e3));
                    %xline(freq_min); xline(freq_max)
                    legend('FFT', 'Top 4 Peaks','','', 'Location', 'best');
                    xlim([1e3, 1e9]);
        end
end
