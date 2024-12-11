clear all

% Load the earthquake catalog
filename = 'usgs1210.csv';
data = readtable(filename, 'VariableNamingRule', 'preserve');

% Inspect the column names
disp(data.Properties.VariableNames);

% Extract latitude, longitude, magnitude, and depth
latitude = data.("latitude");
longitude = data.("longitude");
magnitude = data.("mag");
depth = data.("depth");

% Scale the circle size based on magnitude (optional: multiply for better visualization)
circleSizes = magnitude .^ 2 * 10; % Adjust the multiplier as needed for clarity

% Plot the earthquake events on a map
figure;

scatterHandle = geoscatter(latitude, longitude, circleSizes, depth, 'filled', ...
    'MarkerEdgeColor', 'k'); % Add black edge to the circles

h = colorbar; 
colormap('summer');
h.Label.String = 'Depth (km)';
title('Earthquake Events (Color: Depth, Size: Magnitude)');
set(gca, 'FontSize', 12);
saveas(gcf, 'Earthquake_Events_Map.png');



% set correction to maximum curvature method !!! 
correction = 0.2;

% set the upper bound of the fitted magnitudes!!! (Mc + ??)
upper = 1.5;

% Calculate the b-value using the Gutenberg-Richter relationship
% Define magnitude bins
m_min = floor(min(magnitude) * 10) / 10; % Round down to the nearest 0.1
m_max = ceil(max(magnitude) * 10) / 10;  % Round up to the nearest 0.1
mag_bins = m_min:0.1:m_max;

% Compute the number of events for each bin
num_events = histcounts(magnitude, [mag_bins inf]);

% Maximum Curvature Method: Find the magnitude of completeness (M_c)
[~, max_idx] = max(num_events); % Find the index of the maximum number of events
magnitude_of_completeness = mag_bins(max_idx) + correction; % Magnitude of completeness, usually has a +0.2 to +0.4 correction

% Cumulative number of events
cum_events = cumsum(num_events, 'reverse');

% Perform a linear regression to estimate the b-value (only for magnitudes >= M_c and < M_c+??)
% 1.5 is a experimentally determined value according to your data
selected_bins = (mag_bins >= magnitude_of_completeness) & (mag_bins < magnitude_of_completeness + upper);
mag_bins_fit = mag_bins(selected_bins);
cum_events_fit = cum_events(selected_bins); % Adjust for bin size

log_cum_events_fit = log10(cum_events_fit);
[coefficients, ~] = polyfit(mag_bins_fit, log_cum_events_fit, 1);
b_value = -coefficients(1);

% Plot the Gutenberg-Richter relationship
figure;

% Plot all data points
plot(mag_bins(1:end), log10(cum_events), 'o', 'MarkerFaceColor', 'b');
hold on;

% Highlight the magnitude of completeness (M_c)
xline(magnitude_of_completeness, '--k', ['M_c = ' num2str(magnitude_of_completeness, '%.2f')]);

% Plot the linear fit for M >= M_c
plot(mag_bins_fit, polyval(coefficients, mag_bins_fit), '-r', 'LineWidth', 2);
hold off;

title('Gutenberg-Richter Relationship');
xlabel('Magnitude');
ylabel('Log_{10}(Cumulative Events)');
set(gca, 'FontSize', 12);

equation_str = ['N = ' num2str(coefficients(2), '%.2f') ' - ' num2str(b_value, '%.2f') ' \cdot M'];
x_position = mean(mag_bins); % Adjust X position as needed
y_position = max(log10(cum_events)) - 0.5; % Adjust Y position as needed
text(x_position, y_position, equation_str, 'FontSize', 12, 'Interpreter', 'tex', 'Color', 'red');

% Save the plot
saveas(gcf, 'Gutenberg_Richter_Relationship_with_Mc.png');

% Display results
disp(['Magnitude of Completeness (M_c): ' num2str(magnitude_of_completeness, '%.2f')]);
disp(['b-value: ' num2str(b_value, '%.2f')]);

close all