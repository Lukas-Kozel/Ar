% Define the transfer function
s = tf('s');
P0_nominal = 20 / ((0.2 * s + 1) * (0.4 * s + 1));

% Generate samples for uncertain system
num_samples = 6; % Number of samples
P0_samples = usample(P0_uncertain, num_samples); % Sampled systems

% Frequency range for the plot
f = logspace(-2, 2, 1000); % Frequency from 10^-2 to 10^2 rad/s

% Plot LAFCH (magnitude only) for P0_nominal and samples
figure;
hold on;

% Plot P0_nominal
bodemag(P0_nominal, f, 'r--'); % Nominal system in blue dashed line

% Plot sampled uncertain systems
bodemag(P0_samples(:,:,1), f, 'b-'); % Sample 1 (Red solid line)
bodemag(P0_samples(:,:,2), f, 'g-'); % Sample 2 (Green solid line)
bodemag(P0_samples(:,:,3), f, 'm-'); % Sample 3 (Magenta solid line)
bodemag(P0_samples(:,:,4), f, 'c-'); % Sample 4 (Cyan solid line)
bodemag(P0_samples(:,:,5), f, 'k-'); % Sample 5 (Black solid line)
bodemag(P0_samples(:,:,6), f, 'y-');
% Add grid and title
grid on;

% Add legend for nominal and sampled systems
legend('P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'FontSize', 10, 'Location', 'Best');

hold off;

% Plot LAFCH (magnitude only) for P0_nominal and samples
figure;
hold on;

% Plot P0_nominal
nyquist(P0_nominal, f, 'r--'); % Nominal system in blue dashed line

% Plot sampled uncertain systems
nyquist(P0_samples(:,:,1), f, 'b-'); % Sample 1 (Red solid line)
nyquist(P0_samples(:,:,2), f, 'g-'); % Sample 2 (Green solid line)
nyquist(P0_samples(:,:,3), f, 'm-'); % Sample 3 (Magenta solid line)
nyquist(P0_samples(:,:,4), f, 'c-'); % Sample 4 (Cyan solid line)
nyquist(P0_samples(:,:,5), f, 'k-'); % Sample 5 (Black solid line)
nyquist(P0_samples(:,:,6), f, 'y-'); % Sample 5 (Black solid line)

% Add grid and title
title('Nyquist', 'FontSize', 14);
xlabel('Re', 'FontSize', 12);
ylabel('Im', 'FontSize', 12);

% Add legend for nominal and sampled systems
legend('P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'FontSize', 10, 'Location', 'Best');

hold off;
