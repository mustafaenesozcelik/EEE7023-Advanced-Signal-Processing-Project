clear; clc; close all;

% 1. Simulation Parameters
fs = 10;                % Sampling frequency (10 Hz)
dt = 1/fs;              % Time step (0.1 seconds)
T_total = 50;           % Total simulation time (seconds)
K = T_total * fs;       % Total number of time steps (500 steps)
b_pos = [0; 0];         % Sensor (Beacon) Location [b_x, b_y]

% 2. Process Model (Constant Velocity) Matrices
% State transition matrix (Phi)
Phi = [1 dt 0  0;
       0  1 0  0;
       0  0 1 dt;
       0  0 0  1];

% Process noise covariance matrix (Q)
q_psd = 0.5; % Acceleration noise power spectral density
Q_sub = [dt^3/3, dt^2/2; 
         dt^2/2, dt];
Q = q_psd * blkdiag(Q_sub, Q_sub);

% 3. Ground Truth Generation
x_true = zeros(4, K);
x_true(:,1) = [10; 2; 20; -1]; % Initial state: [p_x; v_x; p_y; v_y]

for k = 2:K
    % Gaussian process noise for actual movement
    w_k = mvnrnd([0 0 0 0], Q)';
    x_true(:,k) = Phi * x_true(:,k-1) + w_k;
end

% 4. Measurement Simulation (Range-Only & Truncated Cauchy Noise)
z_true = zeros(1, K);   
z_meas = zeros(1, K);   

for k = 1:K
    p_x = x_true(1, k);
    p_y = x_true(3, k);
    z_true(k) = sqrt((p_x - b_pos(1))^2 + (p_y - b_pos(2))^2);
end

signal_power = var(z_true);
target_snr_db = 15;
noise_power = signal_power / (10^(target_snr_db / 10));
gamma = sqrt(noise_power) * 0.5; % Cauchy scale parameter

% Generate Cauchy noise
v_k = gamma * tan(pi * (rand(1, K) - 0.5));
max_error = 30; 
v_k(v_k > max_error) = max_error;
v_k(v_k < -max_error) = -max_error;

% Range cannot be a negative distance in the physical world
z_meas = abs(z_true + v_k);

% 5. SIR Particle Filter Initialization
N_particles = 2000;     % Increased particle count for range-only tracking stability
x_est = zeros(4, K);    % State estimates
P_est = zeros(4, 4, K); % Covariance estimates (required for NEES)

P0 = diag([5^2, 1^2, 5^2, 1^2]); % Initial uncertainty
initial_state = x_true(:,1); 
particles = repmat(initial_state, 1, N_particles) + chol(P0)' * randn(4, N_particles);
weights = ones(1, N_particles) / N_particles; 

% 6. SIR Particle Filter Loop
for k = 1:K
    % A. PREDICTION STEP
    if k > 1
        process_noise = chol(Q)' * randn(4, N_particles);
        particles = Phi * particles + process_noise;
    end
    
    % B. UPDATE STEP (WEIGHT CALCULATION)
    p_x_particles = particles(1, :);
    p_y_particles = particles(3, :);
    z_expected = sqrt((p_x_particles - b_pos(1)).^2 + (p_y_particles - b_pos(2)).^2);
    
    innovation = z_meas(k) - z_expected;
    
    % Likelihood using Cauchy PDF
    likelihood = (1 / (pi * gamma)) * (1 ./ (1 + (innovation / gamma).^2));
    
    % Add a tiny value to prevent weights from becoming exactly zero
    likelihood = likelihood + 1e-12; 
    
    weights = weights .* likelihood;
    weights = weights / sum(weights);
    % C. STATE ESTIMATION
    x_est(:, k) = sum(particles .* weights, 2);
    
    % Calculate Covariance for NEES evaluation
    diff = particles - repmat(x_est(:, k), 1, N_particles);
    P_est(:, :, k) = (diff .* repmat(weights, 4, 1)) * diff';
    
    % D. RESAMPLING STEP
    N_eff = 1 / sum(weights.^2);
    if N_eff < N_particles / 2
        cumulative_sum = cumsum(weights);
        step = 1 / N_particles;
        u = (0:step:1-step) + rand() * step;
        
        resampled_particles = zeros(4, N_particles);
        idx = 1;
        for i = 1:N_particles
            while u(i) > cumulative_sum(idx)
                idx = idx + 1;
            end
            resampled_particles(:, i) = particles(:, idx);
        end
        
        particles = resampled_particles;
        weights = ones(1, N_particles) / N_particles;
    end
end

% 7. Performance Metrics (RMSE and NEES)
% Position RMSE Calculation
pos_error_x = x_true(1,:) - x_est(1,:);
pos_error_y = x_true(3,:) - x_est(3,:);
rmse_k = sqrt(pos_error_x.^2 + pos_error_y.^2);
mean_rmse = mean(rmse_k);

% NEES Calculation
nees_k = zeros(1, K);
for k = 1:K
    err = x_true(:,k) - x_est(:,k);
    % Pseudo-inverse is used to prevent singularity issues in covariance
    nees_k(k) = err' * pinv(P_est(:,:,k)) * err;
end
mean_nees = mean(nees_k);

% Output results to Command Window
fprintf('--- Particle Filter Performance ---\n');
fprintf('Mean Position RMSE: %.2f meters\n', mean_rmse);
fprintf('Mean NEES: %.2f\n', mean_nees);

% 8. Comprehensive Visualization
figure('Position', [100, 100, 1200, 800]);

% Trajectory Plot
subplot(2,2,1);
plot(x_true(1,:), x_true(3,:), 'b-', 'LineWidth', 2); hold on;
plot(x_est(1,:), x_est(3,:), 'g--', 'LineWidth', 2);
plot(b_pos(1), b_pos(2), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
title('2D Target Tracking (Range-Only)');
xlabel('X Position (m)'); ylabel('Y Position (m)');
legend('Ground Truth', 'PF Estimate', 'Sensor'); grid on;

% Measurements Plot
subplot(2,2,2);
plot(1:K, z_meas, 'r.', 'MarkerSize', 5); hold on;
plot(1:K, z_true, 'k-', 'LineWidth', 1.5);
title('Truncated Cauchy Measurements vs. True Range');
xlabel('Time Step (k)'); ylabel('Range (m)');
legend('Noisy Measurements', 'True Range'); grid on;

% RMSE Plot
subplot(2,2,3);
plot(1:K, rmse_k, 'm-', 'LineWidth', 1.5);
title(sprintf('Position RMSE (Mean: %.2f m)', mean_rmse));
xlabel('Time Step (k)'); ylabel('RMSE (meters)'); grid on;

% NEES Plot
subplot(2,2,4);
plot(1:K, nees_k, 'c-', 'LineWidth', 1.5); hold on;
yline(4, 'k--', 'Theoretical Mean (4 DoF)', 'LineWidth', 1.5); 
title(sprintf('NEES (Mean: %.2f)', mean_nees));
xlabel('Time Step (k)'); ylabel('NEES'); grid on;