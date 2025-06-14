clc; clear;

%% ---------------- Parameters ---------------- %%
% Tumor growth before treatment
r = 0.1;             % tumor growth rate
T_0 = 0.1;           % initial tumor size
T_max = 5;           % critical threshold for tumor size

% Treatment exposure and model constants
Kg = 0.1;            % growth rate of surviving tumor cells
Kd = 0.4;            % rate of killing surviving tumor cells
d = 0.05;            % clearance rate of dying tumor cells
Exposure = 0.3;      % constant drug exposure level

% Derived parameter for simplified dynamics
gamma = Kg - Kd * Exposure;

%% ---------------- Tumor growth before treatment ---------------- %%
T_before = @(t) T_0 * exp(r * t);    % exponential growth model before treatment

% Find time when tumor hits critical size
t_max = fzero(@(t) T_before(t) - T_max, 10);  

%% ---------------- Binary search for optimal treatment time ---------------- %%
t_min = 0;
eps_val = 0.01;
last_t_treat = NaN;

while (t_max - t_min) > eps_val
    t_treat = (t_min + t_max) / 2;
    S0 = T_before(t_treat);   % viable tumor size at treatment
    D0 = 0;                   % initially no dying cells

    % Simulate analytically after treatment
    t_check = linspace(0, 100, 1000);
    S_check = S0 * exp(gamma * t_check);
    D_check = (Kd * Exposure * S0 / (gamma + d)) * (exp(gamma * t_check) - exp(-d * t_check));
    T_check = S_check + D_check;

    if all(T_check < T_max)
        last_t_treat = t_treat;
        t_min = t_treat;
    else
        t_max = t_treat;
    end
end

%% ---------------- Final Simulation for Best Case ---------------- %%
% Before treatment
t_full_before = linspace(0, last_t_treat, 200);
T_before_vals = T_before(t_full_before);

% After treatment (analytical)
t_post = linspace(0, 100, 1000);
S_post = T_before(last_t_treat) * exp(gamma * t_post);
D_post = (Kd * Exposure * T_before(last_t_treat) / (gamma + d)) * (exp(gamma * t_post) - exp(-d * t_post));
T_post_vals = S_post + D_post;
t_post_shifted = t_post + last_t_treat;

%% ---------------- Early treatment simulation ---------------- %%
t_early = last_t_treat - 2;
if t_early < 0
    warning('t_early < 0, set to 0');
    t_early = 0;
end

S0_early = T_before(t_early);
S_early = S0_early * exp(gamma * t_post);
D_early = (Kd * Exposure * S0_early / (gamma + d)) * (exp(gamma * t_post) - exp(-d * t_post));
T_early_vals = S_early + D_early;
t_early_shifted = t_post + t_early;

%% ---------------- Plot everything ---------------- %%
figure('Position', [400, 100, 1200, 800]);

plot(t_full_before, T_before_vals, 'b', 'LineWidth', 2); hold on;
plot(t_post_shifted, T_post_vals, 'r', 'LineWidth', 2);
plot(t_post_shifted, S_post, 'm--', 'LineWidth', 1.5);  % Plot S(t) for clarity
plot(t_early_shifted, T_early_vals, 'g--', 'LineWidth', 2);

yline(T_max, 'k--', 'Critical size');
yline(T_0, 'k--', 'Initial size');
yline(T_before(last_t_treat), 'k--', ['Optimal treatment start: t = ' num2str(last_t_treat, '%.2f')], 'LineWidth', 2);

xlabel('Time (t)');
ylabel('Tumor size (T)');
title(['S+D model (Analytical): Optimal treatment start — t = ' num2str(last_t_treat, '%.2f')]);
legend('Before treatment', ...
       'After optimal treatment', ...
       'Surviving cells S(t)', ...
       'Early treatment (2 units sooner)', ...
       'Critical size', ...
       'Initial size', ...
       'Best moment');

grid on;
axis([-10 140 -1 6]);
