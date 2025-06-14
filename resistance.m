r = 0.075;
T0 = 2.5;
c0 = 0.1;
c_min = 0.06;
beta = 0.003;
Tcrit = 0.1;


t = linspace(0, 700, 6000);
T = T0 .* exp((r - c_min) .* t + ((c0 - c_min) ./ beta) .* (exp(-beta .* t) - 1));

idx = find(T <= Tcrit, 1);

if ~isempty(idx)
    fprintf('Tumor was killed on %.2f\n', t(idx));
    treatment_success = true;
else
    fprintf('Tumor was not killed\n');
    treatment_success = false;
end

figure;
plot(t, T, 'b', 'LineWidth', 2); hold on;
yline(Tcrit, '--k', 'Critical volume');
xlabel('Time (days)');
ylabel('Tumor Volume (cm^3)');
grid on;

if treatment_success
    xline(t(idx), ':r', 'Killing date', 'LineWidth', 3);
    legend('Tumor volume', 'Critical volume', 'Killing date', 'Location', 'NorthEast');
else
    legend('Tumor volume', 'Critical volume', 'Location', 'NorthEast');
end
