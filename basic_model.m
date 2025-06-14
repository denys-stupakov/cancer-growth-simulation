clc; clear;

syms t;

r = 0.01;
alpha = 0.1;
beta = 10;
mu = 0.8;
T_max = 7;
T_0 = 0.5;
C_0 = 8;

T_before (t) = T_0 * exp(r*t);

t_max_eq = T_before (t) == T_max;
t_min = 0;
t_max = double(solve(t_max_eq, t));
t_max_i = t_max;

t_treat = (t_min + t_max) / 2;
C (t) = mu/beta + (C_0 - mu/beta)*exp(-beta * t);
T_after (t) = exp((r-alpha)/beta*(t-t_treat) + ... 
                  (alpha*beta*C_0 - alpha*mu)/beta^2 * exp(-beta*(t-t_treat)) + ...
                  log(T_before(t_treat)) + ...
                  (alpha*mu-alpha*beta*C_0)/beta^2);

T_before_num = matlabFunction(T_before);
T_after_num = matlabFunction(T_after);


T(t) = piecewise(t < t_treat, T_before(t), t >= t_treat, T_after(t));
eps = 1;

last_t_treat = NaN;

while (t_max - t_min) > eps
    t_treat = (t_max + t_min) / 2;
    T_after(t) = exp((r-alpha)/beta*(t-t_treat) + ... 
                     (alpha*beta*C_0 - alpha*mu)/beta^2 * exp(-beta*(t-t_treat)) + ...
                     log(T_before(t_treat)) + ...
                     (alpha*mu-alpha*beta*C_0)/beta^2);
    T(t) = piecewise(t < t_treat, T_before(t), t >= t_treat, T_after(t));
    S = solve(T(t) == T_max, t, 'Real', true);
    
    if isempty(S)
        last_t_treat = t_treat;
        t_min = t_treat;
    else
        break;
    end
end

t_treat = last_t_treat;
T_after (t) = exp((r-alpha)/beta*(t-t_treat) + ... 
                  (alpha*beta*C_0 - alpha*mu)/beta^2 * exp(-beta*(t-t_treat)) + ...
                  log(T_before(t_treat)) + ...
                  (alpha*mu-alpha*beta*C_0)/beta^2);
T_after_num = matlabFunction(T_after);
T(t) = piecewise(t < t_treat, T_before(t), t >= t_treat, T_after(t));
hold on;
fplot(T_before_num, [0, last_t_treat], 'LineWidth', 2);
fplot(T_after_num, [last_t_treat, 2*t_max], 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Tumor Volume (cm^3)');
xline(t_treat, ':b', 'Treatment date', 'LineWidth', 3);
legend('T(t) before treatment', 'T(t) after treatment', 'Treatment date');
hold off;
