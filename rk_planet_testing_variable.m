BT_struct = butcher_tableau("Bogacki");

tmax = 30;
% Constants
tspan = [0 tmax];
h_ref = .01;
X0 = 1;

orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;
x0 = 8;
y0 = 0;
dxdt0 = 0;
dydt0 = 1.05;

V0 = [x0;y0;dxdt0;dydt0];
t_range = linspace(0,tmax,100);
V_list = compute_planetary_motion(t_range,V0,orbit_params);
gt = V_list(end,:);
    
p = length(BT_struct.C)
error_desired = 1e-3;

gravity_wrapper = @(t, V) gravity_rate_func(t, V, orbit_params);

desired_error_list = logspace(-6, -2, 100);

global_error_list = [];
step_size_list = [];
func_eval_list = [];
failure_list = [];

for i = 1:length(desired_error_list)
    [v_t_list,v_X_list,v_h_avg, v_num_evals, failure_rate] = explicit_RK_variable_step_integration(gravity_wrapper,tspan,V0,h_ref,BT_struct,p,desired_error_list(i));
    step_size_list = [step_size_list, v_h_avg];
    func_eval_list = [func_eval_list, v_num_evals];
    failure_list = [failure_list, failure_rate];

    XB_calc = v_X_list(:,end);
    global_error_list = [global_error_list, norm(XB_calc - gt')];
end


% figure();
% axis equal; axis square;
% axis([-20,20,-20,20])
% hold on;
% plot(0,0,'ro','markerfacecolor','r','markersize',5);
% plot(emp_X_list(1, :),emp_X_list(2, :),'b');
% plot(V_list(:,1),V_list(:,2),'k');
% hold off;


% t_ref = [0, 30];
h_step = logspace(-3, -1.5, 100);
t_range = linspace(tspan(1),tspan(2),100);
t_ref = [tspan(1), tspan(2)];


XA = V0;
% gt = V_list(end,:);

emp_err_list = [];
emp_eval_list = [];

% % Heun's 3rd-order
% h_struct.A = [0, 0, 0; 1/3, 0, 0; 0, 2/3, 0];
% h_struct.B = [1/4, 0, 3/4];
% h_struct.C = [0; 1/3; 2/3];


for i = 1:length(h_step)
    [emp_t_list,emp_X_list,emp_h_avg, emp_num_evals] = explicit_RK_fixed_step_integration(gravity_wrapper,t_ref,V0,h_step(i),BT_struct);
    
    XB_emp = emp_X_list(:,end);
    emp_err_list = [emp_err_list, norm(XB_emp - gt')];

    emp_eval_list = [emp_eval_list, emp_num_evals];
end


filtered_h_step = step_size_list(step_size_list <= 1e-1);
filt_length = length(filtered_h_step);

filtered_global_error = global_error_list(1:filt_length);

emp_p = polyfit(log(h_step), log(emp_err_list), 1)
f_emp = exp(polyval(emp_p, log(h_step)));

v_p = polyfit(log(filtered_h_step), log(filtered_global_error), 1)
f_v = exp(polyval(v_p, log(filtered_h_step)));

figure();
scatter(h_step, emp_err_list, 60,'filled',"r")
hold on
scatter(filtered_h_step, filtered_global_error ,60,'filled',"b")
plot(h_step, f_emp, '-', LineWidth=3, Color="k")
plot(filtered_h_step, f_v, '-', LineWidth=3, Color="k")
hold off
set(gca,'xscale','log')
set(gca,'yscale','log')
legend(["Fixed Timestep","Variable Timestep", "Fit Lines"], Location = "northwest", Interpreter="latex", FontSize=14)
xlabel("Time Step - $h$ (-)", Interpreter="latex", FontSize=14)
ylabel("Global Truncation Error (-)", Interpreter="latex", FontSize=14)
%axis([0 4*pi -1.5 1.5])
title("Global Truncation Error vs Fixed and Variable Step Sizes: Bogacki", Interpreter='latex', FontSize=18)
ax = gca;
ax.FontSize = 16;
ax.XAxis.TickLabelInterpreter = 'Latex';
ax.YAxis.TickLabelInterpreter = 'Latex';
exportgraphics(gca, 'plots/fixed_variable_global_error_timestep_comparison.png', 'Resolution',300)

%% 

filtered_func_eval = func_eval_list(1:filt_length);
v_p = polyfit(log(filtered_func_eval), log(filtered_global_error), 1)
f_v = exp(polyval(v_p, log(filtered_func_eval)));

emp_p = polyfit(log(emp_eval_list), log(emp_err_list), 1)
f_emp = exp(polyval(emp_p, log(emp_eval_list)));

figure();
scatter(emp_eval_list, emp_err_list, 60,'filled',"r")
hold on
scatter(filtered_func_eval, filtered_global_error ,60,'filled',"b")
plot(emp_eval_list, f_emp, '-', LineWidth=3, Color="k")
plot(filtered_func_eval, f_v, '-', LineWidth=3, Color="k")
hold off
set(gca,'xscale','log')
set(gca,'yscale','log')
legend(["Fixed Timestep","Variable Timestep", "Fit Lines"], Location = "northeast", Interpreter="latex", FontSize=14)
xlabel("Function Evaluations (-)", Interpreter="latex", FontSize=14)
ylabel("Global Truncation Error (-)", Interpreter="latex", FontSize=14)
%axis([0 4*pi -1.5 1.5])
title("Global Truncation Error vs Fixed and Variable Function Evaluations: Bogacki", Interpreter='latex', FontSize=18)
ax = gca;
ax.FontSize = 16;
ax.XAxis.TickLabelInterpreter = 'Latex';
ax.YAxis.TickLabelInterpreter = 'Latex';
exportgraphics(gca, 'plots/fixed_variable_global_error_evals_comparison.png', 'Resolution',300)
%% 


figure();
scatter(step_size_list, failure_list*100, 60,'filled',"k")
set(gca,'xscale','log')
set(gca,'yscale','linear')
% legend(["Explicit Midpoint Measured","Heun's 3rd Order Measured", "Runge-Kutta 4th Order Measured", "Fit Lines"], Location = "southeast", Interpreter="latex", FontSize=14)
xlabel("Time Step - $h$ (-)", Interpreter="latex", FontSize=14)
ylabel("Failure Rate ($\%$)", Interpreter="latex", FontSize=14)
%axis([0 4*pi -1.5 1.5])
title("Failure Rate vs Avg. Step Size: Bogacki", Interpreter='latex', FontSize=18)
ax = gca;
ax.FontSize = 16;
ax.XAxis.TickLabelInterpreter = 'Latex';
ax.YAxis.TickLabelInterpreter = 'Latex';
exportgraphics(gca, 'plots/var_st_si_err_ra.png', 'Resolution',300)

%% 

% Constants
tmax = 30;
tspan = [0 tmax];
h_ref = .01;
X0 = 1;

t_range = linspace(tspan(1),tspan(2),100);
V_list = compute_planetary_motion(t_range,V0,orbit_params);
    
    

gravity_wrapper = @(t, V) gravity_rate_func(t, V, orbit_params);
desired_error = 1e-8;

[v_t_list,v_X_list,v_h_avg, v_num_evals, failure_rate] = explicit_RK_variable_step_integration(gravity_wrapper,tspan,V0,h_ref,BT_struct,p,desired_error);
r_list = sqrt(v_X_list(1, :).^2 + v_X_list(2, :).^2);
step_single_error = tspan(1);
step_single_error = [step_single_error, diff(v_t_list)];

% 
figure();
axis equal; axis square;
axis([-20,20,-20,20])
hold on;
plot(0,0,'ro','markerfacecolor','r','markersize',5);
plot(v_X_list(1, :),v_X_list(2, :),'b');
plot(V_list(:,1),V_list(:,2),'k');
hold off;


f=figure();
tcl = tiledlayout(2, 1);

nexttile(tcl);
hold on

plot(t_range, V_list(:, 1), 'k', LineWidth=5);
plot(v_t_list, v_X_list(1, :), 'ro-',markerfaceColor='r',MarkerEdgeColor='r', LineWidth=2, MarkerIndices=1:10:length(v_X_list(1, :)));
xlabel("Time (-)", Interpreter="latex", FontSize=14)
ylabel("x position (-)", Interpreter="latex", FontSize=14)
hold off;

nexttile(tcl);
hold on
plot(t_range, V_list(:, 2), 'k', LineWidth=5);
plot(v_t_list, v_X_list(2, :), 'ro-',markerfaceColor='r',MarkerEdgeColor='r', LineWidth=2, MarkerIndices=1:10:length(v_X_list(2, :)));
xlabel("Time (-)", Interpreter="latex", FontSize=14)
ylabel("y position (-)", Interpreter="latex", FontSize=14)
hold off;

sgtitle("X and Y Positions Over Time", Interpreter='latex', FontSize=15)
hl = legend(["Ground Truth" "Variable Timestep"],Interpreter='latex', FontSize=15);
hl.Layout.Tile = 'South';
% 
% exportgraphics(f, 'plots/variable_step_size_x_y_pos.png', 'Resolution',300)
%% 

f = figure();
plot(r_list, step_single_error,"o", MarkerFaceColor='k', MarkerEdgeColor='k', MarkerIndices=1:100:length(r_list));
hold on
set(gca,'xscale','log')
set(gca,'yscale','log')
ylabel("Time Step - $h$ (-)", Interpreter="latex", FontSize=14)
xlabel("Distance (-)", Interpreter="latex", FontSize=14)
%axis([0 4*pi -1.5 1.5])
title("Avg. Step Size vs Distance from the Sun: Bogacki", Interpreter='latex', FontSize=18)
ax = gca;
ax.FontSize = 16;
ax.XAxis.TickLabelInterpreter = 'Latex';
ax.YAxis.TickLabelInterpreter = 'Latex';
exportgraphics(f, 'plots/dist_from_sun.png', 'Resolution',300)

