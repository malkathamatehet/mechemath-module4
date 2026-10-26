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
h_step = logspace(-3, 0, 100);
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



% emp_p = polyfit(log(h_step), log(emp_err_list), 1)
% f_emp = exp(polyval(emp_p, log(h_step)));
% 
% h_p = polyfit(log(h_step), log(h_err_list), 1)
% f_h = exp(polyval(h_p, log(h_step)));
% 
% rk4_p = polyfit(log(h_step), log(rk4_err_list), 1)
% f_rk4 = exp(polyval(rk4_p, log(h_step)));
% 

figure();
scatter(h_step, emp_err_list, 60,'filled',"r")
hold on
scatter(step_size_list, global_error_list ,60,'filled',"b")
hold off
set(gca,'xscale','log')
set(gca,'yscale','log')
% legend(["Explicit Midpoint Measured","Heun's 3rd Order Measured", "Runge-Kutta 4th Order Measured", "Fit Lines"], Location = "southeast", Interpreter="latex", FontSize=14)
xlabel("Time Step - $h$ (-)", Interpreter="latex", FontSize=14)
ylabel("Global Truncation Error (-)", Interpreter="latex", FontSize=14)
%axis([0 4*pi -1.5 1.5])
title("Global Truncation Error vs Step Size", Interpreter='latex', FontSize=18)
ax = gca;
ax.FontSize = 16;
ax.XAxis.TickLabelInterpreter = 'Latex';
ax.YAxis.TickLabelInterpreter = 'Latex';
% exportgraphics(gca, 'plots/global_error.png', 'Resolution',300)


figure();
scatter(emp_eval_list, emp_err_list, 60,'filled',"r")
hold on
scatter(func_eval_list, global_error_list ,60,'filled',"b")
hold off
set(gca,'xscale','log')
set(gca,'yscale','log')
% legend(["Explicit Midpoint Measured","Heun's 3rd Order Measured", "Runge-Kutta 4th Order Measured", "Fit Lines"], Location = "southeast", Interpreter="latex", FontSize=14)
xlabel("Function Evaluations", Interpreter="latex", FontSize=14)
ylabel("Global Truncation Error (-)", Interpreter="latex", FontSize=14)
%axis([0 4*pi -1.5 1.5])
title("Global Truncation Error vs Function Evaluations", Interpreter='latex', FontSize=18)
ax = gca;
ax.FontSize = 16;
ax.XAxis.TickLabelInterpreter = 'Latex';
ax.YAxis.TickLabelInterpreter = 'Latex';
% exportgraphics(gca, 'plots/global_error.png', 'Resolution',300)


figure();
length(step_size_list)
length(failure_list)
scatter(step_size_list, failure_list, 60,'filled',"r")
hold on
set(gca,'xscale','log')
set(gca,'yscale','linear')
% legend(["Explicit Midpoint Measured","Heun's 3rd Order Measured", "Runge-Kutta 4th Order Measured", "Fit Lines"], Location = "southeast", Interpreter="latex", FontSize=14)
xlabel("Time Step - $h$ (-)", Interpreter="latex", FontSize=14)
ylabel("Failure Rate", Interpreter="latex", FontSize=14)
%axis([0 4*pi -1.5 1.5])
title("Global Truncation Error vs Step Size", Interpreter='latex', FontSize=18)
ax = gca;
ax.FontSize = 16;
ax.XAxis.TickLabelInterpreter = 'Latex';
ax.YAxis.TickLabelInterpreter = 'Latex';
% exportgraphics(gca, 'plots/global_error.png', 'Resolution',300)

%% 

% Constants
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
plot(v_t_list, v_X_list(1, :), 'o',markerfaceColor='r',MarkerEdgeColor='r', LineWidth=2);
xlabel("Time (-)", Interpreter="latex", FontSize=14)
ylabel("x position (-)", Interpreter="latex", FontSize=14)
hold off;

nexttile(tcl);
hold on
plot(t_range, V_list(:, 2), 'k', LineWidth=5);
plot(v_t_list, v_X_list(2, :), 'o',markerfaceColor='r',MarkerEdgeColor='r', LineWidth=2);
xlabel("Time (sec)", Interpreter="latex", FontSize=14)
ylabel("y position (-)", Interpreter="latex", FontSize=14)
hold off;

sgtitle("X and Y Positions Over Time", Interpreter='latex', FontSize=15)
hl = legend(["Ground Truth" "Variable Timestep"],Interpreter='latex', FontSize=15);
hl.Layout.Tile = 'South';
% 
% exportgraphics(f, 'plots/x_y_pos.png', 'Resolution',300)

f = figure();
scatter(step_single_error, r_list, 60,'filled',"r");
hold on
set(gca,'xscale','log')
set(gca,'yscale','log')
% legend(["Explicit Midpoint Measured","Heun's 3rd Order Measured", "Runge-Kutta 4th Order Measured", "Fit Lines"], Location = "southeast", Interpreter="latex", FontSize=14)
xlabel("Time Step - $h$ (-)", Interpreter="latex", FontSize=14)
ylabel("Distance", Interpreter="latex", FontSize=14)
%axis([0 4*pi -1.5 1.5])
title("-", Interpreter='latex', FontSize=18)
ax = gca;
ax.FontSize = 16;
ax.XAxis.TickLabelInterpreter = 'Latex';
ax.YAxis.TickLabelInterpreter = 'Latex';