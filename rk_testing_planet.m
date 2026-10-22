% Explicit Midpoint
emp_struct.A = [0, 0; .5, 0];
emp_struct.B = [0, 1];
emp_struct.C = [0; 1/2];

% Heun's 3rd-order
h_struct.A = [0, 0, 0; 1/3, 0, 0; 0, 2/3, 0];
h_struct.B = [1/4, 0, 3/4];
h_struct.C = [0; 1/3; 2/3];

% Classic 4th-order
rk4_struct.A = [0, 0, 0, 0; 1/2, 0, 0, 0; 0, 1/2, 0, 0; 0, 0, 1, 0];
rk4_struct.B = [1/6, 1/3, 1/3, 1/6];
rk4_struct.C = [0; 1/2; 1/2; 1];

% Constants
tspan = [0 30];
h_ref = .01;
X0 = 1;

orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;
x0 = 8;
y0 = 0;
dxdt0 = 0;
dydt0 = 1.5;

V0 = [x0;y0;dxdt0;dydt0];
t_range = linspace(0,30,100);
V_list = compute_planetary_motion(t_range,V0,orbit_params);
    
    

gravity_wrapper = @(t, V) gravity_rate_func(t, V, orbit_params);


% [emp_t_list,emp_X_list,emp_h_avg, emp_num_evals] = explicit_RK_fixed_step_integration(gravity_wrapper,tspan,V0,h_ref,emp_struct);
% [h_t_list,h_X_list,h_h_avg, h_num_evals] = explicit_RK_fixed_step_integration(gravity_wrapper,tspan,V0,h_ref,h_struct);
% [rk4_t_list,rk4_X_list,rk4_h_avg, rk4_num_evals] = explicit_RK_fixed_step_integration(gravity_wrapper,tspan,V0,h_ref,rk4_struct);
% 
% figure();
% axis equal; axis square;
% axis([-20,20,-20,20])
% hold on;
% plot(0,0,'ro','markerfacecolor','r','markersize',5);
% plot(emp_X_list(1, :),emp_X_list(2, :),'b');
% plot(h_X_list(1, :),h_X_list(2, :),'g');
% plot(rk4_X_list(1, :),rk4_X_list(2, :),'r', LineWidth=3);
% plot(V_list(:,1),V_list(:,2),'k');
% hold off;
% 
% f=figure();
% tcl = tiledlayout(2, 1);
% 
% nexttile(tcl);
% hold on
% 
% plot(t_range, V_list(:, 1), 'k', LineWidth=5);
% plot(emp_t_list, emp_X_list(1, :), 'o',markerfaceColor='r',MarkerEdgeColor='r', LineWidth=2, MarkerIndices=1:200:length(emp_X_list(1, :)));
% plot(h_t_list, h_X_list(1, :),'o', markerfaceColor='b',markeredgeColor='b', LineWidth=2, MarkerIndices=75:200:length(h_X_list(1, :)));
% plot(rk4_t_list, rk4_X_list(1, :), 'o', markerfaceColor='g',markeredgeColor='g', LineWidth=2, MarkerIndices=150:200:length(rk4_X_list(1, :)));
% xlabel("Time (-)", Interpreter="latex", FontSize=14)
% ylabel("x position (-)", Interpreter="latex", FontSize=14)
% hold off;
% 
% nexttile(tcl);
% hold on
% plot(t_range, V_list(:, 2), 'k', LineWidth=5);
% plot(emp_t_list, emp_X_list(2, :), 'o',markerfaceColor='r',MarkerEdgeColor='r', LineWidth=2, MarkerIndices=1:200:length(emp_X_list(2, :)));
% plot(h_t_list, h_X_list(2, :), 'o', markerfaceColor='b',markeredgeColor='b', LineWidth=2, MarkerIndices=75:200:length(h_X_list(2, :)));
% plot(rk4_t_list, rk4_X_list(2, :), 'o', markerfaceColor='g',markeredgeColor='g', LineWidth=2, MarkerIndices=150:200:length(rk4_X_list(2, :)));
% xlabel("Time (sec)", Interpreter="latex", FontSize=14)
% ylabel("y position (-)", Interpreter="latex", FontSize=14)
% hold off;
% 
% sgtitle("X and Y Positions Over Time", Interpreter='latex', FontSize=15)
% hl = legend(["Ground Truth" "Explicit Midpoint" "Heun's 3rd Order" 'Runge-Kutta 4th Order'],Interpreter='latex', FontSize=15);
% hl.Layout.Tile = 'South';
% 
% exportgraphics(f, 'plots/x_y_pos.png', 'Resolution',300)


%% 



t_ref = [0, 3];
h_step = logspace(-2, -1, 100);
t_range = linspace(0,3,100);

V_list = compute_planetary_motion(t_range,V0,orbit_params);


XA = V0;
gt = V_list(end,:);

emp_err_list = [];
h_err_list = [];
rk4_err_list = [];


for i = 1:length(h_step)
    [emp_t_list,emp_X_list,emp_h_avg, emp_num_evals] = explicit_RK_fixed_step_integration(gravity_wrapper,t_ref,V0,h_step(i),emp_struct);
    [h_t_list,h_X_list,h_h_avg, h_num_evals] = explicit_RK_fixed_step_integration(gravity_wrapper,t_ref,V0,h_step(i),h_struct);
    [rk4_t_list,rk4_X_list,rk4_h_avg, rk4_num_evals] = explicit_RK_fixed_step_integration(gravity_wrapper,t_ref,V0,h_step(i),rk4_struct);
    
    XB_emp = emp_X_list(:,end);
    XB_h = h_X_list(:,end);
    XB_rk4 = rk4_X_list(:,end);

    emp_err_list = [emp_err_list, norm(XB_emp - gt')];
    h_err_list = [h_err_list, norm(XB_h - gt')];
    rk4_err_list = [rk4_err_list, norm(XB_rk4 - gt')];

end
% filtered_h_step = h_step(h_step <= 1e-2);
% filt_length = length(filtered_h_step);
% 
% filtered_fe_err = emp_err_list(1:filt_length);
% filtered_mp_err = h_err_list(1:filt_length);
% 
% 
% 
% filtered_fe_nums = fe_evals_list(fe_evals_list >= 1e1);
% filt_fe_num_len = length(filtered_fe_nums);
% filtered_mp_nums = mp_evals_list(mp_evals_list >= 1e1);
% filt_mp_num_len = length(filtered_mp_nums);
% 
% filt_mp_error_2 = h_err_list(1:filt_mp_num_len);
% filt_fe_error_2 = emp_err_list(1:filt_fe_num_len);
% 
% mp_p = polyfit(log(filtered_h_step), log(filtered_mp_err), 1)
% fe_p = polyfit(log(filtered_h_step), log(filtered_fe_err), 1)
% 
% 
% f_func1_mp = exp(polyval(mp_p, log(filtered_h_step)));
% f_func1_fe = exp(polyval(fe_p, log(filtered_h_step)));



emp_p = polyfit(log(h_step), log(emp_err_list), 1)
f_emp = exp(polyval(emp_p, log(h_step)));

h_p = polyfit(log(h_step), log(h_err_list), 1)
f_h = exp(polyval(h_p, log(h_step)));

rk4_p = polyfit(log(h_step), log(rk4_err_list), 1)
f_rk4 = exp(polyval(rk4_p, log(h_step)));


figure();
scatter(h_step, emp_err_list, 60,'filled',"r")
hold on
scatter(h_step, h_err_list,60,'filled',"b")
scatter(h_step, rk4_err_list,60,'filled',"g")

plot(h_step, f_emp, '-', LineWidth=3, Color="k")
plot(h_step, f_h, '-', LineWidth=3, Color="k")
plot(h_step, f_rk4, '-', LineWidth=3, Color="k")
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
exportgraphics(gca, 'plots/global_error.png', 'Resolution',300)