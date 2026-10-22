BT_struct = butcher_tableau("HeunEuler");



% Constants
tref = 10;
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
% t_range = linspace(0,30,100);
% V_list = compute_planetary_motion(t_range,V0,orbit_params);
% 
    

gravity_wrapper = @(t, V) gravity_rate_func(t, V, orbit_params);

XB1_err_list = [];
XB2_err_list = [];
diff_err_list = [];
gt_err_list = [];


for i = 1:length(h_step)
    gt = gravity_wrapper(t_ref + h_step(i), V0);

   [XB1, XB2, num_evals] = RK_step_embedded(gravity_wrapper,t_ref,V0,h_step(i),BT_struct);

    XB1_err_list = [XB1_err_list, norm(XB1 - gt')];
    XB2_err_list = [XB2_err_list, norm(XB2 - gt')];
    diff_err_list = [diff_err_list, norm(XB1 - XB2)];

    gt_err_list = [gt_err_list, norm(gt' - V0)];


end

% filtered_h_step = h_step(h_step <= 1e0);
% filt_length = length(filtered_h_step);
% 
% filtered_fe_err = XB1_err_list(1:filt_length);
% filtered_mp_err = XB2_err_list(1:filt_length);
% 
% filtered_gt_err = func1_gt_err_list(1:filt_length);
% 
% mp_p = polyfit(log(filtered_h_step), log(filtered_mp_err), 1)
% fe_p = polyfit(log(filtered_h_step), log(filtered_fe_err), 1)
% 
% gt_p = polyfit(log(filtered_h_step), log(filtered_gt_err), 1)
% 
% f_func1_mp = exp(polyval(mp_p, log(filtered_h_step)));
% f_func1_fe = exp(polyval(fe_p, log(filtered_h_step)));
% f_func1_gt = exp(polyval(gt_p, log(filtered_h_step)));


figure();
scatter(h_step, XB1_err_list, 150, "filled", "r")
hold on
scatter(h_step, XB2_err_list, 100, "filled", "g")
scatter(h_step, diff_err_list, 50, "filled", "b");
scatter(h_step, gt_err_list, 50, "filled", "k");

% plot(filtered_h_step, f_func1_mp, LineWidth=2, Color="k")
% plot(filtered_h_step, f_func1_fe, LineWidth=2, Color="k")
% plot(filtered_h_step, f_func1_gt, LineWidth=2, Color="k")
hold off
set(gca,'xscale','log')
set(gca,'yscale','log')
% legend(["Explicit Midpoint Measured","Forward Euler Measured", "Ground Truth Measured"], Location = "southeast", Interpreter="latex", FontSize=16)
xlabel("Time Step $h$ (-)", Interpreter="latex", FontSize=16)
ylabel("Local Truncation Error (-)", Interpreter="latex", FontSize=16)
% axis([0 4*pi -1.5 1.5])
title("Local Truncation Error vs. Step Size: Explicit Methods", Interpreter='latex', FontSize=18)
ax = gca;
ax.FontSize = 16;
ax.XAxis.TickLabelInterpreter = 'Latex';
ax.YAxis.TickLabelInterpreter = 'Latex';
% exportgraphics(gca, 'plots/explicit_methods_local_error_testfunc1.png', 'Resolution',300)