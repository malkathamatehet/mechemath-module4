BT_struct = butcher_tableau("DormandPrince");



% Constants
t_ref = 0;
h_step = logspace(-3, 1, 100);
% X0 = 1;

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

X0 = compute_planetary_motion(t_ref, V0, orbit_params);

for i = 1:length(h_step)
    gt = compute_planetary_motion(t_ref + h_step(i), X0, orbit_params);

   [XB1, XB2, num_evals] = RK_step_embedded(gravity_wrapper,t_ref,X0,h_step(i),BT_struct);

    XB1_err_list = [XB1_err_list, norm(XB1 - gt')];
    XB2_err_list = [XB2_err_list, norm(XB2 - gt')];
    diff_err_list = [diff_err_list, norm(XB1 - XB2)];

    gt_err_list = [gt_err_list, norm(gt' - X0)];


end

% filtered_h_step = h_step(h_step <= 1e0);
% filt_length = length(filtered_h_step);
% 
% filtered_fe_err = XB1_err_list(1:filt_length);
% filtered_mp_err = XB2_err_list(1:filt_length);
% 
% filtered_gt_err = func1_gt_err_list(1:filt_length);
% 
XB1_p = polyfit(log(h_step), log(XB1_err_list), 1)
XB2_p = polyfit(log(h_step), log(XB2_err_list), 1)
diff_p = polyfit(log(h_step), log(diff_err_list), 1)


gt_p = polyfit(log(h_step), log(gt_err_list), 1)
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
legend(["XB1 Measured","XB2 Measured", "XB1-XB2 Measured", "Ground Truth Measured"], Location = "northwest", Interpreter="latex", FontSize=16)
xlabel("Time Step $h$ (-)", Interpreter="latex", FontSize=16)
ylabel("Local Truncation Error (-)", Interpreter="latex", FontSize=16)
% axis([0 4*pi -1.5 1.5])
title("Local Truncation Error vs. Step Size: Bogaki", Interpreter='latex', FontSize=18)
ax = gca;
ax.FontSize = 16;
ax.XAxis.TickLabelInterpreter = 'Latex';
ax.YAxis.TickLabelInterpreter = 'Latex';
exportgraphics(gca, 'plots/embedded_local_trunc_error.png', 'Resolution',300)