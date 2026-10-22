%This function computes the value of X at the next time step
%for any arbitrary embedded RK method
%also computes the next step size to use, and whether or not to
%accept/reject the projected value of X(t+h)
%(for variable time step methods)
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%p: how error scales with step size (error = k*hˆp)
%error_desired: the desired local truncation error at each step
%OUTPUTS:
%XB: the approximate value for X(t+h)
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
%h_next: the time-step size at the next iteration
%redo: False if the estimated error was less than error_desired
% True if the estimated error was larger than error_desired
function [X1_list, X2_list, num_evals] = explicit_RK_variable_step_test...
(rate_func_in,tspan,X0,h_ref,BT_struct)
    n_steps = ceil((tspan(2)-tspan(1))/h_ref);

    h_avg = (tspan(2) - tspan(1))/n_steps;
    num_evals = 0;
    XA = X0;
    t_list = linspace(tspan(1), tspan(2), n_steps+1);
    X1_list = [XA];
    X2_list = [XA];
    
    for i = 1:n_steps
        [XB1, XB2, evals] = RK_step_embedded(rate_func_in, t_list(i), XA, h_avg, BT_struct);
        num_evals = num_evals + evals;
        X1_list = [X1_list, XB1];
        X2_list = [X2_list, XB2];
        
        XA = XB2;
    end

end