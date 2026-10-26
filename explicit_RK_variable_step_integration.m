%Runs numerical integration arbitrary RK method using variable time steps
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%p: how error scales with step size (error = k*hˆp)
%error_desired: the desired local truncation error at each step
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals, failure_rate] = explicit_RK_variable_step_integration(rate_func_in,tspan,X0,h_ref,BT_struct,p,error_desired)

    num_evals = 0;
    XA = X0;
    curr_t = tspan(1);
    t_list = [curr_t];
    X_list = [XA];
    curr_h = h_ref;
    h_list = [curr_h];

    t_next = curr_t + curr_h;

    fail = 0;
    
    if t_next > tspan(2)
        t_next = tspan(2);
        curr_h = tspan(2)-curr_t;
    end

    % while curr_t <= tspan(2)
    while curr_t<tspan(2)
        [XB, evals, h_next, redo] = explicit_RK_variable_step(rate_func_in,curr_t,XA,curr_h,BT_struct,p,error_desired);
        num_evals = num_evals + evals;
        % if redo == 1
        %     % curr_h = curr_h / 2;
        % else
        %     X_list = [X_list, XB];
        %     XA = XB;
        %     curr_t = curr_t + curr_h;            
        %     t_list = [t_list, curr_t];
        %     curr_h = h_next;
        %     h_list = [h_list, curr_h];
        % end

        if redo == 0
            X_list = [X_list, XB];
            XA = XB;
            % curr_t = curr_t + curr_h;
            curr_t = t_next;
            t_list = [t_list, curr_t];
            h_list = [h_list, curr_h];
        else
            fail = fail+1;
        end

        curr_h = h_next;

        t_next = curr_t + curr_h;
    
        if t_next > tspan(2)
            t_next = tspan(2);
            curr_h = tspan(2)-curr_t;
        end


        % if tspan(2) - curr_t < h_next
        %     h_next = tspan(2) - curr_t;
        % end
    end
    failure_rate = fail / (fail + length(X_list));
    h_avg = mean(h_list);
end