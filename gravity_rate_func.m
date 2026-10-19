%Rate function describing Newton's law of gravitation
%INPUTS:
%t: the current time
%V: the vector of the position and velocity of the planet
% V = [x_p; y_p; dxdt_p; dydt_p]
%orbit_params: a struct describing the system parameters
% orbit_params.m_sun is the mass of the sun
% orbit_params.m_planet is the mass of the planet
% orbit_params.G is the gravitational constant
%OUTPUTS:
%dVdt: a column vector describing the time derivative of V:
% dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p]
function dVdt = gravity_rate_func(t,V,orbit_params)

    m1 = orbit_params.m_sun;
    m2 = orbit_params.m_planet;
    g = orbit_params.G;
    
    d2xdt2_p = (-m1*g / (V(1)^2 + V(2)^2)^(3/2)) * V(1);
    d2ydt2_p = (-m1*g / (V(1)^2 + V(2)^2)^(3/2)) * V(2);

    dVdt = [V(3);V(4);d2xdt2_p;d2ydt2_p];

end