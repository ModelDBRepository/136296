-- Simple ODE solvers to use with Lua
--

module(..., package.seeall);

--TODO: 
-- * Runge-Kutta-Fehlberg with adaptive step size

-- One of the simplest: Runge-Kutta  4th order (non-stiff ODES)
function rk4(U0, t, h, s, obj)
   local k1, k2, k3, k4
   local U = U0+0
   local U1      -- make a copy

   k1 = obj:derivs(t, U)*h
   k2 = obj:derivs(t + .5*h, U + .5*k1) * h
   k3 = obj:derivs(t + .5*h, U + .5*k2) * h
   k4 = obj:derivs(t + h, U + k3) * h
   
   U1 = (k1 + 2*k2 + 2*k3 + k4) * (1/6)
   return U1, t+h, h
end


rkc_module = require 'rkc' -- Runge-Kutta-Chebyshev for mildly stiff ODEs
rkc = rkc_module.rkc    -- pointer to the default version of integrator
rkc_a = rkc_module.rkc_a


