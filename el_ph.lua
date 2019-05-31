-- Electrophysiology and common stuff for nerve models --

module(..., package.seeall);

-- useful constants
R = 8314.4        -- Universal gaz constant, mJ / mole*K
F = 96485.34      -- Faraday constant Coulomb / mole

-- Diffusion coefficients in aqueous solutions --
-- [cm^2/s]
-- (Hille, B., Ionic Channels of Excitable Membranes)
Daq = {K = 1.96e-5, Na = 1.33e-5, H = 9.31e-5, Ca = 0.79e-5}



local gauss_mem = {}
setmetatable(gauss_mem, {__mode = "v"})  -- make gauss_mem weak

function gaussian(t,p)
   --local key = t..'-'..p.t..'-'..p.w..'-'..p.a
   --if gauss_mem[key] then  
   --   return gauss_mem[key]
   --else
      local x,g,z
      z = 10
      if (t - p.t) * (t - p.t) < z * p.w then
	 g = (t - p.t)^2
	 g = -.5*g / p.w
	 x =  p.a * math.exp(g)
      else x = 0 end
      gauss_mem[key] = x
      return x
   --end
end

-- Some biophysics --

function NaKpump_flux(Ki, Ko, Nai, Nao, p)
   --                                mmole
   -- NOTE: should return fluxes in --------
   --                                cm^2 s  
   -- version 1.0 -- no optimization
   
   local jna, jk
   
   local Km_na = p.a_na + p.b_na * Ki
   local Km_k = p.a_k + p.b_k * Nao
   
   local mm_na = 1 / (1 + Km_na/Nai)^3 
   local mm_k =  1 / (1 + Km_k/Ko)^2

   jna = p.jnamax*mm_na*mm_k --+ arg.p.jrest
   jk = -jna / 1.5
   
   return {jna, jk}
end


-- Basic circuits --

function i2flux(I)
   return I/(1e3*F)
end

function flux2i(flux)
   return flux*1e3*F
end

function reciprocal(x)  return 1/x end

function rsumr(values)
   -- reciprocal of the summed reciprocals --
   local res = 0
   for i,g in ipairs(values) do
   --for i = 1,#values do 
      res = res + 1/g
   end
   return 1/res
end

function sum(values)
   local res = 0
   for _,v in ipairs(values) do res = res+v end
   return res
end

function series_resistances(rrr)
   return sum(rrr)
end

function series_conductances(ggg)
   return rsumr(ggg)
end

function parallel_resistances(rrr)
   return rsumr(rrr)
end

function parallel_conductances(ggg)
   return sum(ggg)
end

-- Electrochemistry functions --

function nernst(ci, co, T, z)
   if not z then z = 1 end -- default charge
   return R*T*math.log(co / ci)/(z*F)
end


-- GHK current equation --
function ghk(ci, co, V, T, z) 
   if not z then z = 1 end -- default charge
   local x = V*F*z / (R*T)
   if math.abs(V) < 1e-7 then
      return (ci - co)*R*T/z
   else
      local y = math.exp(x)
      return x*F*(co-ci*y)/(1 - y)
   end
end


local ghk_mem = {}
setmetatable(ghk_mem, {__mode = "v"})  -- make ghk_mem weak


function make_subkey(vector)
   local key = "-"
   for i,val in ipairs(vector) do
      key = key..'val'..'-'
   end
   return key
end

function ghk_rate_const(p,V)
      local x = V + p[1]
      local rc = (p[2]*x + p[3]) / (1 + p[4]*math.exp(x/p[5]))
      return rc 
end

function myelin_area(D, L)
    --return math.pi*(D+d)*(L + .5*(D-d))
    --return math.pi*(D*L + .5*(D^2 - d^2))
    return math.pi*D*L
end

function ring_area(D, d)
   return .25*math.pi*(D^2 - d^2)
end

function periaxonal_conductance(rho, d, h, l)
   return ring_area(d + 2*h, d)/ (rho*l)
end

function myelin_conductance(gs, D, d, L, N)
   -- two membranes per lamella
   --return .5 * math.pi * gs * L * D / N
   return myelin_area(D, L)*gs/(2*N)
end

function myelin_capacitance(cs, D, d,L, N)
   -- two membranes per lamella
   --return .5 * math.pi * cs * L * D / N 
   return cs *myelin_area(D, L) /(2*N) 
end


