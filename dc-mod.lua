--#!/usr/bin/env luajit 

--[[
-- Trial variant of the double-cable model.
-- Will have model nodes, standard internodes and paranodes
--
-- version with explicit MYSA. More stiff
--
--
--]]

--require "profiler"

module(...,package.seeall)

require 'lib'
require 'el_ph' -- auxilary routines and electrophysiology
require 'ode'   -- ODE solvers

function ficks_flux(C_neigh, C_this, D, L)
   -- expects C in mM, D in cm^2/s, L in cm
   -- returns J in mM*cm/ms
   -- multiply by 1e-3 to convert seconds to milliseconds
   return 1e-3 * D * (C_neigh - C_this)/L 
end

function g_link(g_neigh, g_this)
   return 2/(1/g_this + (g_neigh and 1/g_neigh or 0))
end

format = string.format

---------------------------------------------------

--function ghk_node(dv,s,p)
   ---- ghk Na current
   --local iNa = p.P_Na * el_ph.ghk(p.Na_i, p.Na_o, s.V_a, p.T)*s.m1^3*s.h1 
   --local iNap = p.P_Nap * el_ph.ghk(p.Na_i, p.Na_o, s.V_a, p.T)*s.p^3
   
   -- ghk K  current (slow potassium current)
   --local iK = p.P_K * el_ph.ghk(p.K_i, p.K_o, s.V_a, p.T)*s.n*s.n         
   --local iKs = p.P_Ks * el_ph.ghk(p.K_i, p.K_o, s.V_a, p.T)*s.s
--end

function set_q10c(gating_particles, temp)
   for k,gp in pairs(gating_particles) do
      if not gp.q10c then 
	 gp.q10c = gp.q10^((temp - gp.tref)/10)
      end
   end
   return gating_particles
end

function tanhp1(x)
   return 1 + math.tanh(20*x)
end

function apply_current(t, istim)
   local ia,tx,L
   ia = 0
   for _,iappl in ipairs(istim) do
      local p = iappl
      if (t > p.start - 2*p.width)  and (t < p.stop + 2*p.width) then 
	 tx = (t - p.start) % p.int
	 ia = ia + p.amp*.0625 * tanhp1(tx) * tanhp1((p.width - tx)) * tanhp1((t  - p.start)) * tanhp1((p.stop-t))
      end
   end
   return ia
end

function update_gps(gps, V, s, dv, vid)
   for k,gp in pairs(gps) do
      dv[vid[k]] = 
	 gp.q10c*(gp.alpha(V)*(1 - s[vid[k]]) 
	    - gp.beta(V)*s[vid[k]])
   end
end


function hh_iv(T, z)
   return function (V, Ci, Co)
	     return (V - el_ph.nernst(Ci, Co, T, z))
	  end
end

function ghk_iv(T, z)
   return function (V, Ci, Co)
	     return el_ph.ghk(Ci, Co, V, T, z)
	  end
end


----------------------------
---  Compartment `class' ---
----------------------------

Compartment = lib.inheritFrom(nil)  -- table for a new class


function Compartment:new(phase,params,virtuals,tag)
   local new_inst = {} -- new instance
   setmetatable(new_inst, new_inst)
   new_inst.__index = new_inst
   new_inst.v0 = lib.deepcopy(phase)
   new_inst.p = lib.deepcopy(params)
   new_inst.virt = lib.deepcopy(virtuals)
   new_inst.memoized = {}
   self.base_tag = tag
   return new_inst
end


-- Compartment methods --
function Compartment:type()
    return "compartment"
end

--function Compartment:set_tag(base_tag)
--   self.tag = base_tag..'-'..self.type()
--end

function Compartment:set_right(rcomp) self.right = rcomp end -- set right link
function Compartment:set_left(lcomp) self.left = lcomp end   -- set left link
--function Compartment:apply_phase(phase) self.v = phase end   
function Compartment:apply_parameters(pars) self.p = pars end
function Compartment:erase_memoized() self.memoized = {} end


function Compartment:chain_link(addcomp)
   self:set_right(addcomp)
   addcomp:set_left(self)
end

function Compartment:get_var(key, phase)
   return (self.vid[key] and phase[self.vid[key]]) or self.virt[key]
end

function Compartment:set_var(key,val)
   if self.v0[key] then self.v0[key] = val 
   else self.virt[key] = val
   end
end


function Compartment:cc_dist(neigh)
-- center to center distance
   if neigh then
      return .5*(neigh.p.l + self.p.l)
   else
      return .5*self.p.l
   end
end

function Compartment:axial_current(neigh, phase)
   return ((neigh:get_Vi(phase) - self:get_Vi(phase))*
	g_link(neigh.p.G_a, self.p.G_a))
end

function Compartment:longit_current(neigh, phase, func)
   if neigh then
      return func(self, neigh, phase)
   else
      return 0
   end
end

function Compartment:longit_current2(phase, func)
   return (self:longit_current(self.left, phase,func) + 
	self:longit_current(self.right, phase,func))
end

-- Review dimensions: I have D in cm^2/s, Concentration in mmole/l, which is eqv. to mole/m^3, 
--                    time in ms, lenght in cm, flux in mmole/s
-- 
--                      10^-4 m^2     mole        1           1         mole        
--  (S/V)*D*J*dc/dx  = ----------- * ------- * -------- * -------- =  ---------- = 10^3 mM/ms
--                         s           m^3      10^-2 m    10^-2 m      s * m^3



function Compartment:ion_exchange_maker(neigh, ion_name, D, sect_name)
   if neigh then 
      local dist = self:cc_dist(neigh)
      local csect = self:mean_csect(neigh, sect_name)
      return function (phase)
		local flux -- mM*cm^3/ms
		flux = ficks_flux(phase[neigh.vid[ion_name]], 
				  phase[self.vid[ion_name]], 
				  D, dist)
		return flux * self:mean_csect(neigh, sect_name)
	     end
   end 
   return function(phase) return 0 end
end

function Compartment:min_csect(neigh, key)
   local sect = self.p[key]
   if neigh and neigh.p[key] and  neigh.p[key] < sect then
      sect = neigh.p[key]
   end
   return sect
end

function Compartment:mean_csect(neigh, key)
   local sect = self.p[key]
   if neigh and neigh.p[key] then
      sect = 0.5*(neigh.p[key] + self.p[key])
   end
   return sect
end


function Compartment:Na_exchange(phase)
   return (self.Na_exchange_left(phase) 
	+ self.Na_exchange_right(phase))/self.p.axonal_vol
end

function Compartment:K_exchange(phase)
   return (self.K_exchange_left(phase) 
	+ self.K_exchange_right(phase))/self.p.periax_vol
end

function Compartment:setup()
   --self.Iax_longit = self:Iaxial_maker()
   local p = self.p
   p.axonal_csect = .25 * math.pi * p.d^2
   p.periax_csect = el_ph.ring_area(p.d + 2*p.h, p.d) 

   if self.right then 
      local rp = self.right.p
      rp.periax_csect = el_ph.ring_area(rp.d + 2*rp.h, rp.d) 
      rp.axonal_csect = .25 * math.pi * rp.d^2
   end

   self.Na_exchange_left = 
      self:ion_exchange_maker(self.left,'Na_i',	el_ph.Daq.Na,
			      'axonal_csect')

   self.Na_exchange_right = 
      self:ion_exchange_maker(self.right,'Na_i', el_ph.Daq.Na, 'axonal_csect')

   self.K_exchange_left = 
      self:ion_exchange_maker(self.left,'K_o', el_ph.Daq.K, 'periax_csect')

   self.K_exchange_right = 
      self:ion_exchange_maker(self.right,'K_o', el_ph.Daq.K, 'periax_csect')
end


---------------------------------------


---------------------------------------
-- Space-clamped models
---------------------------------------

SC_simple = lib.inheritFrom(Compartment)

function SC_simple:type()
   return 'scsimp'
end

function set_q10c2(gp, temp)
   if not gp.q10c then
      gp.q10c = gp.q10^((temp - gp.tref)/10)
   end
end

function SC_simple:setup()
   local p = self.p

   for _,current in pairs(p.node_currents) do
      set_q10c(current.gp, p.celsius)
   end

   for _,current in pairs(p.internode_currents) do
      set_q10c(current.gp, p.celsius)
   end
end

function curr_g(current, s, vid)
   local g = current.G
   for key,gp in pairs(current.gp) do
      g = g*s[vid[key]]^gp.order
   end
   return g
end


function SC_simple:deriv(s,dv,t)
   local vid = self.vid
   local V_n = s[vid.V_n]
   local V_int = s[vid.V_int]

   local p = self.p

   local Erest = {
      Na = el_ph.nernst(p.Na_i, p.Na_o, p.T),
      K = el_ph.nernst(p.K_i, p.K_o, p.T)}
   
   local iex = apply_current(t, p.iappl)
   local ilk_i = p.G_lk_i * (V_int - Erest.Na)

   local vnvi = (V_n - V_int)/(p.R_il)--/s[vid.sw])
   local dvn = iex - vnvi

   for _,i in pairs(p.node_currents) do
      dvn = dvn - curr_g(i, s, vid)*(V_n - Erest[i.ion])
      update_gps(i.gp, V_n, s, dv, vid)
   end
   dvn = dvn/(p.C_n + p.C_my)
   
   local dvi = vnvi + dvn*p.C_my - ilk_i
   
   for _,i in pairs(p.internode_currents) do
      dvi = dvi - curr_g(i, s, vid)*(V_int - Erest[i.ion])
      update_gps(i.gp, V_int, s, dv, vid)
   end
   
   dv[vid.V_n] = dvn
   dv[vid.V_int] = dvi/p.C_i

   return dv
end

SC_geom = lib.inheritFrom(Compartment)

function SC_geom:type()
   return 'scgeom' end

function SC_geom:setup()
   local p = self.p

   p.int.G_il = 1/p.int.R_il

   p.node.asurf = math.pi * p.node.l * p.node.d -- axonal surface area

   p.int.asurf  = math.pi * p.int.l * p.int.d
   p.int.asurf_p  = math.pi * p.int.lpara * p.int.d
   p.int.gsurf = p.int.asurf*1.1 -- glial surface facing periaxonal space

   p.node.acsect = .25 * math.pi * p.node.d^2 -- axonal cross-section
   p.int.acsect  = .25 * math.pi * p.int.d^2

   p.node.periax_csect = 
      el_ph.ring_area(p.node.d + 2*p.node.h, p.node.d)

   p.node.gsurf = p.node.wavy*p.node.periax_csect

   p.int.periax_csect = 
      el_ph.ring_area(p.int.d + 2*p.int.h, p.int.d)

   --p.int.glia_csect = 
   --   el_ph.ring_area(p.int.d +  2*(p.int.h + p.int.hg), 
	--			      p.int.h + 2*p.int.h)

   p.node.ax_vol = p.node.l * p.node.acsect -- axonal volume
   p.int.ax_vol = p.int.l * p.int.acsect
   p.int.ax_vol_p = p.int.lpara * p.int.acsect
   --p.int.glial_vol = p.int.l * p.int.glia_csect

   p.node.periax_vol = p.node.l*p.node.periax_csect
   p.int.periax_vol  = p.int.l *p.int.periax_csect
   p.int.periax_vol_p  = p.int.lpara *p.int.periax_csect

   for _,i in pairs(p.node.currents) do
      i.G = i.g*p.node.asurf
   end
   
   for _,i in pairs(p.int.currents) do
      i.G = i.g*p.int.asurf
   end


   
   p.node.G_lk = p.int.g_lk * p.int.asurf  --leakage conductance
   p.int.G_lk  = p.int.g_lk * p.int.asurf  --leakage conductance
   p.int.G_lk_glia = p.int.g_lk_glia * p.int.gsurf
   p.int.G_Kglia = p.int.g_Kglia*p.int.gsurf
   p.int.C_glia = p.int.c_m*p.int.gsurf


   --p.int.G_s=el_ph.myelin_conductance(p.int.g_s, p.int.D, p.int.d, p.int.l, p.int.Nl)

   p.int.C_s=el_ph.myelin_capacitance(p.int.c_s, p.int.D, p.int.d, p.int.l, p.int.Nl)
   p.node.C_m = p.node.asurf*p.node.c_m  -- capacity of the axonal membrane
   p.int.C_m = p.int.asurf*p.int.c_m  -- capacity of the axonal membrane

   -- q10c!
   for _,current in pairs(p.node.currents) do
      set_q10c(current.gp, p.celsius)
   end

   for _,current in pairs(p.int.currents) do
      set_q10c(current.gp, p.celsius)
   end
end


function SC_geom:deriv (s, dv, t)
   local vid = self.vid
   local V_n = s[vid.V_n]
   local V_int = s[vid.V_int]

   local p = self.p
   
   local Erest_n = {
      Na = el_ph.nernst(s[vid.Nai_n], p.Na_o, p.T),
      K = el_ph.nernst(p.K_i, s[vid.Ko_n], p.T)}

   local Erest_p = {
      Na = el_ph.nernst(s[vid.Nai_p], p.Na_o, p.T),
      K = el_ph.nernst(p.K_i, s[vid.Ko_p], p.T)}

   local Erest_i = {
      Na = el_ph.nernst(s[vid.Nai_i], p.Na_o, p.T),
      K = el_ph.nernst(p.K_i, s[vid.Ko_i], p.T)}
   
   local iex = apply_current(t, p.iappl)
   local ilk_n = p.node.G_lk * (V_n - p.node.V_lk)
   
   --local vnvi = (V_n - V_int)/p.int.R_il
   local vnvi = (V_n - V_int)*p.int.G_il*s[vid.sw]

   local fluxes_n = {K = 0, Na = 0}
   local fluxes_i = {K = 0, Na = 0}
   local fluxes_p = {K = 0, Na = 0}

   local nakpf_a_n = el_ph.NaKpump_flux(p.K_i, s[vid.Ko_n], 
					s[vid.Nai_n], p.Na_o,
					p.node.pump_a)
   local nakpf_g_n = el_ph.NaKpump_flux(p.K_i, s[vid.Ko_n], 
					p.Na_i, p.Na_o,
					p.nak_pump_g)
   local nakpf_a_i = el_ph.NaKpump_flux(p.K_i, s[vid.Ko_i], 
					s[vid.Nai_i], p.Na_o,
					p.nak_pump_a)
   local nakpf_a_p = el_ph.NaKpump_flux(p.K_i, s[vid.Ko_p], 
					s[vid.Nai_p], p.Na_o,
					p.nak_pump_a)

   local nakpf_g_i = el_ph.NaKpump_flux(p.K_i, s[vid.Ko_i], 
					p.Na_i, p.Na_o,
					p.nak_pump_g)
   local nakpf_g_p = el_ph.NaKpump_flux(p.K_i, s[vid.Ko_p], 
					p.Na_i, p.Na_o,
					p.nak_pump_g)
   
   local ipump_n = el_ph.flux2i(nakpf_a_n[1])*p.node.asurf/3
   local ipump_i = el_ph.flux2i(nakpf_a_i[1])*p.int.asurf/3
   local ipump_p = el_ph.flux2i(nakpf_a_p[1])*p.int.asurf_p/3
   local ipump_g = el_ph.flux2i(nakpf_g_i[1])*p.int.gsurf/3

   exch_Na_np = 2*ficks_flux(s[vid.Nai_p], s[vid.Nai_n], 
			   el_ph.Daq.Na,
			   0.5*p.int.lpara)*p.node.acsect

   exch_Na_pi = 2*ficks_flux(s[vid.Nai_i], s[vid.Nai_p], 
			   el_ph.Daq.Na, 0.25*p.int.l)*p.int.acsect

   exch_K_np = 2*ficks_flux(s[vid.Ko_p], s[vid.Ko_n], 
			  el_ph.Daq.K, 0.5*p.int.lpara)*p.int.periax_csect

   exch_K_pi = 2*ficks_flux(s[vid.Ko_i], s[vid.Ko_p], 
			  el_ph.Daq.K, 0.25*p.int.l)*p.int.periax_csect


   fluxes_n.K = (nakpf_a_n[2]*p.node.asurf +
		 nakpf_g_n[2]*p.node.gsurf + 
		 exch_K_np)

   fluxes_p.K = (nakpf_a_i[2]*p.int.asurf_p + 
		 nakpf_g_i[2]*p.int.asurf_p*1.1 +
		 exch_K_pi - exch_K_np)

   fluxes_i.K = (nakpf_a_i[2]*p.int.asurf + 
		 nakpf_g_i[2]*p.int.gsurf -
		 exch_K_pi)

   fluxes_n.Na  = nakpf_a_n[1]*p.node.asurf - exch_Na_np
   fluxes_p.Na  = nakpf_a_i[1]*p.int.asurf_p + exch_Na_np - exch_Na_pi
   fluxes_i.Na  = nakpf_a_i[1]*p.int.asurf + exch_Na_pi

   local dvn = iex - (vnvi + ilk_n + ipump_n)
   
   local x
   for _,i in pairs(p.node.currents) do
      x = curr_g(i, s, vid)*(V_n - Erest_n[i.ion])
      fluxes_n[i.ion] =  fluxes_n[i.ion] + el_ph.i2flux(x)
      dvn = dvn - x
      update_gps(i.gp, V_n, s, dv, vid)
   end
   
   dvn = dvn/(p.node.C_m + p.int.C_s)
   
   
   local ilk_i = (1-p.int.naparaw)*p.int.G_lk * (V_int - Erest_i.Na)
   local ilk_p = p.int.naparaw*p.int.G_lk * (V_int - Erest_p.Na)

   local dvi = vnvi + dvn*p.int.C_s - (ilk_i + ilk_p + ipump_p + ipump_i)

   fluxes_i.Na = fluxes_i.Na + el_ph.i2flux(ilk_i)
   fluxes_p.Na = fluxes_p.Na + el_ph.i2flux(ilk_p)
  
   local x,xi,xp
   for iname,i in pairs(p.int.currents) do
      x = curr_g(i, s, vid)
      local v = p.int.kparaw*Erest_p[i.ion] + (1-p.int.kparaw)*Erest_i[i.ion] 
      xi = (1-p.int.kparaw)*(V_int - Erest_i[i.ion])*x
      xp = p.int.kparaw * (V_int - Erest_p[i.ion])*x
      if iname == "iKf" then
	 self.virt.iKfp = xp
	 self.virt.iKfi = xi
      end
      fluxes_i[i.ion] =  fluxes_i[i.ion] + el_ph.i2flux(xi)
      fluxes_p[i.ion] =  fluxes_p[i.ion] + el_ph.i2flux(xp)
      update_gps(i.gp, V_int, s, dv, vid)
      dvi = dvi - xi - xp
   end

   self.virt.ipump_p = ipump_p

   local iKg = p.int.G_Kglia*(s[vid.V_g] - Erest_i.K) 
   dv[vid.V_g] = -(iKg + 
		   p.int.G_lk_glia*(s[vid.V_g] - p.int.V_lk_glia)
		+ ipump_g)/p.int.C_m

   fluxes_i.K = fluxes_i.K + el_ph.i2flux(iKg)
   
   dv[vid.V_n] = dvn
   dv[vid.V_int] = dvi/p.int.C_m

   dv[vid.Nai_n] = -fluxes_n.Na/p.node.ax_vol   

   dv[vid.Nai_p] = -fluxes_p.Na/p.int.ax_vol_p + 
      (p.Na_i - s[vid.Nai_p])/self.p.int.tau_passive_Na

   dv[vid.Nai_i] = -fluxes_i.Na/p.int.ax_vol + 
      (p.Na_i - s[vid.Nai_i])/self.p.int.tau_passive_Na

   dv[vid.Ko_n] = fluxes_n.K/p.node.periax_vol + 
      (p.K_o - s[vid.Ko_n])/self.p.node.tau_passive
  
   dv[vid.Ko_p] = fluxes_p.K/p.int.periax_vol_p +
      (p.K_o - s[vid.Ko_p])/self.p.int.tau_passive_K

   dv[vid.Ko_i] = fluxes_i.K/p.int.periax_vol + 
      (p.K_o - s[vid.Ko_i])/self.p.int.tau_passive_K

   dv[vid.sw] = ((p.ko_sig(s[vid.Ko_p]) - s[vid.sw]) /  p.int.tau_sw)

   return dv
end




--------------------------

-----------------------------
---  Nerve `class' ----
-----------------------------

Nerve = {} 

function Nerve:parameter_setup(obj)
    local surface,volume,cross_section
    local p
    for i,comp in ipairs(obj.layout) do
        comp:setup()
    end
end

function Nerve:set_state(obj) 
   obj.phase0 = setmetatable({}, lib.vecmt)
   obj.derivatives = setmetatable({}, lib.vecmt)

    -- create a sorted list of variable names
    obj.sorted_keys = {} -- with virtuals
    obj.sorted_vkeys = {} -- phase only
    for i, comp in ipairs(obj.layout) do
       local x = {}
       local xv = {} -- without virtuals (phase only)
       for index,value in pairs(comp.v0) do
	  table.insert(x,index)
	  table.insert(xv, index)
       end
       if comp.virt then
	  --io.stderr:write(string.format("\n comp.virt: \n"))
	  for index,value in pairs(comp.virt) do
	     table.insert(x,index)
	  end
       end
       table.sort(x)
       table.sort(xv)
       table.insert(obj.sorted_keys,x)
       table.insert(obj.sorted_vkeys,xv)
    end

    -- setup variable indices for each compartment
    local j = 1
    for z,comp in ipairs(obj.layout) do
       comp.vid = {}
       for _, key in ipairs(obj.sorted_vkeys[z]) do
	  comp.vid[key] = j
	  obj.phase0[j] = comp.v0[key]
	  j = j + 1
       end
    end
 end


function Nerve:erase_memoized()
   for z, comp in ipairs(self.layout) do 
      comp:erase_memoized() 
   end
end

function Nerve:derivs(t, phase) -- calculates derivatives for a state
   --local dvs = setmetatable({}, lib.stmt)
--   print("Nerve.derivs: phase size = ", #phase)

   self:erase_memoized()

   local dvs = self.derivatives

   for n,comp in ipairs(self.layout) do -- now calculate derivatives
      dvs = comp:deriv(phase, dvs, t)
   end

   return dvs
end

function Nerve:run(t0, h, t1, stf)    
    -- numerically integrate from t0 to t1,
    -- using step h and stepper function stf
   local st,ks,n
   n = 0
   local save_step = 5e-3
   local max_h = 2
    
    ks = save_step > h and  math.floor(save_step/h) or 1
    local kp = ks*10

    local stepper = ode[stf]
    
    local t = t0
    local tnext

    -- phase we start from:
    local U = self.phase0 + 0
--    print("Nerve.run: U[1] = ", U[1])
    --local U1 = U + 0
    
    while t <= t1 do
       U, tnext, h = stepper(U, t, h, _s_, self)
       if math.mod(n,ks) == 0 then
	  if math.mod(n,kp) == 0 then
	     io.stderr:write(string.format('\rmodel time: %3.3f ms',t))
	  end
	  io.stdout:write(self:pr_state(U, t))
	  --io.stdout:write(self:pp_state(t))
	  io.stdout:flush()
       end
       n = n+1
       if h > max_h then 
	  h = max_h
	  tnext = t+h
       end
       t = tnext
    end
    --io.stderr:write('\n')
    --io.stdout:write('\n')
    self.phase0 = U + 0
    return t, U
end

function Nerve:pr_state(phase, t)
   local x = string.format('%3.5e',t)
   for z,comp in ipairs(self.layout) do
      for i,key in ipairs(self.sorted_keys[z]) do
	 x = x..string.format('\t%3.5e', comp:get_var(key, phase))
      end
   end
   return (x..'\n')
end

function Nerve:luaprint_state(phase, t)
   --prints state in a lua-loadable format
   --x = string.format("{time = %3.5e,\n" ,time)
   x = "return {"
   for z, comp in ipairs(self.layout) do
      x = x..string.format('{')
      for i,key in ipairs(self.sorted_keys[z]) do
	 x = x..string.format(' %s = %3.5e,', key, comp:get_var(key, phase))
      end
      x = x.."},\n"
   end
   return x.."}\n"
end

function Nerve:pythonprint_state(phase, t)
   --prints state in a python-loadable format
   --x = string.format("{time = %3.5e,\n" ,time)
   x = "["
   for z, comp in ipairs(self.layout) do
      x = x..string.format('{')
      for i,key in ipairs(self.sorted_keys[z]) do
	 x = x..string.format(" '%s': %3.5e,", key, comp:get_var(key, phase))
      end
      x = x.."},\n"
   end
   return x.."]\n"
end


function Nerve:pp_state(phase, t) -- pretty print state
   local x = string.format("(('time . %3.5e)", t)
   local j = 1
   for z,comp in ipairs(self.layout) do
      x = x..string.format("\n (%s-%d ", self.layout_rule.tag, z)
      for i,key in ipairs(self.sorted_keys[z]) do
	 x = x..string.format(" (%s-%s-%d %3.5e)", 
			      comp:type(), key, j, comp:get_var(key, phase))
	 j = j+1
      end
      x = x..")"
   end
   return (x..'))\n')
end

function Nerve:print_parameters()
   local x = ""
   for z,comp in ipairs(self.layout) do
      for key,value in pairs(comp.p) do
	 if type(value) == "number" then
	    x = x..string.format('%s-%s = \t%3.5e\n', comp:type(), key, value)
	 end
      end
   end
   return (x..'\n')
end

function Nerve:gp_print()
   local x = ''
   local j = 2
   local curr_type
   local types = {}
   for z,comp in ipairs(self.layout) do
      curr_type = comp:type()
      if not types[curr_type] then types[curr_type] = 1 end
      for i,key in ipairs(self.sorted_keys[z]) do
	 x = x..string.format("%s_%s%d_%s = %d\n", 
			      self.layout_rule.tag,
			      comp:type(), types[curr_type], key, j)
	 
	 j = j + 1
      end
      types[curr_type] = types[curr_type] + 1      
   end
   return x..'\n'
end
   
function Nerve:luaload_state(fname)
   local f, sts, x, j
   local state
   f = assert(loadfile(fname))
   state = f()
   j = 0
   for z, comp in ipairs(self.layout) do
      for key, val in pairs(state[z]) do
	 comp:set_var(key,val)
      end
   end
   Nerve:set_state(self)
end

function Nerve:construct_from_rule(rule, comps)
   local layout = {}
   local prev
   local node_count = 0
   
   local i =1
   for _,name in ipairs(rule.rule) do
      --if known[name] then
      local comp_type  = comps[name].type and comps[name].type or name
      local comp_class = getfenv()[comp_type]
      layout[i] = comp_class:new(comps[name].state,
				 comps[name].pars, 
				 comps[name].virtuals, '-')
      if i > 1 then
	 layout[i-1]:chain_link(layout[i])
      end
      if layout[i]:type() == 'node' then
	 node_count = node_count + 1
	 layout[i].p.number = node_count
      end
      i = i+1
   end
   return layout
end

function Nerve:new(compartments,layout_rule) -- class constructor
    local newn = {} -- new nerve instance
    setmetatable(newn,newn)
    newn.__index = self
    
    newn.layout_rule = layout_rule
    newn.layout = Nerve:construct_from_rule(layout_rule, compartments, newn.known)
    Nerve:set_state(newn)
    Nerve:parameter_setup(newn)
    return newn
end

