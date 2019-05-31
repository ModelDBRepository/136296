module(...,package.seeall)

require "chlib"
require "el_ph"
require 'lib'

--- Refs:
--- SE-model: Swarz-Eikhof model [cited from Rattay and Aberham, 1993]
--- MRG: McIntyre, Richardson, Grill model [MRG, 2002]
--- MRGAxon: MRG, but original NEURON files


local global_variable_values = {
   K_o = 3.2,   -- mM  -- to check
   Na_i = 12,  -- mM
   V_a = -87,
}




local temp = 37
local global_parameters = {
   dNai_switch = 1, -- either 0 or 1
   dKo_switch = 1,
   rho_a = 0.07,        -- KOhm*cm, axonal resistivity
   T = temp + 273,  -- temperature
   celsius = temp,
   V_leak = -80, 
   c_m = 2, -- uF/cm^2, membrane capacitance per cm^2
   Na_o = 140, -- mM
   Na_i = 9.5,  -- mM
   K_i = 155,  -- mM   -- to check
   K_o = 3.2,  -- mM  -- to check
   g_s = 1.0,  -- mS, transmyelin conductivity, per lamella membrane
   --ko_sig = lib.logistic(200, 12, 4.2, 1), -- (interesting)
   --ko_sig = lib.logistic(137, 2.5, 6, 1), -- (normal?)
   --ko_sig = lib.logistic(90, 2.5, 6.6, 1), --[poster]
   ko_sig = lib.logistic(4, 3.5, 7.0, 1), --[paper sc1!]
   --ko_sig = lib.logistic(4, 3.5, 7.0, 1),
   --ko_sig = lib.hill(4, 25, 7.0, 1),
   --ko_sig = lib.logistic(4, 4.0, 6.5, 1),
   nak_pump = {
      --jnamax = 0.15e-6, -- mmole/(cm^2*s) -- in kidney
      jnamax = 0.03e-6, -- mmole/(cm^2*s) -- 
      a_na = 0.5, -- mM
      a_k = 0.1, -- mM
      --a_k = 5.4, -- mM
      b_na = 2.4e-2,
      b_k = 5.4e-3,
   },
   nak_pump_a = {
      jnamax =0.07e-6, --mmole/(cm^2*s) 
      a_na = 9,
      b_na = 1e-2,
      a_k = 0.1,
      b_k = 1e-3,
   },
   nak_pump_g = {
      jnamax =0.07e-6, --mmole/(cm^2*s) 
      a_na = 0.1,
      b_na = 1e-2,
      a_k = 4.5,
      b_k = 1e-3,
   }
}

segments={} -- compartment descriptions


segments.SC_mine = {
   -- Geometry-based (2nd version) model
   -- Some parameters are from HANA-2005
   type = "SC_geom",
   pars = {
      stimulated = true,
      node = {
	 c_m = 2, --uF/cm^2
	 l = 1e-4, d = 3.5e-4,
	 h = 9e-4, -- cm, height of the restricted space (was: 8)
	 g_lk = 2, V_lk = -80,
	 wavy = 9, -- account for 'wavyness' of glial membrane in node
	 pump_a = global_parameters.nak_pump_a,
	 tau_passive = 2e4, -- passive K exchange in the node char. time (ms)
	 currents = {
	    iNat = {
	       ion = "Na",
	       gp = {m_n  = chlib.fsc.m.srb, h_n = chlib.fsc.h.srb},
	       g  = 2800},
	    iNap = {
	       ion = "Na",
	       gp = {p_n = chlib.psc.p.hana},
	       g = 24},
	    iKs = {
	       ion = "K",
	       gp = {s_n = chlib.spc.s.my},
	       --g = 112,
	       g = 90}}},
      int = { -- internode
	 c_m = 1, c_s = 0.1, l= 0.129, d = 8.8e-4,
	 lpara = 80e-4,
	 D = 10e-4, h = 4e-7, g_lk = 4.2e-3, --g_s = 1.0,
	 Nl = 140,
	 R_il = 50e3, -- internodal leakage, not yet based on
	 --geometry!!!
	 kparaw = 0.58, -- доля каналов в параноде.
	 naparaw = 0.06,
	 tau_passive_K = 2e3, -- passive K exchange char. time (ms)
	 tau_passive_Na = 5e3, -- passive Na exchange char. time (ms)
	 tau_sw = 1000, -- tau for sw changes
	 g_Kglia = 2, -- total V-independent K conductance in Glion
	 g_lk_glia = 0.01,
	 V_lk_glia = -85,
	 currents = {
	    iKs = {
	       ion = "K",
	       gp = {s_i = chlib.spc.s.my},
	       g = 0.23},
	    iKf = {
	       ion = "K",
	       gp = {n_i = chlib.fpc.n.hana},
	       g = 3}}},
   },
   state = {V_n = -80, V_int = -80, V_g = -90,
	    m_n = 0.2, h_n = 0.2, p_n = 0.2,  -- nonsense values
	    s_n = 0.1, s_i = 0.2, n_i = 0.2,
	    Ko_n = 3.2, Ko_p = 3.2, Ko_i = 3.2, Nai_n = 9, Nai_p = 9,
	    Nai_i = 9,
	    sw = 1.0,
	 },
   virtuals = {ipump_p = 0, iKfp = 0, iKfi = 0}
}

segments.SC_mine.pars.node.pump_a.jnamax = global_parameters.nak_pump_a.jnamax*1.0

segments.Hana = { -- Based on HAA-2005 paper
   type = "SC_simple",
   pars = {
      stimulated = true,
      C_n  = 0.22e-6, -- uF
      C_my = 0.17e-6,
      C_i =  379e-6,
      G_lk_i =  1.7e-6,
      R_il = 41e3,
      node_currents = {
	 iNat = {
	    ion = "Na",
	    gp = {m_n  = chlib.fsc.m.hana, h_n = chlib.fsc.h.hana},
	    G  = 276e-6},
	 iNap = {
	    ion = "Na",
	    gp = {p_n = chlib.psc.p.hana},
	    --G = 7.1e-6},
	    G = 2.5e-5},
	 iKs = {
	    ion = "K",
	    gp = {s_n = chlib.spc.s.hana},
	    G = 17.4e-6},
	 iKf = {
	    ion = "K",
	    gp = {n_n = chlib.fpc.n.hana},
	    G = 4.1e-6}
      },
      internode_currents = {
	 iKs = {
	    ion = "K",
	    gp = {s_i = chlib.spc.s.hana},
	    G = 87.1e-6}
      },
   },
   state = {
      V_n = -86, V_int = -86,
      h_n = 0.2, m_n = 0.2, p_n = 0.2,s_n = 0.2,
      n_n = 0.2, s_i = 0.2,
   },
   virtuals = {ipump=0},
}

segments.Hana2 = { 
   -- Based on HANA-2005 paper, other currents, closer
   -- to my own model
   -- As well as Hana1 is a strangely hyperpolarized
   type = "SC_simple",
   pars = {
      stimulated = true,
      C_n  = 0.22e-6, -- uF
      C_my = 0.17e-6,
      C_i =  379e-6,
      --
      G_lk_i =  1.7e-6,
      --
      R_il = 41e3,
      --R_il  = 120e3,
      --
      node_currents = {
	 iNat = {
	    ion = "Na",
	    gp = {m_n  = chlib.fsc.m.srb, h_n = chlib.fsc.h.srb},
	    G  = 276e-6},
	 iNap = {
	    ion = "Na",
	    gp = {p_n = chlib.psc.p.my},
	    G = 7.1e-6},
	 iKs = {
	    ion = "K",
	    gp = {s_n = chlib.spc.s.my},
	    G = 17.4e-6},
	 iKf = {
	    ion = "K",
	    gp = {n_n = chlib.fpc.n.hana},
	    G = 0.1e-6}
      },
      internode_currents = {
	 iKs = {
	    ion = "K",
	    gp = {s_i = chlib.spc.s.hana},
	    G = 87.1e-6},
	 iKf = {
	    ion = "K",
	    gp = {n_i = chlib.fpc.n.hana},
	    G = 100e-6}}

   },
   state = {
      V_n = -86, V_int = -86,
      h_n = 0.2, m_n = 0.2, p_n = 0.2,s_n = 0.2,
      n_n = 0.2, s_i = 0.2, n_i = 0.2,
   },
   virtuals = {ipump=0},
}




sc1 = {tag='sc1', description = "A space-clamped model",
	     rule = {'SC_mine'}}

hana1 = {tag='hana1', description = "A space-clamped model (Hana)",
	       rule = {'Hana'}}
hana2 = {tag='hana2', description = "A space-clamped model (Hana)",
	       rule = {'Hana2'}}



--for k, gp in pairs(segments.Node.pars.gating_particles) do
   --io.stderr:write(string.format('k is %s\n', k))
--   if segments.Node.state[k] then 
--      segments.Node.state[k] = gp.inf(global_variable_values.V_a)
--   end
--end


for _,segment in pairs(segments) do
   for key,val in pairs(global_variable_values) do
      if segment.state[key] then -- only for those with such a variable
	 segment.state[key] = val
      end
   end
end

for _,segment in pairs(segments) do 
   for key,val in pairs(global_parameters) do
      if not segment.pars[key] then
	 segment.pars[key] = val
      end
   end
end

numerics = {
   --stepper = 'rkc'
   stepper = 'rkc_a'
}

