-- Ion currents library --

module(..., package.seeall);

require 'lib'

function particle_tau (alpha, beta)
   return function (V) return 1 / (alpha(V) + beta(V)) end
end

function particle_inf (alpha, beta)
   return function (V) return 1 / (1 + beta(V) / alpha(V)) end
end

function gating_rater1(a,b,c)
-- returns a rate function, like \alpha or \beta
   return function (V)
	     local x = V - b
	     if math.abs(x) < 1e-7 then 
		return -a*c
	     else 
		return a*x/(1 - math.exp(x/c))
	     end
	  end
end


function gating_rater2(a,b,c)
   return function (V)
	     return a/(1 + math.exp((V-b)/c))
	  end
end


------------------------------------------------------------
---- Fast K current activation -----------------------------
------------------------------------------------------------

fpc = {n={}}

-- MRG:
fpc.n.mrg = {order=4, q10 = 3.0, tref = 20}
fpc.n.mrg.alpha = gating_rater1(7.97e-3, -83.2, -1.1)
fpc.n.mrg.beta = gating_rater1(-1.41e-2, -66, 10.5)


fpc.n.srb = { --fast K activation [Swarz, Reid, Bostock]
   q10 = 3.0,
   tref = 20,
   order = 4,
   alpha = gating_rater1(7.98e-3, -93.2, -1.1),
   beta = gating_rater1(-1.42e-2, -76.0, 10.5),
}

fpc.n.hana = lib.deepcopy(fpc.n.srb)

fpc.n.my = { --fast K activation [me]
   q10 = 3.0,
   tref = 20,
   order = 4,
   alpha = gating_rater1(7.98e-3, -87.2, -1.1),
   beta = gating_rater1(-1.42e-2, -70, 10.5),
}


for _,gv in pairs(fpc.n) do
   gv.inf = particle_inf(gv.alpha, gv.beta)
   gv.tau = particle_tau(gv.alpha, gv.beta)
end

------------------------------------------------------------
---- Fast sodium current activation ------------------------
------------------------------------------------------------
fsc = {m = {},h = {}}

fsc.m.se = { -- fast Na current activation
   -- from Scwharz-Eikhof model
   order = 3,
   q10 = 2.2,
   tref = 37,
   alpha = gating_rater1(1.87, 25.4-80, -6.06),
   beta = gating_rater1(-3.97, 21-80, 9.41),
}

fsc.m.mrg = { --fast Na current activation, [McIntyre :MRGaxon/Axnode.mod]
   q10 = 2.2,
   tref = 20,
   order = 3,
   alpha = gating_rater1(1.86, -20.4, -10.3),
   beta = gating_rater1(-0.086, -25.7, 9.16),
}

fsc.m.srb = { -- fast Na current activation, [Swarz, Reid, Bostock]
   q10 = 2.2,
   tref = 20,
   order = 3,
   alpha = gating_rater1(1.86, -18.4, -10.3),
   beta = gating_rater1(-0.086, -22.7, 9.16)
}

fsc.m.hana = lib.deepcopy(fsc.m.srb)

------------------------------------------------------------


------------------------------------------------------------
---- Fast sodium current inactivation ----------------------
------------------------------------------------------------
fsc.h.mrg = { --fast Na current inactivation [McIntyre :MRGaxon/Axnode.mod]
   q10 = 2.9,
   tref = 20,
   order = 1,
   alpha = gating_rater1(-0.062, -114, 11),
   beta = gating_rater2(2.3, -31.8, -13.4),
}

fsc.h.srb = { --fast Na current inactivation [Swarz, Reid, Bostock]
   q10 = 2.9,
   tref = 20,
   order = 1,
   alpha = gating_rater1(-0.0336, -111, 11),
   beta = gating_rater2(2.3, -28.8, -13.4),
}

fsc.h.hana = lib.deepcopy(fsc.h.srb)

fsc.h.se = { -- fast Na current inactivation
   -- alpha and beta are taken from Scwharz-Eikhof model
   q10 = 2.9,
   order = 1,
   tref = 37,
   alpha = gating_rater1(-0.55, -(80+27.74), 9.06),
   beta = gating_rater2(22.6, 56-80, -12.5)
}

----------------------------------------
for _,gv in pairs(fsc) do
   for _,gvv in pairs(gv) do
      gvv.inf = particle_inf(gvv.alpha, gvv.beta)
      gvv.tau = particle_tau(gvv.alpha, gvv.beta)
   end
end

------------------------------------------------------------


------------------------------------------------------------
---- Persitent sodum currents ------------------------------
------------------------------------------------------------

psc = {p={}}
psc.p.mrg = { -- persistent sodium current activation
   -- [McIntyre etal 2002 :MRGaxon/Axnode.mod]
   q10 = 2.2,
   tref = 20,
   order = 3,
   alpha = gating_rater1(0.01, -27, -10.2),
   beta = gating_rater1(-2.5e-4, -34, 10),
}

psc.p.fh = { -- persistent Na current activation
   -- [FH, cit. from Rattay-Aberham]
   q10 = 2.2,
   tref = 20,
   order = 2,
   alpha = gating_rater1(6e-3, 40-80, -10),
   beta = gating_rater1(-9e-2, -(80+25), 20),
}

psc.p.hana = { -- persistent Na current activation
   -- [SRB]
   -- [Hennings, Arendt-Nielsen, Andersen, 2005]
   q10 = 2.2,
   tref = 20,
   order = 3,
   alpha = gating_rater1(0.93, -38.4, -10.3),
   beta = gating_rater1(-0.043, -42.7, 9.16),
}

psc.p.my = { -- persistent Na current activation
   -- [Hennings, Arendt-Nielsen, Andersen, 2005] -> me
   q10 = 2.2,
   tref = 20,
   order = 3,
   alpha = gating_rater1(0.1, -38.4, -10.3),
   beta = gating_rater1(-0.02, -30.7, 9.16),
}

for _,gv in pairs(psc.p) do
   gv.inf = particle_inf(gv.alpha, gv.beta)
   gv.tau = particle_tau(gv.alpha, gv.beta)
end

------------------------------------------------------------
---- Slow potassium currents -------------------------------
------------------------------------------------------------

spc = {s={}}

spc.s.mrg = { -- Ks current [McIntyre etal 2002]
   q10 = 3.0,
   tref = 36,
   order = 1,
   alpha = gating_rater2(0.3, -53, -5), -- strange
   beta = gating_rater2(0.03, -90, -1), -- strange
}

spc.s.srb = { -- Ks current, [Swarz, Reid, Bostock, 1995]
   q10 = 3.0,
   tref = 20,
   order = 1,
   alpha = gating_rater1(1.22e-3, -12.5, -23.6),
   beta = gating_rater1(-7.39e-4, -80.1, 21.8),
}

spc.s.hana = { -- Ks current, [Hennings, Arendt-Nielsen, Andersen, 2005]
   q10 = 3.0,
   tref = 20,
   order = 1,
   alpha = gating_rater1(1.22e-2, -12.5, -16.9),
   beta = gating_rater1(-7.36e-4, -80.1, 12.6),
}


spc.s.my = { -- Ks current, 
   -- [Hennings, Arendt-Nielsen, Andersen, 2005] -> me
   q10 = 3.0,
   tref = 20,
   order = 1,
   alpha = gating_rater1(1.2e-3, -19, -18),
   beta = gating_rater1(-7.4e-4, -82, 16),
}


spc.s.se = { -- (slow?) potassium current from Swarz-Eikhof model
   q10 = 3.0,
   tref = 37,
   order = 2,
   alpha = gating_rater1(0.13, 35-80, -10),
   beta = gating_rater1(-0.32, 10-80, 10),
}

for _,gv in pairs(spc.s) do
   gv.inf = particle_inf(gv.alpha, gv.beta)
   gv.tau = particle_tau(gv.alpha, gv.beta)
end
------------------------------------------------------------
-- References:
--
-- Hennings, Arendt-Nielsen, Andersen. Breakdown of accomodation in
-- nerve: a possible role for persistent sodium current,
-- _Theoretical Biology and Medical Modelling_ *2*:16 (2005)
--
-- McIntyre, Richardson, Grill, 
-- _J. Neurophysiol._, *87*:995--1006, (2002)
--
-- Rattay, Aberham, Modelling axonal membranes for functional
-- electrical simulation. _IEEE Trans. Biomed. Eng._,
-- *40*(12):1201--1209, (1993)
--
-- Swarz, Reid, Bostock, Action potentials and membrane currents in
-- the human node of Ranvier, _Eur. Biophys. J_, *430*:283--292 (1995)