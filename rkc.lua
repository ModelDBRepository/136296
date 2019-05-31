-- (almost) direct translation from rkc.f [http://www.netlib.org/ode/rkc.f]
-- Implementation of Runge-Kutta-Chebyshev explicit ODE solver
--
-- Plan:
-- 1. basic stage: fixed timestep and s (order of Chebyshev polynomials) [DONE]
-- 2. error estimation & control     [in progress]
-- 3. automatic time step adjustment [almost-DONE]
-- 4. automatic s
-- 5. spectral radius
--

--require "profiler"

require 'lib'

module(..., package.seeall);

function rkc_v1(t,h,s,obj) -- first version
    local c,mu,mut,nu,gammat = calc_coefs(s)
    local phase = obj.phase       
    local dv = obj.derivatives

    local W0 = phase 
    local F0 = obj:derivs(t,phase)
    local W1 = phase + mut[1]*h*F0
    local Wm1,Wm2 = W1,W0
    local W = W0+0
    
    local thr = 1e3
    
    for j = 2,s do
        if math.abs(absmax_phase(W)) > thr then
            io.stderr:write('j: ',j,', mu: ',
                mu[j],', nu: ', nu[j], ', mut: ',
                mut[j],', gammat: ', gammat[j],'\n')
            io.stderr:write(obj:pr_state(t,W))
            error(string.format("the solution has ran away: %f",
				absmax_phase(W) ))
        end
        W = (1 - mu[j] - nu[j])*W0 + mu[j]*Wm1  + nu[j]*Wm2 
        W = W + mut[j]*h*(obj:derivs(t+c[j-1]*h,Wm1) + gammat[j]*F0)
        Wm2 = Wm1 -- wrap around
        Wm1 = W
    end
    --io.stderr:write('\n')
    return W
end


function calc_coefs(s) 
-- calculate coefficients for number of stages s. 
-- Return vectors. Used in rkc_v1
    local coefs = {}
    local epsilon,w0,w1,tmp1,tmp2,arg
    local b = {}
    local mu = {}
    local nu = {}
    local mut = {} -- mu with tilda
    local c = {}   -- step size adjustments
    local gammat = {} -- gamma with tilda
    local T,dT,d2T = {},{},{} -- Chebyshev polynomial and its derivatives
    
    epsilon =  2/13
    w0 = 1 + epsilon/(s*s)    
    
    tmp1 = w0*w0 - 1        -- temporary var
    tmp2 = math.sqrt(tmp1)  -- temporary var
    arg = s*math.log(w0 + tmp2)
    
    -- w1 = T's(w0)/T''s(w0)
    w1 = math.sinh(arg)*tmp1 / (math.cosh(arg)*s*tmp2 - w0 *math.sinh(arg))

    --?--
    b[0] = .25/(w0*w0)
    b[1] = b[0]

    mut[1] = w1*b[1]  
    
    T[0] = 1 -- todo: check!
    T[1] = w0
    dT[0] = 0
    dT[1] = 1
    d2T[0] = 0
    d2T[1] = 0
    

    c[0]= 0           --
    c[1] = mut[1]     -- as in rkc.f

    for j = 2,s do
        T[j] =   2*w0*T[j-1] - T[j-2]
        dT[j] =  2*w0*dT[j-1] - dT[j-2] + 2*T[j-1]
        d2T[j] = 2*w0*d2T[j-1] - d2T[j-2] + 4*dT[j-1]

        b[j] = d2T[j] / (dT[j] * dT[j])
        mu[j] = 2 * w0*b[j] / b[j-1]
        nu[j] = -b[j] / b[j-2]
        mut[j] = mu[j] * w1 / w0
        gammat[j] = -1 + T[j-1]*b[j-1]

        c[j] =  mu[j]*c[j-1] + nu[j]*c[j-2] + mut[j]*(1 + gammat[j]) -- as in the rkc.f code
        --c[j] = w1*d2T[j]/dT[j] -- as in the paper (B.P. Sommeijer et. al, 1997), not in the rkc.f code
    end
    --c[1] = c[2]/dT[2]
    --c[s] = 1
    
    return c,mu,mut,nu,gammat

end


function rkc_v3(W0, t,h,s,obj) -- third version (speed optimisation)
    local eps,w0,w1,tmp1,arg
    eps = 2/13
    w0  = 1+ eps/(s*s)
    tmp1 = w0*w0 - 1
    arg = s*math.log(w0 + math.sqrt(tmp1))
    w1 = math.sinh(arg)*tmp1 / (math.cosh(arg)*s*math.sqrt(tmp1) - 
			  w0 *math.sinh(arg))

    local cfs_2 = { -- coefs_{j-2} 
        b = .25/(w0*w0),
        T = 1, dT = 0, d2T = 0, c = 0,
    }

    local cfs_1 = { -- coefs_{j-1}
        b = cfs_2.b, T = w0, dT = 1, d2T = 0,
        mut = w1*cfs_2.b,
        c = w1*cfs_2.b, -- as in rkc.f
    }


    local dv

    local F0 = obj:derivs(t, W0) + 0 -- derivs with the original phase
    local Wm1 = W0 + F0*(cfs_1.mut*h)
    local Wm2 = W0 + 0
    local W = W0 + 0 --setmetatable({}, lib.vecmt)

    
    local cfs = {}
    for j = 2,s do
        cfs.T   = 2*w0*cfs_1.T - cfs_2.T
        cfs.dT  = 2*w0*cfs_1.dT - cfs_2.dT + 2*cfs_1.T
        cfs.d2T = 2*w0*cfs_1.d2T - cfs_2.d2T + 4*cfs_1.dT

        
        cfs.b = cfs.d2T / (cfs.dT * cfs.dT)
        cfs.mu = 2*w0*cfs.b/cfs_1.b
        cfs.nu = -cfs.b / cfs_2.b
        cfs.mut = cfs.mu * w1 / w0
        cfs.gat = -1 + cfs_1.T * cfs_1.b

        --as in the rkc.f
        --cfs.c = cfs.mu*cfs_1.c + cfs.nu*cfs_2.c + cfs.mut*(1 + cfs.gat) 

        -- as in the paper (B.P. Sommeijer et. al, 1997), 
	-- not as in the rkc.f code
        cfs.c = w1*cfs.d2T/cfs.dT 

        dv = obj:derivs(t + cfs_1.c*h, Wm1) -- calc derivs with another phase
	for n,val in ipairs(W0) do
	   W[n] = ((1 - cfs.mu - cfs.nu)*val + 
		cfs.mu*Wm1[n] + cfs.nu*Wm2[n] +
		cfs.mut*h*(dv[n] + cfs.gat*F0[n]))
	   ----------
	   Wm2[n] = Wm1[n]
	   Wm1[n] = W[n]
	end

        --wrap around --TODO: redefine as a tail recursion :)
	--(cp_vector target source) [desctructive for target]
        --lib.cp_vector(Wm2,Wm1)
        --lib.cp_vector(Wm1,W)
        lib.cp_simp_table(cfs_2,cfs_1)
        lib.cp_simp_table(cfs_1,cfs)
	
	--io.stderr:write(string.format("::: %3.20e :::\n\n", 
	--			      lib.sum_sq_vect(W)))

     end
     return W, t+h, h
end


function Estnm1(U0, U1, dv0, dv1, h)
   local E = setmetatable({}, lib.vecmt)
   for i,val in ipairs(U0) do
      E[i] = ((U0[i] - U1[i])*2 + h*(dv0[i] + dv1[i]))*2/5
   end
   return E
end

function Est2(a,b,h)
   return (a*12 + 6*h*b)/15
end

local hprev = 0.0

function rkc_a(U0, t, h, s, obj)
   
   -- rkc with adjustable timestep [gives ca 3x speedup]
   --profiler.start()

   local recipr_atol = 1/1e-3

   local nmax = 1e3
   local j,n = 0,0

   local dv0 = obj:derivs(t, U0) + 0
   local dv1
   local U1, E,  err
   local rec = 0.85 -- h recalculating coefficient
   
   h = h/rec
   repeat
      U1 = rkc_v3(U0, t, h, s, obj)
      dv1 = obj:derivs(t+h, U1)
      E = Estnm1(U0, U1, dv0, dv1, h) * recipr_atol
      err,n = lib.sum_sq_vect(E)
      err = math.sqrt(err/n)
      h = h*rec
      j = j+1
   until (err <= 1 or j > nmax) 
   
   h = h/rec
   local tnext = t + h
   --profiler.stop()
   --io.stderr:write(string.format('-------------\n'))
return U1, tnext, h
end


rkc = rkc_v3 -- pointer to the active verstion
