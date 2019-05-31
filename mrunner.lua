#!/usr/local/bin/luajit -O2

-- Run models from this file

require 'lib'   -- small utility library
require 'el_ph' -- auxilary routines and electrophysiology
require 'ode'   -- ODE solvers

--require "profiler"

_s_ = 200

--_s_ = 200 -- [global] s parameter for RKC algorithm
--local s = 95 -- minimal for h = 1e-3   |
--local s = 215  -- minimal for h = 5e-3 | Mysa0
--local s = 275  -- minimal for h = 1e-2 |
--local s = 30 --  minimal for mysa2 and h = 1e-2, 1ms->1s
--local s = 50
--local s  = 10 -- minimal for mysa2 and h = 1e-3; 1ms->3s


function safe_write(fname, str)
-- Safely write a string to a file
   local f = assert(io.open(fname,'w'))
   f:write(str)
   f:close()
end

function mymessage(str)
--wrap around stderr:write
   io.stderr:write(str..'\n')
end

function patch_table(t_to, t_from, n)
   local v = t_from[n]
   if type(v) == "string" then
      if type(t_from[n+1]) ~= "string" then
	 t_to[v] = t_from[n+1]
      else
	 patch_table(t_to[v], t_from,n+1)
      end
   end
end


--- The main function ---
function main()
   local t0, t1, h, stf, conf
   -- Default values:
   t0, h, t1 = 0, 1e-6, 10

   local tstep = 200


   stf = 'rkc_a' -- stepper function

   conf_file  = 'dc-conf'
   model_file  = 'dc-mod'
   
   local sstate_fname  -- save state file name
   local lstate_fname  -- load state file name
   local gp_fname   -- gnuplot dictionary
   
   local verbose = 0
   local loadp = false
   local savep = false
   local arg_layout 

   local patchstr

   local print_paramsp = false

   local stim = {}
    -- User defined values:
    for key,val in pairs(arg) do -- PATTERNS
       if val == '-c' then
	  conf_file = arg[key+1]     -- config file
       elseif val == '-m' then
	  model_file = arg[key+1]
       elseif val == '-h' then
	  h = tonumber(arg[key+1]) -- time step
       elseif val == '-t' then
	  t0 = tonumber(arg[key+1]) -- start time (ms)
       elseif val == '-T' then
	  t1 = tonumber(arg[key+1]) -- end time (ms)
       elseif val == '-i' then
	  stf = arg[key+1]          -- stepper function 
       elseif val == '-e' then
	  stim = assert(loadstring(arg[key+1]))()
       elseif val == "-P" then
	  patchstr = arg[key+1]
       elseif val == '-p' then
	  print_paramsp = true
       elseif val == '-s' then
	  _s_ = tonumber(arg[key+1])
       elseif val == '-w' then
	  tstep = tonumber(arg[key+1])
       elseif val == '-L' then
	  arg_layout = arg[key+1]
       elseif val == '-v' then 
	  verbose = verbose + 1
       elseif val == '-l' then
	  loadp = true
	  if arg[key+1] then 
	     lstate_fname = arg[key+1]
	  end
       elseif val == '-S' then
	  savep = true
	  if arg[key+1] then
	     sstate_fname  = arg[key+1]
	  end
       end

    end
    
    conf = require(conf_file)
    if patchstr then
       patch = assert(loadstring(patchstr))()
       patch_table(conf,patch,1)
    end
    model = require(model_file)


    for name, segment in pairs(conf.segments) do
       if segment.pars.stimulated then
	  segment.pars.iappl = stim
       end
    end

    if arg_layout then 
       layout_rule = conf[arg_layout]
    end
    
    nn = model.Nerve:new(conf.segments,layout_rule)
    
    if loadp then nn:luaload_state(lstate_fname) end
     
 
    if not sstate_fname then
       sstate_fname = nn.layout_rule.tag..'-state.lua'
    end

    --- Write the Gnuplot dictionary ---
    gp_fname = nn.layout_rule.tag..'-dict.gp'
    safe_write(gp_fname, nn:gp_print())

    if verbose > 0 then
       io.stderr:write(nn:pp_state(t0)..'\n')
    end

    if print_paramsp then
       io.stderr:write(nn:print_parameters())
    end
    
    --- Run the simulation ---
    local t = t0
    if profiler then profiler.start() end
    while t <= t1 do
       local phase
       local t2 = (t+tstep > t1) and t1 or t+tstep
       t,phase = nn:run(t, h, t2, stf)
       --- Save state ---
       if savep then 
	  safe_write(sstate_fname, nn:luaprint_state(phase, tend))
       end
    end
    io.stderr:write('\n')
    if profiler then profiler.stop() end

 end

main()
