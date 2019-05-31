-- General purpose functions
module(..., package.seeall);


-- basic inheritance 
function inheritFrom(baseClass)
    local new_class = {}
    --local class_mt = {__index = new_class}
    new_class.__index = new_class

    function new_class:new(phase,params,virtuals,tag)
        local new_inst = {}
        setmetatable(new_inst, new_class)
        --new_inst.__index = self
        new_inst.v0 = deepcopy(phase)
        new_inst.p = deepcopy(params)
	new_inst.virt = deepcopy(virtuals)
	new_inst.memoized = {}
	self.base_tag = tag
        return new_inst
    end

    -- key to inheritance
    if baseClass then
        setmetatable (new_class, {__index=baseClass})
    end

    return new_class
 end
-----------


-- Utility functions --

function logistic(a,k,s,f)
   return function(x)
	     return f + a/(1+math.exp(-k*(x-s)))
	  end
end

function hill(a,k,s,f)
   return function(x)
	     return f + a/(1+ (x/s)^-k)
	  end
end


function sum_pow_phase(phase, pow)
   local res = 0
   local n = 0
   for i, node in ipairs(phase) do
      for j,val in pairs(node) do
	 res = res + val^pow
	 n = n+1
      end
   end
   return res,n
end

function sum_sq_vect(vect)
   local res,n = 0,0
   for i, val in ipairs(vect) do
      res = res + val*val
      n = n+1
   end
   return res, n
end

function cp_phase(target,source)
    for n,comp in ipairs(source) do
        for name,val in pairs(comp) do
            target[n][name] = val
        end
    end
end

function absmax_vect(vect)
    local x = -1
    for i,val in ipairs(vect) do
       if math.abs(val) > x then x = math.abs(val) end
    end
    return x
 end

function cp_simp_table(target,source)
    for key,value in pairs(source) do
        target[key] = value
    end
end


function cp_vector (target,source)
   if not target then 
      target = setmetatable({}, vecmt) 
   end
   for j,val in ipairs(source) do
   --for j = 1, # source do
      target[j] = source[j]
   end
   return target
end

------------------------------------------------------------------
-- metatable for state addition and multiplication
--

function checktypes (a, b)      -- check if both a,b are tables and if not, put the table first
        local btf = true
        if not (type (a) == "table" and type(b) == "table") then
            btf = false -- only one table
            if type(a) ~= "table" then 
                a,b = b,a  -- put the table first
            end
        end
        return btf,a,b
    end

-- Metatable to add and multiply states TODO: rewrite to make general
-- i.e. 1) pointwise multiplication or summation of the tables of equal sizes and
--      2) pointwise multiplication or summation of table and a constant

stmt = {
    __add = function(a, b) -- addition
            local btf
            btf,a,b = checktypes(a,b)

            local res = {}
            for k,node in ipairs(a) do
                res[k] = {}
                for index,value in pairs(node) do
                    res[k][index] = value + (btf and b[k][index] or b)
                end
            end
            return setmetatable(res,stmt)
    end,
    __sub = function(a, b) -- substraction
            --local btf
            --btf,a,b = checktypes(a,b)

            local res = {}
            for k,node in ipairs(a) do
                res[k] = {}
                for index,value in pairs(node) do
		   res[k][index] = value - b[k][index] --(btf and b[k][index] or b)
                end
            end
            return setmetatable(res,stmt)
    end,
    __mul = function(a, b) -- multiplication
        local btf
        btf,a,b = checktypes(a,b)
        local res = {}
        for k,node in ipairs(a) do
            res[k] = {}
            for index,value in pairs(node) do
                res[k][index] = value * (btf and b[k][index] or b)
            end
        end
        return setmetatable(res,stmt)
    end,
}

vecmt = {
   __add = function(a, b) -- addition
	      local btf
	      btf,a,b = checktypes(a,b)
	      
	      local res = {}
	      for k,value in ipairs(a) do
		 res[k] = value + (btf and b[k] or b)
	      end
	      return setmetatable(res,vecmt)
	   end,
   __sub = function(a, b) -- substraction
	      local btf
	      btf,a,b = checktypes(a,b)
	      
	      local res = {}
	      for k,value in ipairs(a) do
		 res[k] = value - (btf and b[k] or b)
	      end
	      return setmetatable(res,vecmt)
	   end,
   __mul = function(a, b) -- multiplication
	      local btf
	      btf,a,b = checktypes(a,b)
	      
	      local res = {}
	      for k,value in ipairs(a) do
		 res[k] = value * (btf and b[k] or b)
	      end
	      return setmetatable(res,vecmt)
	   end,
}


function deepcopy(object) -- returns a deep copy of the table
    local lookup_table = {}
    local function _copy(object)
        if type(object) ~= "table" then
            return object
        elseif lookup_table[object] then
            return lookup_table[object]
        end
        local new_table = {}
        lookup_table[object] = new_table
        for index, value in pairs(object) do
            new_table[_copy(index)] = _copy(value)
        end
        return setmetatable(new_table, getmetatable(object))
    end
    return _copy(object)
end
