--[[
  Monte Carlo Reaction Simulation (RS) extension to SIMION charged particle simulator

  This program implements a monte carlo model for the simulation of chemical
  reaction kinetics involving charged particles (typically molecular ions) in the
  SIMION charged particle simulator.
  Additionally this program allows the simulation of reaction kinetics in an
  ideally stirred reactor in a standalone mode without SIMION.

  Copyright (C) 2012 - Department of Physical and Theoretical
  Chemistry of the university of Wuppertal, Germany


  This file is part of Monte Carlo Reaction Simulation (RS)

  RS is free software: You may redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.


  The distribution of this program shipping with SIMION as an example is licensed exclusively
  under the SIMION 8.1 license.
  For derivative works based on the version distributed with SIMION,
  the GPL license applies as stated above.

  ------------
  RS_Simulation.lua

  Reaction Simulation extension to SIMION / SDS: Main simulation class

  This is the main class of the chemical reaction dynamic simulation extension to the SDS implementation of SDS.
  The management of the simulated particles and the core simulation code (calculation of the reaction probability etc.)
  are located here

  The simulation code is encapsulated in its own class, all methods of the RS simulation object are exposed via
  functions in the "Simulation" namespace

  @author Walter Wissdorf
  @version: 0.4.4

--]]

local Parser       = require 'RS_ConfigFileParser'
local ParticleList = require 'RS_ParticleList'
local Particle     = require 'RS_Particle'

--TODO: Calculate the ratio of ill events (count all reaction events and calculate ratio = ill events / all events)


local Simulation = {} -- namespace of the Simulation class

	--[[
		Constructor:
		(constructor according to "Programming in Lua" online literature source)
		constructs a new instance of the simulation class
		@param: conf_file = the filename of the configuration file which contains the simulation configuration.
				This file is imported and parsed by a config file parser

		@param: simon_callbacks = the callback function table
				This is a table with necessary functions in the actual SIMION user program (in the case of RS probably
				a modified version of the SDS - statistical diffusion simulation - code).
				The necessary callback functions are:
					- updateIonMass(newMass) updates the mass of the actual ion with the newMass
					- updateIonCharge(newCharge) updates the charge of the actual ion with newCharge
					- updateIonColor(newColorIndex) updates the color index of the actual ion with newColorIndex
		@param: rateConstConversionFactor = a factor by which the given rate constants are divded to allow conversion of the basic time unit
		@param: logIllEvents (boolean) = if true every single ill reaction event will be logged to the log / terminal

		@return: a new instance of the RS_Simulation class

	--]]
	function Simulation:new(conf_file, simion_callbacks, rateConstantConversionFactor, logIllEvents)

		local sim = {} --create a new instance of the simulation class

		--set class metatable and sets the base class (SubstanceTable) as index
		--(setups the inheritance of the new instance to the base class)
		setmetatable(sim, self)
		self.__index = self --set index
		-- --------------------------------------------------------------------

		--some constant member fields: ----------------------------------------
		sim.randomWalkScale = 0.01
		sim.searchRadius = 0.01
		sim.debug = false

		-- simulation initialization: -------------------------------------------------------------------------------------
		-- (setup of internal state)

		sim.logIllEvents = logIllEvents

		sim.simion_cb = simion_callbacks
		
		-- Pseudorandom number generator (PRNG), generating uniformly
		-- random number in range [1, 0).
		-- We prefer SIMION's random number generator since ANSI rand()
		-- is not ensured to be a good PRNG.
		-- Note: simion.seed(n) can seed SIMION's random number generator.
		sim.random = (simion or {}).rand or math.random

		-- open and parse config file:
		Parser:openFile(conf_file)
		sim.s, sim.r = Parser:parse(rateConstantConversionFactor) -- sim.s => the substance table; sim.r => the reactions table

		sim.ri= {} -- the substance specific independent reactions dict (independent == only one discrete educt)
		sim.rd= {} -- the substance specific dependent reactions dict (dependent == at least two discrete educts)
		sim.ri_statProbs = {} --the substance specific independent static reaction probabilities (saved for fast access in the main simulation method performance)
		-- sim.ri, sim.rd and sim.ri_statProbs have the following structure:
		--   -they are dictionaries, with the discrete substances as keys
		--   -for every discrete substance it contains a table (list) which contains:
		--   -- the independent reactions (sim.ri)
		--   -- the dependent reactions (sim.rd)
		--   -- the static reaction probabilities of the independent reactions (sim.ri_statProbs)

		-- initialize internal data structures and variables:
		sim.p = ParticleList:new() -- the particle list (a double linked list of all particles)
		sim.ionMap = {} -- a map from the SIMION ion indices to the particles in the kinetic simulation
		sim.c = {} -- the concentrations dictionary, a dictionary which contains the absolute numbers of the discrete particles (per substance)
		sim.nsteps = 0 -- the actual time step of the simulation
		sim.illEvents = 0 -- the number of "ill" reaction events (actually an ill reaction event is if the reaction probability is >1)
		sim.sumTimestep = 0 --the sum of all time steps
		sim.startTime = 0 -- the time of simulation start



		-- init substance specific reaction tables ----------------------------------------------------------------
		-- -- init dictionaries
		for i, subst in sim.s:getSubstanceIterator() do
			sim.ri[subst] = {}
			sim.rd[subst] = {}
			sim.ri_statProbs[subst] = {}
		end

		-- -- iterate through all reactions
		for _, r in ipairs(sim.r) do
			if r:isIndependent() == true then
				-- the reaction is independent from other discrete educts: therefore insert it to the independent reactions table
				local ed = next(r:getDiscreteEducts(), nil)    -- the reac. is independent=> only one element in the discrete educts table
				table.insert(sim.ri[ed], r)
				table.insert(sim.ri_statProbs[ed], r:getStaticProb())
			else
				-- the reaction is dependent on more than one discrete educt:
				--   insert it to the dependent reaction table for every single discrete educt:
				for educt, factor in pairs(r:getDiscreteEducts()) do
					table.insert(sim.rd[educt], r)
				end
			end
		end
		-- ---------------------------------------------------------------------------------------------------------


		-- test / debug print of independent simulations ----------------
		if sim.debug then
			for subst, ri in pairs(sim.ri) do
				if #ri > 0 then
					print(subst:getName())
					for j, reac in ipairs(ri) do
						reac:printReaction()
					end
				end
			end
		end
		-- --------------------------------------------------------------


		-- init concentrations dict ---- --------------------------------
		for _, subst in sim.s:getDiscreteSubstanceIterator() do
			sim.c[subst] = 0
		end
		-- --------------------------------------------------------------

		return sim --everything done, return the new instance of the simulation class
	end


	--[[
	 Print Reaction Tables:
	 Prints the actual reaction tables of the individual substances
	 in a formatted manner to the console / logfile
	--]]
	function Simulation:printReactionTables()
		print("simulation reaction tables: --------------------------------")
		for _, subst in self.s:getSubstanceIterator() do --iterate trough all substances
			if subst:getType() ~= "discrete" then --if the substance is not discrete => move on
				break
			end

			-- print out the substance name and the independent and dependent reaction tables
			print("\n\n  >"..subst:getName())

			print("  >>> independent: ")
			for _, r in ipairs(self.ri[subst]) do
				r:printReaction()
			end

			print("  >>> dependent: ")
			for _, r in ipairs(self.rd[subst]) do
				r:printReaction()
			end

		end
		print("---------------------------------------------------------\n")
	end


	--[[
	 Print simulation state:
	 Prints the state of the simulation to the console / logfile
	 The state of the simulation is:
	     - the actual time step of the simulation
	     - the actual concentrations of the simulated discrete species
	     - if in debug mode: the positions and the actual substances of all simulated particles
	--]]
	function Simulation:printSimulationState()
		print("simulation state: ------------------------------------")
		print("sum timestep: ".. self.sumTimestep)


		for i, subst in self.s:getDiscreteSubstanceIterator() do
			print("  >"..subst:getName().."  >".. self.c[subst])
		end
		print("------------------------------------------------------\n")


		if self.debug then
			for part in self.p:getListIterator() do
				print("x:"..part.x.." y:"..part.y.." substance:"..part.s:getName())
			end
		end
	end


	--[[
	 Get concentration string:
	 Returns the concentrations (absolute numbers in this case) of the individual simulated species
	 in a formatted string, which can be written to the simulation log etc.

	 @returns: a string with the actual numbers of the simulated species
	--]]
	---
	function Simulation:getConcentrationString()
		local str = tostring(self.sumTimestep)..";"
		for i, subst in self.s:getDiscreteSubstanceIterator() do
			str = str .. self.c[subst]..";"
		end
		str = str .."\n"
		return str
	end


	--[[
	 Print Configuration:
	 Prints the configuration (Substances, Reactions) of the simulation to the console / log
	--]]
	function Simulation:printConfiguration()
		print("parsed substances: -----------------------------------")
		local bufStr
		for k, subst in self.s:getSubstanceIterator() do
			bufStr = subst:getName().." "..subst:getType()
			if subst:getType() ~= "discrete" then
				bufStr = bufStr .. " "..subst:getStaticConcentration()
			end
			print(bufStr)
		end
		print("------------------------------------------------------\n")


		print("parsed reactions: ------------------------------------")
		for i, r in ipairs(self.r) do
			r:printReaction()
		end
		print("------------------------------------------------------\n")
	end


	--[[
	 Print Particles:
	 Prints the spatial positions and the actual substance names of all particles in the simulation
	 to the log / console
	--]]
	function Simulation:printParticles()
		print('particle positions: ')
		for index, p in pairs(self.ionMap) do
			print(index, p.x,p.y,p.z, p.s:getName())
		end
	end


	--[[
	 Random Walk:
	 Performs a random walk of the positions of all particles of the simulation.
	 The mean random walk distance is scaled with the randomWalkScale member variable.
	 (The method is typically used to update the particle positions in a one pot simulation (without SIMION / SDS) )
	 --]]
	function Simulation:randomWalk()
		local dx=0
		local dy=0

		--iterate trough all particles
		for particle in self.p:getListIterator() do
			-- random walk / position update (wrap around if outside domain )
			dx = (self.random()-0.5)* self.randomWalkScale
			dy = (self.random()-0.5)* self.randomWalkScale
		
			particle.x = particle.x + dx
			particle.y = particle.y + dy

			-- wrap around ...
			if particle.x > 1 then
				particle.x = particle.x - 1
			end

			if particle.x < 0 then
				particle.x = 1 + particle.x
			end

			if particle.y > 1 then
				particle.y = particle.y - 1
			end

			if particle.y < 0 then
				particle.y = 1 + particle.y
			end
			-- ----------------
		end
	end


	--[[
	 Add particle:
	 adds a new discrete particle to the chemical reaction simulation

	 @param: particle = the particle to add to the simulation particle table
	                    (an instance of the RS_Particle class)

	 @param: index = the external particle index of the particle which can be used to access the new particle
	                 (typically this would be the index of the particle in SIMION)
	--]]
	function Simulation:addParticle(particle, index)
		self.p:insert(particle)
		self.c[particle.s] = self.c[particle.s]+1
		self:setP(index, particle) --update the external index / particle map
	end


	--[[
	 Destroy Particle:
	 destroys / removes a particle from the reaction simulation

	 @param: particle = the particle to destroy / remove from the particle list
	                    (an instance of the RS_Particle class)
	--]]
	function Simulation:destroyParticle(particle)
		local subst = particle.s --the substance which the particle is belonging to
		self.p:remove(particle)--remove the particle from the particle list
		self.c[subst] = self.c[subst] -1 --update concentration table
	end


	--[[
	 Distance:
	 gets the distance between two simulated particles

	 @param: p1,p2 = two particles (instances of RS_Particle) for which the distance between them is calculated

	 @return: spatial distance between p1 and p2
	--]]
	--[TODO: implement for three spatial dimensions]
	--local function distance(p1, p2)
	--	return math.sqrt( (p1.x-p2.x)^2 + (p1.y-p2.y)^2 )
	--end


	--[[
	 Add Products:
	 adds the products of a reaction (the particles resulting from a chemical reaction) to the
	 simulation (to the simulation particle table)

	 @param: reaction = a reaction (instance of RS_Reaction) for which the products are added
	 @param: xPos,yPos,zPos = the spatial position of the reaction event in x,y,z direction
	--]]
	function Simulation:addProducts(reaction, xPos,yPos,zPos)
		-- add products
		for product, factor in pairs(reaction:getProducts()) do -- at the moment only discrete products are allowed
			if product:getType() == "discrete" then --if the product is discrete (actually enforced by the config file parser)
				for i=1, factor do
					self:addParticle(Particle:new(product, xPos,yPos,zPos))
				end
			end
		end
	end


	--[[
	 React:
	 the core function of the simulation: performs the reactions of a given particle

	 Notice:
	 Currently only independent reactions (reactions with only one discrete particle on the educt) side are
	 considered, the dependent reactions are parsed but ignored silently in the simulation loop actually

	 (This method was profiled and optimized for execution speed, therefore the implementation is not as straightforward
	 as in the parts which are not that critical for the simulation speed)

	 @param: index = the (external) index of the ion / SIMION particle which reacts
	 @param: dt = the time since the last reaction event / simulation timestep
	 @param: KE = the kinetic energy of the particle
	 @param: pressure = the pressure (in units of the norm pressure for the "k" values)
	--]]

	-- TODO: pressure dependence of the reactions:
	-- an elementary reaction scales linearly with the pressure?
	-- an elementary reaction scales linearly with the collision frequency (!)
	-- Check the dependence of the collision frequency in dependence to the pressure

	function Simulation:react(index, KE, dt)

		-- Only independent reactions are currently considered!

		-- local dead = false -- flag if the particle is dead already, only needed if we consider ion-ion reactions

		-- The basic scheme of the react method is:
		--   run independent reactions:
		--     for all independent reactions of the actual particle do:
		--       check if a random value is below the static-reaction probability * the time step length (dt)
		--       if yes
		--          react => create product particles, destroy educt particles
		--
		local p = self.ionMap[index] --get the particle associated with the ion with index "index" didn't use Simulation:P(index) because of performance

		local iReactions = self.ri[p.s] --independent reactions of the actual particle
		local iReactionsStaticProbs = self.ri_statProbs[p.s]
		for i=1, #iReactions do  -- didn't use "for i,iReac in ipairs(self.ri[p.s]) do" because of performance

			local iReac = iReactions[i]  -- iReac is the i-th independent reaction of the substance (p.s) of particle p

			-- now check if we react or not:
			--  get a random value "r" if r is below the reaction probability => react
			local r = self.random()

			-- If we consider reactions with an activation energy, we have to modify the reaction probability
			-- by an reasonable model of the effect of a significantly nonthermal kinetic energy.
			-- The code is currently deactivated, because we actually don't need it and we can save the performance for

			-- if iReac:getActivationEnergy() and KE > iReac:getActivationEnergy() then
				--if the KE is higher than the activation energy (and we have an activation energy set)
				--we definitively react
			--	prob = 1
			--else
				--local prob = iReac:getStaticProb() * dt
			--end
			--instead we silently ignore the kinetic energy of the simulated particle:
			local prob = iReactionsStaticProbs[i] * dt

			if r < (prob) then

				if prob >= 1 then --if the probability is >1 the time step was too long, and the reaction event was ill
					self.illEvents = self.illEvents +1
					if self.logIllEvents== true then
						-- log the ill events to the log / terminal
						print("ill reac. event, p:"..prob.." reac:"..iReac:getStringRepresentation())
					end
				end

				-- despite the possible illnes of the reaction event: now REACT! => remove actual particle
				-- (independent reaction => therefore only on discrete educt and this must be the actual particle)
				self:destroyParticle(p)

				-- add /append product particles (actually, there is only one possible discrete product, because
				-- the number of particles cannot increase (no ion generation in SIMION possible at present) and particle destruction
				-- reactions are not yet implemented, an independent reaction has only one discrete educt.
				-- Actually the config file parser enforces only one discrete products)

				local ps = iReac:getDiscreteProductsTable()[1] --the product substance
				local pp = Particle:new(ps, p.x,p.y,p.z)--the product particle
				self:addParticle(pp, index) --addParticle updates the external index / particle mapping also

				if self.simion_cb then --if the simion callback functions are defined (if we are doing a simulation connected to simion):
					self.simion_cb.updateIonMass(ps:getMass())
					self.simion_cb.updateIonCharge(ps:getCharge())
					self.simion_cb.updateIonColor(self.s:getSubstanceIndex(ps))
				end

				--end the independent reactions loop => the actual educt particle does not longer exist
				break
				-- [TODO (when dependent reactions are implemented):]
				-- a "dead" flag has to be set => no dependent reactions could happen now, the particle is dead
			end
		end
		-- [TODO implement dependent reactions]
		-- --------------------------------------------------------
	end


	--[[
	 Advance Timestep:
	 Udpates the internal simulation / simulation time step statistics. If the internal statistics are used,
	 this method have to be called once for every simulation time step.

	 @param: dt = the actual time step length
	--]]
	function Simulation:advanceTimestep(dt)
		self.nsteps = self.nsteps +1
		self.sumTimestep = self.sumTimestep + dt
	end


	--[[
	 P: (get particle)
	 returns the particle in the kinetic simulation associated with a ion (with the index "index") in SIMION

	 @param: index = the (external) index of the particle to return

	 @return: the particle (instance of class RS_Particle) with the given external index
	--]]
	function Simulation:P(index)
		return self.ionMap[index]
	end


	--[[
	 Set P:
	 sets the external ion index "index" of the particle "particle"

	 @param: index = the new external index
	 @param: particle = the particle (instance of RS_Particle class) to set the external index for
	--]]
	function Simulation:setP(index, particle)
		self.ionMap[index] = particle
	end


	--[[
	 Remove P (Remove Particle):
	 removes a particle (with the external index "index") from the internal Particle index map

	 @param: index = the (external) index of the particle to remove
	--]]
	function Simulation:removeP(index)
		self.ionMap[index] = nil
	end


	--[[
	 Start Runtime Log:
	 starts the clock for the simulation run time measurement
	--]]
	function Simulation:startRuntimeLog()
		self.startTime = os.time()
	end


	--[[
	 Stop Runtime Log:
	 stops the clock for the simulation run time measurement and prints the result of the
	 runtime measurement to the logfile / the console
	--]]
	function Simulation:stopRuntimeLog()
		local runtime = os.time() - self.startTime
		print('Runtime: '..runtime.. ' startTime:'.. self.startTime)
	end

--- end of class definition "Simulation" --- --- --- --- --- --- --- --- --- --- --- --- --- ---

return Simulation
