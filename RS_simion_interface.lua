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
  RS_simion_interface.lua

  Reaction Simulation (RS) extension to SIMION / SDS:
  Interface to SIMION / SDS

  This script is an interface between the SIMION user program
  (probably a modified version of the SDS user program) and RS.

  It
     - creates an instance of the simulation object (RS_Simulation) and controls it
     - provides some convenience methods
     - stores some callback functions which are provided by the SIMION user program which are called from the RS code

  @author Walter Wissdorf
  @version: 0.4.4
--]]

local Simulation = require "RS_Simulation"
local Particle   = require "RS_Particle"

local Interface = {} -- a local namespace for the SIMION interface
local simion_cb = {} -- a table containing all callback functions back into SIMION
local sim = {} -- a  RS simulation object


--[[
 Test method
 simple test method, if the interfacing to SIMON and the import of the local namespaces works
--]]
function Interface:test()
	print ("Kinetic SDS interface")
end


--[[
 Interface initialization:
 generates a RS simulation instance and inits the interface to it with some callback function from SIMION / (SDS)
 @param: cp_updateIonColor = a callback function from SIMION to update the ion color of the ion
 @param: cb_updateIonMass = a callback function from SIMION to update the ion mass of the ion
 @param: cb_updateIonCharge = a callback function from SIMION to update the ion charge of the ion
 @param: logIllEvents (boolean) = switch if ill events should be logged (if true ill events are logged to the log file / terminal)
--]]
function Interface:init(cb_updateIonColor, cb_updateIonMass, cb_updateIonCharge, logIllEvents)

	-- create simion callback function table
	simion_cb.updateIonColor = cb_updateIonColor
	simion_cb.updateIonMass = cb_updateIonMass
	simion_cb.updateIonCharge = cb_updateIonCharge

	--init simulation with configuration file located in the SIMION project directory
	--the conversion factor for the parsed rate constants is 1e6 because we have 1e6 microseconds / second and
	--the basic time unit in SIMION is microsecond while the rate constants are given in units of seconds
	sim = Simulation:new("./RS.conf", simion_cb, 1e6, logIllEvents)
	sim:printConfiguration() --print the configuration to the logfile
end


--[[
 Add particle
 adds a new simulation particle to RS
 @param: particleIndex = the (reference) index of the new particle in the SIMION simulation
 @param: substanceIndex = the index of the substance this new particle belongs to
 @param: x,y,z = the initial position of the particle in x,y,z direction
--]]
function Interface:addParticle(particleIndex, substanceIndex, x,y,z)
	local subst = sim.s:getDiscreteSubstanceByIndex(substanceIndex) -- the substance which is indexed by "substanceIndex" in the SIMION Part

	-- add to particle list:
	local p = Particle:new(subst, x,y,z) --create a new particle
	sim:addParticle(p, particleIndex)

	-- update the ion mass and the ion charge in SIMION (call the callback functions passed from SIMION):
	simion_cb.updateIonMass(subst:getMass())
	simion_cb.updateIonCharge(subst:getCharge())
end


--[[
 Update particle position
 updates the spatial position of a given particle
 @param: index = the (reference) index of the particle in the SIMION simulation
 @param: x,y,z = the new position of the particle
--]]
function Interface:updateParticlePosition(index, x,y,z)
	sim:P(index).x = x
	sim:P(index).y = y
	sim:P(index).z = z
end


--[[
 React a particle
 primary phase of the RS simulation:
 perform the chemical reactions of a given particle
 @param: index = the SIMION index of the particle which reacts
 @param: KE = the actual kinetic energy of the particle
 @param: dt = the time since the last reaction event / simulation timestep
--]]
function Interface:react(index, KE, dt)
	sim:react(index, KE, dt)
end

--[[
 Advance Timestep:
 Udpates the internal simulation / simulation time step statistics. If the statistics of RS are used,
 this method have to be called once for every simulation time step.

 @param: dt = the actual time step length
--]]
function Interface:advanceTimestep(dt)
	sim:advanceTimestep(dt)
end


--[[
 Print State
 prints (verbosely) the actual state of RS to the log file / console
--]]
function Interface:printState()
	self:printSimulationState()
	self:printParticles()
end

--[[
 Setups the chemical species concentration logging (the parameters of the logging mechanism are set)
 @param stepwidth = the number of time steps (calls of interface:react) after which a logging event takes place
 @param maxTimestep = the maximum number of time steps (calls of interface:react) to which logging is performed
--]]
function Interface:setupTimestepLogging(stepwidth, maxTimestep)
	self.stepwidth = stepwidth
	self.maxTimestep = maxTimestep
end

--[[
 Log time step
 logs the actual time step state in terms of the actual concentrations of the chemical species
--]]
function Interface:logTimestep()
	if sim.nsteps % self.stepwidth == 0 and sim.nsteps < self.maxTimestep then
		-- print the time step number and then the individual discrete substance concentrations
		local str = tostring(sim.nsteps)..","..tostring(sim.sumTimestep)..","
		for i, subst in sim.s:getDiscreteSubstanceIterator() do
			str = str .. tostring(sim.c[subst])..","
		end
		print(str)
	end
end


--[[
 Start runtime log
 starts the logging of the runtime
--]]
function Interface:startRuntimeLog()
	sim:startRuntimeLog()
end


--[[
 Stop runtime log
 stops the logging of the runtime, the result of the runtime measurement is written to the log file
--]]
function Interface:stopRuntimeLog()
	sim:stopRuntimeLog()
end


--[[
 Prints Simulation Statistics
 prints a overview of some statistics of the simulation (number of time steps, mean time step length,
 number of ill events) to the log file
--]]
function Interface:printSimulationStatistics()
	print("#"..sim.nsteps.." ".. sim.sumTimestep/sim.nsteps.." ill events: "..sim.illEvents)
end

--return the interface namespace
return Interface
