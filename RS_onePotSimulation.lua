--[[
  Monte Carlo Reaction Simulation (RS) extension to SIMION charged particle simulator

  This program implements a monte carlo model for the simulation of chemical
  reaction kinetics involving charged particles (typically molecular ions) in the
  SIMION charged particle simulator.
  Additionally this program allows the simulation of reaction kinetics in an
  ideally stirred reactor in a standalone mode without SIMION.

  Copyright (C) 2012 - Physical and Theoretical Chemistry /
  Institute of Pure and Applied Mass Spectrometry
  of the University of Wuppertal, Germany


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
  RS_onePotSimulation

  This script implements a simple "one pot reaction simulation", a monte
  carlo simulation of chemical reaction dynamics in an ideal stirred reactor.

  This script is primarily used to test the reaction simulation extension to
  SIMION / SDS without an actual ion trajectory simulation, but if the ideal
  stirred reactor assumption is valid in a simulation problem, the script
  could also be used for productive simulations.

  This script is a simple template, to show how the Reaction Simulation
  extension can be used directly, without SIMION.
  More complex usages could easily derived from it.

  Invocation / command-line arguments:
   lua RS_onePotSimulation.lua nSimSteps maximumDt nParticles outputFileName

   nSimSteps = number of simulation time steps to simulate
   maximumDt = the maximum length of a simulation time step (the time step
               (lengths are equally random distributed below this value)
   nParticles = the number of particles to simulate
   outputFileName = name of the simulation results file / simulation log file

  @author Walter Wissdorf
  @version: 0.4.4
--]]

local Simulation = require "RS_Simulation"
local Particle   = require "RS_Particle"

-- The pepperfish lua profiler can be used to profile the Reaction Simulation
-- code in depth:
-- http://lua-users.org/wiki/PepperfishProfiler
local useProfiler = false

if useProfiler then
	require "profiler"
end


-- check if all commandline arguments are given,
-- if not, give usage hint and terminate:
if #arg ~= 4 then
	print("Wrong argument number. ")
	print("Usage information: ")
	print("lua RS_onePotSimulation.lua nSimSteps maxdt nParticles outputFileName")
	os.exit()
end

-- parse the given commandline arguments
local nSimSteps = tonumber(arg[1])
local maxdt = tonumber(arg[2])
local nParticles = tonumber(arg[3])
local outputFileName = arg[4]
local outputFile = io.open(outputFileName, "w")


-- start the profiler
if useProfiler then
	profiler = newProfiler('call', 100)
	profiler:start()
end
local tStart = os.time() -- start runtime measurement

-- create and init a new RS.Simulation object - an actual simulation instance:
-- (the time basic unit is arbitrary in one pot simulations, therefore the
-- reaction rate conversion factor is 1).
-- init simulation with reaction system
local sim = Simulation:new("./RS.conf", nil, 1, true)
sim:printConfiguration()

-- add a bunch of particles, all initialized with the first existing substance:
local function addParticle(i, substanceIndex)
	local subst = sim.s:getDiscreteSubstanceByIndex(substanceIndex)
	local p = Particle:new(subst, sim.random(), sim.random(), sim.random())
	sim:addParticle(p, i)
end
for i = 1, nParticles do
	addParticle(i, 1)
end

sim:printSimulationState()


-- The main simulation loop:
-- perform the reaction method of the simulation object for every particle
-- in every time step
for ts=1, nSimSteps do
	-- time step lengths are equally random distributed with an given maximum value
	local dt = sim.random()*maxdt

	sim:advanceTimestep(dt)
	for p = 1, nParticles do
		sim:react(p, 0, dt)
	end

	-- print the actual simulation state and
	-- log the concentrations of the simulated chemical species to the output file
	sim:printSimulationState()
	outputFile:write(sim:getConcentrationString())
end
-- -------------------------

-- write some statistics to the output file
outputFile:write(
	" ill events: "..sim.illEvents ..
	" mean dt: "..tostring(sim.sumTimestep / sim.nsteps)
)
print("ill events:"..sim.illEvents)


-- stop runtime measurement and print used wall time:
local tStop = os.time()
print("time used: "..tStop-tStart)

-- if we use the profiler, get the results and write them to the result file:
if useProfiler then
	profiler:stop()
	profiler:report(outputFile)
end

outputFile:close()
