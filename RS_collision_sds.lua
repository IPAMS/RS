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
  Modified version of SDS: Chemical reactions simulation in SIMION / SDS

  This version of the SDS code implements chemical kinetic simulation
  (monte carlo simulation of chemical kinetics of charged particles)
  in SIMION.
  
  Version of the RS extension:
  @version: 0.4.4
 --]]
simion.workbench_program()

--the SDS model can also be imported from its installation position in the SIMION folder
local SDS = simion.import('collision_sds.lua', '8.1.1.15.20120807 noinstall')
assert(rawget(SDS, '_VERSION'), 'please upgrade collision_sds.lua')

-- import kinetic simulation interface methods
local RS = simion.import 'RS_simion_interface.lua'

-- Reaction Simulation specific adjustables:

-- Logging frequency.  If > 0, this gives the period (number of time steps)
-- between logging events.
-- Each log event includes concentrations of all discrete substances.
adjustable RS_log_ts_period = 1

-- maximum number of timesteps for which the concentration log is produced
adjustable RS_log_ts_max = 600

-- whether ill reaction events should be logged (0=no, 1=yes)
adjustable RS_log_ill_events = 0

adjustable SDS_enable

--[[
  Updates the color index of the actual ion.
  Is a callback function called by RS to update simulation state.
--]]
local function RS_updateIonColor(colorindex)
  ion_color = colorindex
end

--[[
  Updates the electrical charge of the actual ion / charged particle.
  Is a callback function called by RS to update simulation state.
--]]
local function RS_updateIonCharge(ionCharge)
  ion_charge = ionCharge
end
-- ---------------------------------------------------

local function getcoords_mm() return ion_px_mm, ion_py_mm, ion_pz_mm end


local first = true
local SDS_initialize_original = SDS.segment.initialize
function SDS.segment.initialize()
  SDS_initialize_original()
  if SDS_enable == 0 then return end

  if first then   first = false
    -- Initialize the kinetic simulation, pass the callback methods to the reaction simulation
    RS:init(RS_updateIonColor, SDS.update_ion_mass, RS_updateIonCharge, RS_log_ill_events ~= 0)
  end -- first

  -- Init ion mass.
  -- SDS.update_ion_mass(ion_mass) is not called here because the ion mass
  -- is initialized by the kinetic simulation script, and SDS.update_ion_mass
  -- is called from RS_simion_interface.lua (RS:addParticle()).

  -- add the actual charged particle / ion to the kinetic simulation
  RS:addParticle(ion_number, ion_color, getcoords_mm())

  -- if this is the first particle, init the chemical species concentration logging
  RS:setupTimestepLogging(RS_log_ts_period, RS_log_ts_max)
end

local SDS_other_actions_original = SDS.segment.other_actions
local last_tof
function SDS.segment.other_actions()
  SDS_other_actions_original()
  if SDS_enable == 0 then return end

  -- Update the spatial particle position in the reaction simulation and
  -- perform the chemical reactions of the particle:
  RS:updateParticlePosition(ion_number, getcoords_mm())
  local E = math.sqrt(ion_dvoltsx_mm^2 + ion_dvoltsy_mm^2 + ion_dvoltsz_mm^2) -- V/mm
  local ke = SDS.ions_local_mfp_mm[ion_number] * E  -- V
  RS:react(ion_number, ke, ion_time_step)

  -- for every (grouped) simulation time step: update the statistics of RS,
  -- and log the actual concentrations / the actual state of RS
  if ion_time_of_flight ~= last_tof then
    last_tof = ion_time_of_flight
    --if sim_grouped == 0 then error 'RS:advanceTimestep only works properly if Grouped flying.  Improve?' end
    RS:advanceTimestep(ion_time_step)
    if RS_log_ts_period > 0 then
      RS:logTimestep()
    end
  end
end

local SDS_terminate_run_original = SDS.segment.terminate_run or function()end
function SDS.segment.terminate_run()
  SDS_terminate_run_original()

  -- print the simulation statistics one time
  if sim_grouped == 0 then
	  print("Warning: Non grouped fly")
  end
  RS:printSimulationStatistics()
end

-- install segments.
SDS.install()


return SDS
