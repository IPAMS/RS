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
  RS_Particle.lua

  Data-Type class for chemical reaction dynamic simulation extension to the SDS algorithm (RS):
  This class implements a simulated particle (which are stored in the particle table of the RS extension)

  @author Walter Wissdorf
  @version: 0.4.4
--]]
local Particle = {} --define namespace

	--[[
	 Constructor (according to basic pattern described in "Programming in Lua"):
	 Constructs a new Particle instance

	 @param: substance = the substance to which the new particle belongs to
	 @param: x,y,z = x,y,z world position of the particle (optional, is initialized randomly if ommitted)

	 @return: a new instance of the Particle class
	--]]
	function Particle:new(substance, x,y,z)
		local new_instance = {} --create a new instance

		--set class metatable and sets the base class (Substance) as index (setups inheritance)
		setmetatable(new_instance, self)
		self.__index = self --set index

		new_instance.s = substance
		new_instance.x = x or 0
		new_instance.y = y or 0
		new_instance.z = z or 0

		new_instance.next = nil -- link to the next node in the linked list
		new_instance.before = nil -- link to the preceding node in the linked list

		return new_instance
	end
--- end of class definition "Particle" --- --- --- --- --- --- --- --- --- --- --- --- --- ---

return Particle --add the particle namespace to the RS namespace
