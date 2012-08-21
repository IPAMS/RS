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
  RS_Substance

  Data-Type class for chemical reaction dynamic simulation extension to the SDS algorithm (RS):

  This class implements a substance which is one of the partners in a simulated chemical reaction

  There are three types of substances:
	 - isotropic = neutral isotropic distributed
	 - discrete = substances which are modeled by discrete simulation particles
	 - field = non isotropic neutral substances (not fully implemented right now)

  @author Walter Wissdorf
  @version: 0.4.4
--]]

local Substance = {} --define namespace

	-- "static" field: possible types of substances
	Substance.types = {isotropic = 'isotropic', discrete = 'discrete', field= 'field'}


	--[[
	 Constructor (according to basic pattern described in "Programming in Lua"):
	 constructs a new Substance instance

	 @param: name = the name of the substance
	 @param: type = the type of the substance (one of the fields in Substance.types)

	 @return: a new Substance instance
	--]]
	function Substance:new(name, type)

		local t = self.types[type] --check if "type" an allowed type

		if t ~= nil then --create the new instance
			local new_instance = {} --create a new instance

			--set class metatable and sets the base class (Substance) as index (setups inheritance)
			setmetatable(new_instance, self)
			self.__index = self --set index

			--set member fields (probably self-explanatory)
			new_instance.name = name
			new_instance.type = t
			new_instance.static_concentration = 0    -- for isotropic stubstances
			new_instance.mass = 0   --for discrete substances
			new_instance.charge = 0 --for discrete substances

			return new_instance
		else
			error('illegal substance type')
		end
	end


	--[[
     Set Static Concentration:
     sets the static concentration of this substance (only relevant for isotropic substances)

     @param: newConcentration = the new concentration of this (isotropic) substance
	--]]
	function Substance:setStaticConcentration(newConcentration)
			self.static_concentration = newConcentration
	end


	--[[
     Get Static Concentration:
     gets the static concentration of this substance (only defined for isotropic substances)

     @return: the static concentration of this substance
	--]]
	function Substance:getStaticConcentration()
		return self.static_concentration
	end


	--[[
	 Get Name:
	 gets the name of the substance

	 @return: the name of the substance
	--]]
	function Substance:getName()
		return self.name
	end


	--[[
	 Get Type:
	 gets the type of the substance

	 @return: the type of the substance (one of the fields in Substance.types)
	--]]
	function Substance:getType()
		return self.type
	end


	--[[
	 Set Mass:
	 sets the mass of the substance (the mass is only relevant for discrete substances)

	 @param: newMass = the new mass of this (discrete) substance
	--]]
	function Substance:setMass(newMass)
		self.mass = newMass
	end


	--[[
	 Get Mass:
	 gets the mass of the substance (the mass is only relevant and probably only defined for discrete substances)

	 @return: the mass of the discrete substance
	--]]
	function Substance:getMass()
		return self.mass
	end


	--[[
	 Set Charge:
	 sets the charge of the substance (the charge is only relevant for discrete substances)

	 @param: charge = the new charge of the (discrete) substance
	--]]
	function Substance:setCharge(charge)
		self.charge = charge
	end


	--[[
	 Get Charge:
	 gets the charge of the substance (the charge is only defined for discrete substances)

	 @return: the charge of the (discrete) substance
	--]]
	function Substance:getCharge()
		return self.charge
	end


--- end of class definition "Substance" --- --- --- --- --- --- --- --- --- --- --- --- --- ---

return Substance
