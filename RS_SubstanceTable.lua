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
  RS_SubstanceTable

  Datatype class for kinetic simulation extension to the SDS algorithm

  This class implements a data structure in which all substances of the
  kinetic simulation are stored in a convenient manner.

  The substances are stored in a numerically indexed array, but can also
  be retrieved by their names.
  In addition, the discrete substances are addressable separately.

  @author Walter Wissdorf
  @version: 0.4.4
--]]


local ST = {} --define namespace


	--[[
	 Constructor (according to basic pattern described in "Programming in Lua"):
	 constructs a new SubstanceTable instance

	 @return: a new SubstanceTable instance
	--]]
	function ST:new()

		local new_instance = {} --create a new instance

		--set class metatable and sets the base class (SubstanceTable) as index (setups inheritance)
		setmetatable(new_instance, self)
		self.__index = self --set index

		--set member fields
		new_instance.s = {} -- the (indexed) substance table
		new_instance.m = {} -- map from substanceName => substance index in s
		new_instance.ds= {} -- list (indexed) from discrete substance => substance index in s
		return new_instance
	end


	--[[
	 Add Substance:
	 add a substance to the SubstanceTable. The substances are identifiable by their names. If there is already a substance
	 with a given name in the SubstanceTable, this existing substance is silently overwritten.

	 @param: substance = a substance (an instance of the Substance class)
	 @param: substanceName = a name / identifier of the substance
	--]]
	function ST:addSubstance(substance, substanceName)
		if self.m[substanceName] then
			local buf = self.s[self.m[substanceName]] -- the old substance entry in the substance table

			-- if the substance is a discrete substance and the old was a discrete => do nothing
			-- (discrete substance table is up to date after update of the substance in the substance table)

			-- if the substance is a non discrete substance, and the old was discrete => remove substance from discrete table
			if buf:getType() == 'discrete' and substance:getType() ~= 'discrete' then
				-- remove from discrete table
				table.remove(self.ds, self.m[substanceName])
			end

			-- if the substance is discrete and the old was non discrete => add substance to discrete substances table
			if buf:getType() ~= 'discrete' and substance:getType() == 'discrete' then
				-- add to discrete table
				table.insert(self.ds, self.m[substanceName])
			end

			-- the substance identifier is already existing, update the substance silently in the substance table...
			self.s[self.m[substanceName]] = substance

		else
			-- add a new substance
			table.insert(self.s, substance)
			self.m[substanceName] = #self.s

			-- if substance is a discrete substance: add the substance index to discrete substance table
			if substance:getType() == 'discrete' then
				table.insert(self.ds, #self.s)
			end
		end
	end


	--[[
     Get Substance By Name:
     gets a Substance by its name / identifier

     @param: substanceName = the name of the Substance

     @return: the substance with the given substance name
	--]]
	function ST:getSubstanceByName(substanceName)
		return self.s[self.m[substanceName]]
	end


	--[[
	 Get Substance Name Iterator:
	 gets an iterator which iterates trough all substances stored in the SubstanceTable.
	 The iterator returns the Substance names and the according Substances

	 @return: an iterator of the substance names and the substances
	--]]
	function ST:getSubstanceNameIterator()
		local iter, map, substName = pairs(self.m)
		local index
		return
			function ()
				substName, index = iter(map, substName)
				return substName, self.s[index]
			end
	end


	--[[
     Get Substance By Index:
     gets a substance by its numerical index

     @param: index = the numerical index of the substance

     @return: the Substance with the given numerical index
	--]]
	function ST:getSubstanceByIndex(index)
		return self.s[index]
    end


	--[[
	 Get Substance Index:
	 gets the index of a given substance

	 @param: substance = the substance to return the index for

	 @return: the index of the given substance
	--]]
    function ST:getSubstanceIndex(substance)
        return self.m[substance:getName()]
    end


	--[[
	 Get Discrete Substance By Index:
	 gets a discrete substance by its index

	 @param index = The numerical index (in the set of discrete substances) of the discrete substance to return

	 @return: the discrete substance with the given index
	--]]
	function ST:getDiscreteSubstanceByIndex(index)
		return self.s[self.ds[index]]
	end


	--[[
	 Get Substance Iterator:
	 gets an iterator of all substances in a SubstanceTable
	 The iterator returns the numerical index and the according substance.

	 @return: an iterator of all substances
	--]]
	function ST:getSubstanceIterator()
		local i = 0

		return
			function ()
				i = i+1
				if i <= #self.s then
					return i, self.s[i]
				end
			end
	end


	--[[
	 Get Discrete Substance Iterator:
	 gets an iterator of all discrete substances in a SubstanceTable
	 The iterator returns the numerical index (in the set of discrete substances) and the according substance.

	 @return: an iterator of all discrete substances
	--]]
	function ST:getDiscreteSubstanceIterator()
		local i = 0

		return
			function ()
				i = i+1
				if i <= #self.ds then
					return i, self.s[self.ds[i]]
				end
			end
	end


	--[[
	 Get Number of Substances

	 @return: The number of all substances in the SubstanceTable
	--]]
	function ST:getNumberOfSubstances()
		return #self.s
	end


	--[[
	 Get number of discrete substances

	 @return: The number of discrete substances in the SubstanceTable
	--]]
	function ST:getNumberOfDiscreteSubstances()
		return #self.ds
	end


	--[[
	 Print Substance Table:
	 prints the substance table in a readable manner to the console / log file
	--]]
	function ST:printSubstanceTable()
		for i, s in self:getSubstanceIterator() do
			print(i.."-->"..s:getName().." -->"..s:getType())
		end
	end


--- end of class definition "SubstanceTable" --- --- --- --- --- --- --- --- --- --- --- --- --- ---

return ST
