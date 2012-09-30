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
  RS_Reaction

  Data-Type class for chemical reaction dynamic simulation extension to the SDS algorithm (RS):

  This class implements a chemical reaction which consists of a set of
  educts, a set of products, stochiometric factors to the individual substances
  and a rate constant.

  The reaction is assumed to be invariant, therefore there are no methods to change the reaction parameters
  or the educts / products tables

  @author Walter Wissdorf
  @version: 0.4.4
--]]


local Reaction = {} --define namespace


	--[[
	 Constructor (according to basic pattern described in "Programming in Lua"):
	 Constructs a new Reaction instance.

	 The attending educts and products are passed to the Reaction constructor in dictionary data structures,
	 which have the following general layout:

	   {substance_1 = stoechiometricFactor_1, substance_2 = stoechiometricFactor_2, ...}

	 Example: for the reaction A + B => 2C the dictionaries / tables would be:
	   educts=  {A = 1, B = 1}
	   products = {C = 2}

	  (the substances are instances of the "RS_Substance" class)

	 @param: educts = educts dictionary / table
	 @param: products = products dicttionary / table
	 @param: rateConstant = the rate constant of the reaction (base time unit is microseconds)
	 @param: e_A = (optional) activation energy of the reaction

	 @return: a new RS_Reaction instance
	--]]
	function Reaction:new(educts, products, rateConstant, e_A)
		local ni = {} --create a new instance

		--set class metatable and sets the base class (Substance) as index (setups inheritance)
		setmetatable(ni, self)
		self.__index = self --set index

		-- init member variables:
		ni.educts = educts or {}
		ni.products = products or {}
		ni.discreteEducts = {} -- the table of discrete educts
		ni.discreteProductsTable = {} -- a simple list of the discrete products

		ni.rateConstant = rateConstant or -1.0 --the rate constant of the reaction
		ni.activationEnergy = e_A or nil
		ni.staticProb = ni.rateConstant --the static (invariant part of) reaction probability
										--(product of all isotropic educt concentrations and the reaction rateConstant)
		ni.indep = false -- is set to true if the reaction is independent => only dependent on one discrete educt

		--init the discrete products table
		for product, factor in pairs(ni.products) do -- actually only discrete products are allowed by the config file parser,
			                                        -- but the reaction implementation should be general
			if product:getType() == "discrete" then --if the product is discrete (actually enforced by the config file parser)
				-- ni.discreteProducts[product] = factor
				for i=1, factor do
					table.insert(ni.discreteProductsTable, product)
				end
			end
		end


		-- search discrete educts, calculate static reaction probability
		local nDiscrete = 0 --number of found discrete educt reaction partners
		for educt, factor in pairs(ni.educts) do
			if educt:getType() == 'discrete' then
				nDiscrete = nDiscrete + factor
				ni.discreteEducts[educt] = factor
			elseif educt:getType() == 'isotropic' then
				ni.staticProb = ni.staticProb * (educt:getStaticConcentration() ^factor)
			end
		end

		-- if we have only one discrete educt, the whole reaction is classified as independent
		-- (from other discrete species)
		if nDiscrete == 1 then
			ni.indep = true
		end

		return ni
	end


	--[[
	 Set Educts:
     sets a new educts table (definition see above in the comments to the class constructor)

     @param: et = the new educts table
	--]]
	function Reaction:setEducts(et)
		self.educts = et
	end


	--[[
	 Get discrete Educts:
	 returns the table of discrete educts (educts which are represented by simulation particles) of the reaction

	 @return: the discrete educts table of the reaction (table definition according to the class constructor)
	--]]
	function Reaction:getDiscreteEducts()
		return self.discreteEducts
	end


	--[[
 	 Get Products:
 	 gets the products table of the reaction

 	 @return: the products table of the reaction (table definition according to the class constructor)
	--]]
	function Reaction:getProducts()
		return self.products
	end


	--[[
	 Get Discrete Products Table:
	 gets the table (a linear list, one table entry per discrete product particle)
	 of the discrete products of the reaction

	 @return: the table of discrete products of the reaction
	--]]
	function Reaction:getDiscreteProductsTable()
		return self.discreteProductsTable
	end


	--[[
	 Get Rate Constant:
	 gets the rate constant of the Reaction

	 @return: the rate constant of the Reaction
	--]]
	function Reaction:getRateConstant()
		return self.rateConstant
	end


	--[[
	 Get Activation Energy:
	 gets the activation energy of the Reaction

	 @return: the activation energy of the Reaction
	--]]
	function Reaction:getActivationEnergy()
		return self.activationEnergy
	end


	--[[
	 Is Independent:
	 gets the independence state of the Reaction

	 @return: true if the Reaction is independent
	--]]
	function Reaction:isIndependent()
		return self.indep
	end


	--[[
	 Get Static Prob:
	 gets the static reaction probability of the Reaction.
	 The static reaction probability is the invariant part of the reaction probability:
	 The product of all isotropic educt concentrations and the reaction rateConstant

	 @return: the static reaction probability
	--]]
	function Reaction:getStaticProb()
		return self.staticProb
	end

    --[[
     Gets a string representation of the reaction

     @return: a string which describes the whole reaction system
    --]]
	function Reaction:getStringRepresentation()
		local str =
				Reaction:getSubstTableString(self.educts) ..
				" => "..
				Reaction:getSubstTableString(self.products) ..
				" | k ="..self.rateConstant ..
				" | static Prob=" .. self.staticProb

		return str
	end

	--[[
	 Print Reaction:
	 prints a summary of the Reaction to the console / log file
	--]]
	function Reaction:printReaction()
		print(self:getStringRepresentation() )
		-- additional information, was used for debugging:
		-- print("Discrete Educts: "..Reaction:getSubstTableString(self.discreteEducts))
		-- print("is independent:"..tostring(self:isIndependent()))

	end


	--[[
     Get Substance Table String:
     creates a formatted string with the contents of the given substance table (help-method for the printing / logging
     methods)

     @param: r = the substance table to construct the formatted string for
     @return: the formatted string representation of the substance table "r"
	--]]
	function Reaction:getSubstTableString(r)
		local str = ""
		local first = true
		for subst, fact in pairs(r) do
			--the first entry needs no preceeing "+"
			if first then
				first = false
			else
				str = str.." + "
			end

			str = str..fact.." "..subst:getName()
		end

		return str
	end


--- end of class definition "Reaction" --- --- --- --- --- --- --- --- --- --- --- --- --- ---

return Reaction
