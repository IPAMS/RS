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
  RS Config file parser

  This file contains a config file parser for the config files for SDS
  (statistical diffusion simulation) with monte carlo chemical reactions
  simulation.

  This class implements a table in which all substances of the kinetic
  simulation are stored in a convenient way.

  @author Walter Wissdorf
  @version: 0.4.4
--]]

local SubstanceTable = require 'RS_SubstanceTable'
local Substance      = require 'RS_Substance'
local Reaction       = require 'RS_Reaction'


local P = {} -- create a local namespace

-- define some constants:
-- -- the key strings for the sections of the config file:
P.keySubs = "[substances]"  -- the substances part of the config file
P.keyReac = "[reactions]"   -- the reactions part of the config file


--[[
 Open File:
 opens and reads the contents of a given config file
 @param filename = the filename of the file to open
--]]
function P:openFile(filename)
	self.f = assert(io.open(filename, "r")) -- self.f = the file handler 
	self.t = self.f:read("*all")            -- self.t = the input text
	self.f:close()
end


--[[
 Print Input:
 prints the input of the config file parser (the contents of the opened file)
--]]
function P:printInput()
	print(self.t)
end


--[[
 Split:
 Splits a string based on a separator string or pattern
 (implementation based on implementation found in lua wiki)

 @param: str = the input string
 @param: inSplitPattern = the pattern which is used as delimiter for the split operation
 @param (optional): outResults = table to be filled with the results
 @return: an array of pieces of the string
 @return: an array of the matched keys
--]]
function P:split(str, inSplitPattern, outResults )
    local outMatches = {}
    local outResults = outResults or {}

	local theStart = 1
	local theSplitStart, theSplitEnd = string.find(str, inSplitPattern, theStart)
	table.insert(outMatches, "")
	while theSplitStart do
		table.insert(outResults, string.sub(str, theStart, theSplitStart-1 ))
		table.insert(outMatches, string.sub(str, theSplitStart, theSplitEnd))
		theStart = theSplitEnd + 1
		theSplitStart, theSplitEnd = string.find(str, inSplitPattern, theStart)
	end
	table.insert(outResults, string.sub(str, theStart))
	return outResults, outMatches
end


--[[
 Parse Substances:
 parses the substance part of the configuration file
 @param: inp = the input string to parse (substances section of a config file in a string)
 @return: the parsed substances packed in a RS_SubstanceTable object
--]]
function P:parseSubstances(inp)
	local result = SubstanceTable:new()
	
	print("parsing substances...")
	local k = 1
	for line in string.gmatch(inp, "[%w_-]+[ \t]*%a+[ \t]*[%deE%.%_%-%+]*[ \t]*[%deE._-]*[ \t]*\r*\n")
	do --match all substance lines
		--extract parts of a substance line:
		local name, sType, optNumber1, optNumber2
			= string.gmatch(line, "([%w_-]+)[ \t]*(%a+)[ \t]*([%deE%.%_%-%+]*)[ \t]*([%deE._-]*)")()
			--gmatch gives an iterator, but we have only one substance per line
		--check / process input:
		local subst = Substance:new(name, sType)
		if sType == 'isotropic' and optNumber1 ~= '' then
			subst:setStaticConcentration(tonumber(optNumber1))
		elseif sType == 'isotropic' then
			print ("Isotropic substance without concentration value found, the conc. is treated as 0")
		end

		if sType == 'discrete' and optNumber1 ~= '' and optNumber2 ~= '' then
			subst:setMass(tonumber(optNumber1))
			subst:setCharge(tonumber(optNumber2))
		elseif sType =='discrete' then
			error ("Discrete substance and mass and / or charge value missing")
		end
		
		result:addSubstance(subst, name) --add the parsed substance to the substance table
		k = k+1
	end
	print("complete...")
	return result
end


--[[
 Parse Reactions:
 parses the reactions part of a given configuration file
 @param: inp = the input string to parse (reaction section of a config file in a string)
 @param: st = a substance table with the substances occurring in the reactions to parse
 @param: rateConstConvFactor = a factor by which the given rate constants are
   divided to allow conversion of the basic time unit (typical rate constants
   are given with an basic time unit of seconds, but the basic time unit in
   SIMION is microsecond)
 @return: an array with the parsed reactions, packed in RS_Reaction objects
--]]
function P:parseReactions(inp, st, rateConstConvFactor)
	print("parsing reactions...")
	local result = {}


	--for every line of the input...
	local k = 0
	for line in string.gmatch(inp,"[^\r\n]+") do 
		line = line:gsub("[ \t]+","") --remove whitespace (spaces are unnecessary in our case)
		-- split the line at the ";" symbol (delimiter between reaction definition,
		-- rate constant and E_a (optional activation energy) )
		local parts = P:split(line,";")
		if #parts ~= 2 and #parts ~= 3 then
			error("wrong format in reaction line "..k)
		end
		
		local formulaStr = parts[1] --the formula string is the first value (left of the ";")
		local constantStr = parts[2] -- the rate constant is second value (after the first ";")
		local e_AStr = parts[3] -- the string of the optional activation energy (nil if the activation energy is not set)
		-- -----------------
		
		-- now parse the formula string 
		local educts = {}
		local products = {}
		local reacParts = P:split(formulaStr,"=>") -- the "=>" symbol is the delimiter between educt and product side
		local edStr = reacParts[1]
		local prStr = reacParts[2]

		-- parse the educt side
		for i,edStr in ipairs(P:split(edStr,"+")) do
			local mult,substName = P:parseReactionPartnerString(edStr)
			local subst = st:getSubstanceByName(substName) --the substance object with the parsed substance name
			if subst then
				if educts[subst] then -- the substance is already present on the educt side
					educts[subst]= educts[subst] + tonumber(mult) -- the present substance multiplier is increased
				else
					educts[subst]= tonumber(mult) -- new substance entry on the educt side
				end

			else
				error("Educt "..substName.." is not in substance table")
			end
		end
				
		-- parse the product side
		for i,prStr in ipairs(P:split(prStr,"+")) do
			local mult,substName = P:parseReactionPartnerString(prStr)
			local subst = st:getSubstanceByName(substName) --the substance object with the parsed substance name
			if subst then
				if products[subst] then -- the substance is already present on the product side
					products[subst] = products[subst] + tonumber(mult) --increase substance multiplier
				else
					products[subst] = tonumber(mult) -- new substance entry on the product side
				end
			else
				error("Product "..substName.." is not in substance table")
			end
		end
		-- ---------------------------
		
		-- now parse the rate constant string
		-- (the rate constant string is typically given with a basic time unit of seconds)
		-- to allow conversion of the time unit, divide the rate constant by a conversion factor
		local rateConstant = tonumber(constantStr) / rateConstConvFactor
		-- ---------------------------

		-- now parse the activation energy string
		local e_A
		if e_AStr then
			e_A  = tonumber(e_AStr)
		end
		-- ---------------------------
		
		local r = Reaction:new(educts,products,rateConstant,e_A)
		table.insert(result,r)
		
		k = k+1
	end
	print("complete...")
	
	return result
end


--[[
 Parse Reaction Partner String:
 parses a reaction partner string representation, which comprises of a
 substance name with a stochiometric multiplier
 the multiplier could be omitted and is assumed to be 1 in this case
 @param: str = the input string to parse
 @return: the parsed stochiometric multiplier
 @return: the parsed substance name
--]]
function P:parseReactionPartnerString(str)
	local mult,substName
	if str:find("%d") == 1 then -- if first char of str is a number
		mult, substName = str:match("(%d+)([%w_-]+)") --parse with multiplicator
	else
		mult = 1
		substName = str:match("([%w_-]+)")
	end

	return mult,substName
end


--[[
 Check Simulation Configuration:
 Checks if a given simulation configuration (consisting
 of a table of substances (of type "Substance", see RS_Substance.lua,
 packed in a RS_SubstanceTable, see RS_SubstanceTable.lua))
 and a list of reactions (of type "Reaction", see RS_Reaction.lua also)

 Currently the following tests are performed:
  -"field" and "isotropic" substances are only allowed to be on the left side of a reaction

 @param: substances = the substances in the configuration, packed in a RS_SubstanceTable object
 @param: reactions = array of reactions (RS_Reaction objects)
--]]
function P:checkSimulationConfiguration(substances,reactions)
	for i,r in ipairs(reactions) do 
		-- check if field and isotropic is only on the left side of the reactions
		local flag = false

		for rp,_ in pairs(r:getProducts()) do --iterate trough all products of the reaction
			if rp:getType() == Substance.types["isotropic"] or
			   rp:getType() == Substance.types["field"]
			then
				flag = true
			end
		end

		if flag == true then
			--print a warning message:
			print("Warning: 'isotropic' or 'field' substance on product side of reaction "..
			      i..", those concentrations are considered static!")
		end
		-- ----------------------------------------------------------------------	
	end
end


--[[
 Parse:
 parses an opened configuration file

 @param: rateConstConvFactor = a factor by which the parsed rate constants
         are divided to allow conversion of the basic time unit
 @return: the parsed substance table (as RS_SubstanceTable containing RS_Substance objects)
 @return: an array of reactions (array of RS_Reaction objects)
--]]
function P:parse(rateConstantConvFactor)
	local parts,matches = P:split(self.t,"%[%w+%]")
	local substanceTable = P:parseSubstances(parts[2])
	local reactions = P:parseReactions(parts[3],substanceTable,rateConstantConvFactor)
	
	P:checkSimulationConfiguration(substanceTable,reactions)
	return substanceTable, reactions
end


return P
