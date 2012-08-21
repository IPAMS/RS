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
  RS_ParticleList

  Data-Type class for chemical reaction dynamic simulation extension to the
  SDS algorithm (RS):

  This class implements a simple bidirectional linked list for the
  storage of the simulated particles.
  This implementation has has O(1) costs for element deletion, while the
  original Lua implementation of table operations has O(n) costs for element
  deletion.
  Because of the probably very high number of element deletion operations
  in the RS code, the utilization of a simple Lua table as particle list
  is not feasible.

  @author Walter Wissdorf
  @version: 0.4.4
--]]

local ParticleList = {} -- define namespace


	--[[
	 Constructor (according to basic pattern in "Programming in Lua"):
	 Constructs a new ParticleList instance

	 @return: a new instance of the ParticleList class
	--]]
	function ParticleList:new()
		local new_instance = {} -- create a new instance

		-- set class metatable and sets the base class (Substance) as
		-- index (setups inheritance)
		setmetatable(new_instance, self)
		self.__index = self -- set index

		new_instance.start = nil -- the first node of the linked list
		new_instance.nNodes = 0 -- the number of nodes in the linked list
		return new_instance
	end


	--[[
	 Insert:
	 Inserts a node to the start of the linked list / particle table

	 @param: node = the new list node to add to the linked list / particle
	         table (probably an instance of RS.Particle)
	--]]
	function ParticleList:insert(node)

		-- --insert it at the start position:
		-- back link the old start node (if it exists) to the new start node:
		if self.start then
			self.start.before = node
		end

		-- link the old start node as follower to the new start node:
		node.next = self.start

		-- link the list start to the new start node
		self.start = node

		self.nNodes = self.nNodes + 1
	end


	--[[
	 Remove:
	 Removes an element from the linked list / particle table

	 @param: node = the element to remove (probably an instance of RS.Particle
	--]]
	function ParticleList:remove(node)

		-- link and backlink the before node to the next node
		if node.before then
			node.before.next = node.next
		else
			-- it is the first one => update start node
			self.start = node.next
		end

		if node.next then
			node.next.before = node.before
		end

		self.nNodes = self.nNodes - 1
	end


	--[[
	 Get List Iterator:
	 returns an iterator for the linked list / particle table (which can be
	 used in "for" loops
	 to iterate trough the linked list / particle table)

	 @return: a new iterator for the linked list / particle table
	--]]
	function ParticleList:getListIterator()
		local actualNode = self.start
		return
			function ()
				if actualNode then
					local result = actualNode
					actualNode = actualNode.next
					return result
				end
			end
	end


	--[[
	 Get Number of Nodes:
     returns the number of nodes in the linked list / particle table

     @return: the number of nodes
	--]]
	function ParticleList:getNumberOfNodes()
		return self.nNodes
	end


	--[[
	 Print List:
     prints the list for debugging purposes
	--]]
	function ParticleList:printList()
		for particle in self:getListIterator() do
			print(particle.s:getName(), particle.x, particle.y, particle.z)
		end
	end

--- end of class definition "ParticleList" --- --- --- --- --- --- --- ---

return ParticleList
