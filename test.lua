-- test.lua - SIMION workbench user program illustrating use of
-- the chemical reaction model extension to SDS (SIMON SDS / RS)

simion.workbench_program()

-- Load SDS / RS user program: 
local SDS = simion.import("RS_collision_sds.lua")

function SDS.init()
  --[optionally:] Plot gas flow. (uncomment the following lines to plot the gas flow,
  --contourlib81.lua from the SIMION examples must be present)
  --local CON = simion.import 'contourlib81.lua'  -- [3][4]
  --CON.plot{func=SDS.velocity,  npoints=20, z=0, mark=true} --[2]
end
