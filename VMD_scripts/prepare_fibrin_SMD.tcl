
set allatoms [atomselect top all] 

$allatoms set beta 0
set fixedatom [atomselect top "resname CYS and (resid 135 165 193) and x <0"]
$fixedatom set beta 1

$allatoms set occupancy 0
set smdatom [atomselect top "resname CYS and (resid 135 165 193) and x >0"]
$smdatom set occupancy 1

$allatoms writepdb fibrin_mod.pdb

set smdpos [lindex [$smdatom get {x y z}] 0]	 
set fixedpos [lindex [$fixedatom get {x y z}] 0]	 
vecnorm [vecsub $smdpos $fixedpos] 