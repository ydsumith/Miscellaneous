# this script fills a region with replica of surfactant

set inputfile "ALFO"
set outputfile "solution"

mol delete all
mol new "$inputfile.psf"
mol addfile "$inputfile.pdb"

set topguy [atomselect top "all"]
set molid [molinfo top get id]
set minimax [measure minmax $topguy] 
set vec [vecsub [lindex $minimax 1] [lindex $minimax 0]] 
set unit_L [lindex $vec 0]
set unit_B [lindex $vec 1]
set unit_H [lindex $vec 2]
puts "initial\nunit_L = $unit_L\nunit_B = $unit_B\nunit_H = $unit_H"
set center [measure center $topguy] 
puts "cellOrigin $center" 


pbc set {13.0 12.0 22.0} -all; 	# set L B H of basix cell
set cell [pbc get -now]; 		# just display pbc info 
$topguy delete 

TopoTools::replicatemol top 5 5 5

mol delete $molid; # delete basis file

set finalstructure [atomselect top "all"]
$finalstructure writepdb "$outputfile.pdb"
package require topotools
topo writelammpsdata "$outputfile.data" full