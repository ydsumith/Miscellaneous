set inputfile glycan
set seg1 O7
set seg2 O6
set seg3 O5
set seg4 O4
#------------dont modify below---------
draw delete all
package require psfgen

mol new "$inputfile.psf"
mol addfile "$inputfile.pdb"

set molid [molinfo top get id]
mol delrep 0 $molid
mol color SegName
mol representation VDW 1.000000 12.000000
mol selection name FE and segname $seg1
mol material Opaque
mol addrep $molid
mol color SegName
mol representation VDW 1.000000 12.000000
mol selection name FE and segname $seg2
mol material Opaque
mol addrep $molid
mol color SegName
mol representation VDW 1.000000 12.000000
mol selection name FE and segname $seg3
mol material Opaque
mol addrep $molid
mol color SegName
mol representation VDW 1.000000 12.000000
mol selection name FE and segname $seg4
mol material Opaque
mol addrep $molid
mol color name
mol representation CPK 1.4 0.5 8 8
mol selection resname HEME
mol material Opaque
mol addrep $molid

set all [atomselect top "all"] 
set FE1 [atomselect top "name FE and segname $seg1"] 
set FE2 [atomselect top "name FE and segname $seg2"] 
set FE3 [atomselect top "name FE and segname $seg3"] 
set FE4 [atomselect top "name FE and segname $seg4"] 

set FE1xyz [measure center $FE1 weight mass] 
set FE2xyz [measure center $FE2 weight mass] 
set FE3xyz [measure center $FE3 weight mass] 
set FE4xyz [measure center $FE4 weight mass] 

draw color blue
draw line $FE1xyz $FE2xyz
draw line $FE1xyz $FE4xyz

set FE12 [vecnorm [vecsub $FE1xyz $FE2xyz]]; # normal along 1-2
set FE14 [vecnorm [vecsub $FE1xyz $FE4xyz]]; # normal along 1-4
set normal [vecnorm [veccross $FE12 $FE14]];

set matrix [transvecinv $normal]
$all move $matrix 
set FE1xyz [measure center $FE1 weight mass]; #measure again for updated
$all moveby [vecscale -1.0 $FE1xyz]; # move all with FE1 to origin
set FE1xyz [measure center $FE1 weight mass] 
set FE2xyz [measure center $FE2 weight mass] 
set FE3xyz [measure center $FE3 weight mass] 
set FE4xyz [measure center $FE4 weight mass] 

draw color red
draw line $FE1xyz $FE2xyz
draw line $FE1xyz $FE4xyz

set FE12 [vecnorm [vecsub $FE1xyz $FE2xyz]]; # normal along 1-2
set FE14 [vecnorm [vecsub $FE1xyz $FE4xyz]]; # normal along 1-4
set costheta [expr "[vecdot $FE14 {0 0 1}]"]
set theta [expr "acos($costheta)"]
set matrix [transaxis x $theta rad]
$all move $matrix

set FE1xyz [measure center $FE1 weight mass] 
set FE2xyz [measure center $FE2 weight mass] 
set FE3xyz [measure center $FE3 weight mass] 
set FE4xyz [measure center $FE4 weight mass] 

draw color green
draw line $FE1xyz $FE2xyz
draw line $FE1xyz $FE4xyz
$all writepdb "${inputfile}_rot.pdb"
puts "final file written to ${inputfile}_rot.pdb"