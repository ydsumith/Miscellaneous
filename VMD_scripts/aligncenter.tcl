proc align {{molid top}} {  
	 set all [atomselect $molid all] 
	 
	 set com [center_of_mass $all]
	 puts "center of mass = $com"
	 $all moveby [vecscale -1.0 $com]; # move all to origin
	 
	 set minmax [measure minmax $all] 
	 set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]] 
	 puts "cellBasisVector1 [lindex $vec 0] 0 0" 
	 puts "cellBasisVector2 0 [lindex $vec 1] 0" 
	 puts "cellBasisVector3 0 0 [lindex $vec 2]" 
	 set center [measure center $all] 
	 puts "cellOrigin $center" 
	 set xoffset [expr [lindex $center 0]*-1]
	 set yoffset [expr [lindex $center 1]*-1]
	 set zoffset [expr [lindex $center 2]*-1]
	 puts "$xoffset $yoffset $zoffset "
	 
	 puts "Enter name of outputfile:"
	 set outputfile [gets stdin]
	 $all writepdb "${outputfile}_centered.pdb"
	 puts "Success!! check ${outputfile}_centered.pdb for details"
	 $all delete 
} 

proc center_of_mass {selection} {
        # some error checking
        if {[$selection num] <= 0} {
                error "center_of_mass: needs a selection with atoms"
        }
        # set the center of mass to 0
        set com [veczero]
        # set the total mass to 0
        set mass 0
        # [$selection get {x y z}] returns the coordinates {x y z} 
        # [$selection get {mass}] returns the masses
        # so the following says "for each pair of {coordinates} and masses,
		#  do the computation ..."
        foreach coord [$selection get {x y z}] m [$selection get mass] {
           # sum of the masses
           set mass [expr $mass + $m]
           # sum up the product of mass and coordinate
           set com [vecadd $com [vecscale $m $coord]]
        }
        # and scale by the inverse of the number of atoms
        if {$mass == 0} {
                error "center_of_mass: total mass is zero"
        }
        # The "1.0" can't be "1", since otherwise integer division is done
        return [vecscale [expr 1.0/$mass] $com]
}