proc get_strain {} { 

 set all [atomselect top "all not water"] 

 set nf [molinfo top get numframes]
 
 set outfile [open "oxygenated_strain.out" w]
 puts $outfile "frame L B H xc yc zc"
 
 for {set i 0} {$i < $nf } {incr i} {
	 $all frame $i
	 set minmax [measure minmax $all] 
	 set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]] 
	 set center [measure center $all] 
	 
	 puts $outfile "$i [lindex $vec 0] [lindex $vec 1] [lindex $vec 2] $center" 
	 #puts "cellBasisVector2 0 [lindex $vec 1] 0" 
	 #puts "cellBasisVector3 0 0 [lindex $vec 2]" 
	 #puts "cellOrigin $center" 
 }
 $all delete 
 close $outfile
} 

