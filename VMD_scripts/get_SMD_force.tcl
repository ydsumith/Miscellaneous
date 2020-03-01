set nx 0.02529
set ny 0.99585 
set nz -0.087379

# fs
#set timestep 2
# Ang/fs
#set SMDvel 0.0001

puts "using normal vectors: $nx $ny $nz"

### Loop over all lines of the log file
set file [open out_2314321.pbs.scm.log r]
set output [open out_force.dat w]

puts $output "step Disp (nm) Force(nN)"

while { [gets $file line] != -1 } {

### Determine if a line contains SMD output. If so, write the
### timestep followed by f(dot)n to the output file
  if {[lindex $line 0] == "SMD"} {
	if {[lindex $line 1] == "0"} {
			set x0 [lindex $line 2]
			set y0 [lindex $line 3]
			set z0 [lindex $line 4]
		}
	 set displacement [expr sqrt( ([lindex $line 2]-$x0 )**2 + ([lindex $line 3]-$y0 )**2 + ([lindex $line 4]-$z0 )**2)/10 ]
	 set force [expr ($nx*[lindex $line 5]+ $ny*[lindex $line 6] + $nz*[lindex $line 7])/1000 ]
	 # old stuff
	 #puts $output "[lindex $line 1] $displacement [expr $nx*[lindex $line 5]+ $ny*[lindex $line 6] + $nz*[lindex $line 7]] [expr [lindex $line 1]*$SMDvel/10/$timestep]"
	 puts $output "[lindex $line 1] [format "%.3f" $displacement] [format "%.3f" $force]"
    }
}

### Close the log file and the output .dat file
close $file
close $output