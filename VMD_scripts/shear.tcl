
proc calcforces {step unique timestep F0 v d0 nx ny nz nx2 ny2 nz2 } {
	#peform only every 10 steps
	if { $step % 10 == 0} {
		# center of the system intially
		set xc 0 
		set yc 0  
		set zc 0
		set RC 15; # (Ang) influence maximum three layers of atoms

		set xproximity 	70; # limit max spread along x as xproximity*2
		set yproximity 	70; # limit max spread along y as yproximity*2
		set zproximity 	70; # limit max spread along z as zproximity*2
		
		set outputfile "carboxy_xcompress";
		
		set curtime 	[expr $step*$timestep]; # t in ps
		set d0vtRC		[expr ($d0 - $v*$curtime - $RC)]

		if { $step == 0} {
			set compressfile [open "$outputfile.log" w];
			puts $compressfile "step\t Natms \t Force(nN) \t upshift(nm) \t downshift(nm)";
			close $compressfile;
		}

		set compressfile [open "$outputfile.log" a]; 

		set counter 0
		set forcemag 0
		set avgupshift 0
		set avgdownshift 0

		while {[nextatom]} {   
			#set mass   [getmass]   ;			# mass (atomic units)
			set rvec [getcoord] ;				# get the atom's coordinates
			set charge [getcharge];
			set roundcharge [expr int(1000*$charge)];			
			
			foreach { x y z } $rvec { break }; 	# get the components of the vector
			set rx [expr $x-$xc];
			set ry [expr $y-$yc];
			set rz [expr $z-$zc];
			
			if { $roundcharge != -833 && $roundcharge != 416} {
				# find the |(ri-rc).n|
				set rircdotn [expr abs($rx)*$nx + abs($ry)*$ny + abs($rz)*$nz];
				
				# check |(ri-rc).n| > (d0 - vt - RC)
				if { $rircdotn > $d0vtRC } {	
					set direction [expr ($nx*$rx/abs($rx) + $ny*$ry/abs($ry) + $nz*$rz/abs($rz))];
					set factor [expr $F0*$direction];
					set forceX [expr -$factor*$nx2]
					set forceY [expr -$factor*$ny2]
					set forceZ [expr -$factor*$nz2]
					set force [expr abs($forceX)+abs($forceY)+abs($forceZ)]
					addforce "$forceX $forceY $forceZ"
					
					set forcemag [expr $forcemag + $force]
					incr counter
					set disp [expr ($x*$nx2 + $y*$ny2 + $z*$nz2)]
					if {$direction > 0} {
						set avgupshift [expr $avgupshift + $disp]
					} else {
						set avgdownshift [expr $avgdownshift + $disp]
					}
				}			
			}
			if { abs($rx) > $xproximity } {
				set forceX [expr -1*$rx/abs($rx)]; 
				set forceY 0
				set forceZ 0
				addforce "$forceX $forceY $forceZ"
			}
			if { abs($ry) > $yproximity } {
				set forceX 0
				set forceY [expr -1*$ry/abs($ry)]
				set forceZ 0
				addforce "$forceX $forceY $forceZ"
			}
			if { abs($rz) > $zproximity } {
				set forceX 0
				set forceY 0
				set forceZ [expr -1*$rz/abs($rz)]
				addforce "$forceX $forceY $forceZ"
			}
		}  
		set forcemag [expr $forcemag * 69.7/1000];
		set avgupshift [expr $avgupshift/10];
		set avgdownshift [expr $avgdownshift/10];
		puts $compressfile  "$step \t $counter \t $forcemag \t $avgupshift \t $avgdownshift"
		close $compressfile
	} else {
		while {[nextatom]} { 
		} 
    }	
}