
proc calcforces {step unique timestep F0 v d0 nx ny nz} {
	#peform only every 10 steps
	if { $step % 10 == 0} {
		# center of the system intially
		set xc 0 
		set yc 0  
		set zc 0
		set RC 8; # (Ang) influence maximum three layers of atoms

		set xproximity 	70; # limit max spread along x as xproximity*2
		set yproximity 	70; # limit max spread along y as yproximity*2
		set zproximity 	70; # limit max spread along z as zproximity*2
		
		set outputfile "compressorlog";
		
		set curtime 	[expr $step*$timestep]; # t in ps
		set d0vtRC		[expr ($d0 - $v*$curtime - $RC)]

		if { $step == 0} {
			set compressfile [open "$outputfile.log" w];
			puts $compressfile "step\t Natms \t Force";
			close $compressfile;
		}

		set compressfile [open "$outputfile.log" a]; 

		set counter 0
		set forcemag 0

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
					set factor [expr ($rircdotn-$d0vtRC)*($rircdotn-$d0vtRC)*$F0];
					set forceX [expr -$factor*$nx*$rx/abs($rx)]
					set forceY [expr -$factor*$ny*$ry/abs($ry)]
					set forceZ [expr -$factor*$nz*$rz/abs($rz)]
					set force [expr abs($forceX)+abs($forceY)+abs($forceZ)]
					addforce "$forceX $forceY $forceZ"
					
					set forcemag [expr $forcemag + $force]
					incr counter
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
		puts $compressfile  "$step \t $counter \t $forcemag"
		close $compressfile
	} else {
		while {[nextatom]} { 
		} 
    }	
}