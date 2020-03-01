# Written by Sumith YD, 9/8/2016

# Make a L x B platinum plate with single layer
# Write to $output.psf and $output.pdb

puts "usage: \nmake_plate L(A) B(A) seg(2chars) x_shift(A) y_shift(A) z_shift(A) output(string)"
proc make_plate {L B seg x_shift y_shift z_shift output} {

  # spacing
  set a 1.50
  set half_a [expr $a/2.0]
  set m [expr int($L/$half_a)]
  set n [expr int($B/$half_a)]
  set natom [expr ($m*$n)/2]
  set segcounter 1
  
  set fd_top [open temp_platinum.top w]
  puts $fd_top "0 0"			
  puts $fd_top "MASS     1 PT   195.08400 P"
  #puts $fd_top "AUTO ANGLES DIHE"
  puts $fd_top "RESI PLA          0.00"
  
  set fd_pgn [open temp_platinum.pgn w]
  puts $fd_pgn "topology temp_platinum.top"
  puts $fd_pgn "segment $seg$segcounter {residue 1 PLA}"
  
  # Define the atom names for each platinum
  set ind 0
  for {set i 0} {$i<$natom} {incr i} {
	if {$ind == 1000} {
		set ind 0
	}
	if {$ind > 99} {
       set name($i,0) P$ind
    } elseif {$ind > 9} {
       set name($i,0) P0$ind
    } else {
       set name($i,0) P00$ind
    }
    incr ind
  }
  # calculate the coordinates
  set ind 0
  set atomcntr 1
  for {set i 0} {$i<$m} {incr i} {
	for {set j 0} {$j<$n} {incr j} {
		if {($i+$j)%2 == 0} {
			set xpos [expr $i*$half_a]
            set ypos [expr $j*$half_a]
			set zpos 0.001
			
			if {$atomcntr == 1001} {
				set atomcntr 0
				incr segcounter
				puts $fd_pgn "segment $seg$segcounter {residue 1 PLA}"
			}
			if {$segcounter == 1} {
			puts $fd_top "ATOM $name($ind,0) PT      0.00"
			}			
			puts $fd_pgn "coord $seg$segcounter 1 $name($ind,0) {[expr ($xpos+$x_shift)] [expr ($ypos+$y_shift)] [expr ($zpos+$z_shift)]}"

			incr ind
			incr atomcntr
		}		
	}
  }
  
  puts $fd_pgn "writepsf $output.psf"
  puts $fd_pgn "writepdb $output.pdb"
  
  close $fd_top
  close $fd_pgn
  
  exec psfgen temp_platinum.pgn
  #exec rm temp_platinum.top
  #exec rm temp_platinum.pgn
}
