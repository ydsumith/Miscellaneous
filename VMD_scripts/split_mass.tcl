# copyright 2017 Sumith Yesudasan
# script to split the fibrinogen along x-axis
# then calculate the mass of each fragment

set N 10
set outputfile coarse

set all [atomselect top all]
set minmax [measure minmax $all] 
set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]] 
set L [lindex $vec 0]
set B [lindex $vec 1]
set H [lindex $vec 2]
set xmin [lindex [lindex $minmax 0] 0]
set xmax [lindex [lindex $minmax 1] 0]

set dx [expr $L/$N]

puts "L = $L" 
puts "B = $B" 
puts "H = $H" 
puts "xmin = $xmin, xmax = $xmax"
puts "\nsplitting at every $dx Ang"

set totmass 0
for {set i 1} {$i <= $N} {incr i} {
	set xstart [expr $xmin+($i-1)*$dx]
	set xend [expr $xmin+$i*$dx]
	set current [atomselect top "x > $xstart and x < $xend"]
	set mass 0
	foreach frag_mass [$current get mass] {
		set mass [expr {$mass + $frag_mass}]
	}
	puts "$i,  $xstart to $xend, mass = $mass"
	set totmass [expr {$totmass + $mass}]
}
puts "\nTotal Mass = $totmass"
