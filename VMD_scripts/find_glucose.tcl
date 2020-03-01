mol new solvate.psf
mol addfile run1.dcd first 750 last -1 step 10 waitfor -1

set nf [molinfo top get numframes] 

set seltext1 "same residue as resname AGLC and within 2.1 of protein"
set seltext2 "same residue as resname AGLC and within 2.3 of protein"
set seltext3 "same residue as resname AGLC and within 2.5 of protein"
set seltext4 "same residue as resname AGLC and within 2.7 of protein"
set seltext5 "same residue as resname AGLC and within 2.9 of protein"

set outfile [open "data_set1.dat" a] 

puts $outfile "frame\twith2.1\twith2.3\twith2.5\twith2.7\twith2.9"

for {set i 0} {$i < $nf} {incr i} { 
puts "frame $i"
set selected1 [atomselect top $seltext1 frame $i] 
set selected2 [atomselect top $seltext2 frame $i] 
set selected3 [atomselect top $seltext3 frame $i] 
set selected4 [atomselect top $seltext4 frame $i] 
set selected5 [atomselect top $seltext5 frame $i] 

set number1 [expr "[$selected1 num]/24"]
set number2 [expr "[$selected2 num]/24"]
set number3 [expr "[$selected3 num]/24"]
set number4 [expr "[$selected4 num]/24"]
set number5 [expr "[$selected5 num]/24"]

puts $outfile "$i\t$number1\t$number2\t$number3\t$number4\t$number5"

}
close $outfile

mol delete all