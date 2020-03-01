##
## mol load pdb mypdb.pdb
## set sel [atomselect top protein]
##
## set matrix [parsematrix mypdb.pdb]  
##   get the move matrix according to the BIOMT REMARK in the pdb file
##   returns how many monomers are there, and return a list of 4x4 matrices
## mono2poly -o outname -chain {A B C D} $sel $matrix  
##  {A B C D} represent there are 4 monomers, if {} is used, chain is null
##  default of outname is poly , default of chain is {} 

proc parsematrix {orig_file}  {
  set infile [open  $orig_file r]
  set mtnum 0
  while {[gets $infile line]>=0} {
    set title [lindex $line 0]
    set biomt [lindex $line 2]
    if { $title == "REMARK" && [string match BIOMT? $biomt] } {
      set linenum [string index $biomt  end]
      if {$linenum == 1} {
        incr mtnum 1
        set matri($mtnum) {}
      }
      set lineelement [lrange $line 4 7]
      lappend matrix($mtnum) $lineelement
    }
  }
  close $infile

  set matrixlist {}
  foreach name [lsort -integer [array names matri]] {
    lappend matrix($name) {0.0 0.0 0.0 1.0}
    lappend matrixlist $matrix($name)
  }

  if {[llength matrixlist]==0} {
    error "There is no BIOTM REMARK information in this pdb file"
  }

  puts "this is a $mtnum polymer"
  return $matrixlist
}


proc mono2poly_Usage {} {
  puts "you should input:"
  puts "mono2ploy -chain chainlist -o outfilename sel matrix"
  error ""
}


proc mono2poly {args} {
  set cmdlinelength [llength $args]
  if {$cmdlinelength!=2 && $cmdlinelength<4} {
    mono2poly_Usage
  }
  set sel [lindex $args [expr $cmdlinelength-2]]
  set matrix [lindex $args end]
  set polynum [llength $matrix]
  puts $polynum
  set output poly
  set chainlist {}
  if {$cmdlinelength>=4} {
    set args [lrange $args 0 [expr $cmdlinelength-3]]
    set i 0
    while {$i<[llength $args]} {
      set opt [lindex $args $i]
      switch -exact -- $opt {
        -o { 
          incr i
          set output [lindex $args $i]
          incr i 
          continue
        }
        -chain { 
          
          incr i
          set chainlist [lindex $args $i]
          if { [llength $chainlist] != $polynum && [llength $chainlist]!=0} {
            error "chain list should have same element number with matrix or empty list"
          }
          incr i 
          continue
        }
        default {
          error "Unknown option $opt"
        }
      }
    }
  }
  $sel writepdb ${output}.org.pdb
  set out [open ${output}.pdb w]

  set i 0
  foreach mat $matrix {
    mol new ${output}.org.pdb waitfor all
    set mono [atomselect top all]
    $mono move $mat 
    if { [llength $chainlist] == $polynum } {
      $mono set chain [lindex $chainlist $i]
    } elseif { [llength $chainlist] == 0 } {
      set chainstr "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
      set chainchar [string index $chainstr $i]
      puts "Warning: making up sequential chain IDs, no chain list provided"
      puts "Using chain code '$chainchar'"
      $mono set chain $chainchar
    } else {
      puts "Warning: not setting monomer chains IDs, mismatched chain list size"
    }

    $mono writepdb ${output}.${i}.pdb
    set channel [open ${output}.${i}.pdb r]
    while {[gets $channel line]>=0} {
      if {[lindex $line 0]=="ATOM"} {
        puts $out $line
      }
    }
    close $channel
    file delete ${output}.${i}.pdb
    mol delete top
    incr i 
  }
  close $out
  file delete ${output}.org.pdb
  mol load pdb ${output}.pdb
  set all [atomselect top all]
  $all writepdb ${output}.pdb
  mol delete top
  return "${output}.pdb file have finished"
}

