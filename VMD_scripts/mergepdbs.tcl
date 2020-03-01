##
## Load and merge all of the PDB files in a given directory
## Usage: source mergepdbs.txt
## Output is written to "mergedpdb.psf" and "mergedpdb.pdb"
##

package require psfgen

# this is a hack, just to get the topology file
package require membrane
set topologyfile [format "%s/plugins/noarch/tcl/membrane1.0/top_all27_prot_lipid.inp" $env(VMDDIR)]
# topology /path/to/top_all27_prot_lipid.inp
topology $topologyfile
set nseg 1
foreach pdb [lsort [glob *.pdb]] {
  set segid V$nseg 
  segment $segid { 
    first NONE
    last NONE
    pdb $pdb 
  } 
  coordpdb $pdb $segid
  incr nseg
} 
guesscoord
writepsf mergedpdb.psf
writepdb mergedpdb.pdb

