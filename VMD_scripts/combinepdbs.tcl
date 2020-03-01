package require psfgen 
readpsf platbot.psf 
readpsf plattop.psf 
readpsf solvate.psf 

coordpdb platbot.pdb 
coordpdb plattop.pdb 
coordpdb solvate.pdb 

writepsf all.psf 
writepdb all.pdb 