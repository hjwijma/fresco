#
# FarEnoughZone.mcr
# ===================
#
# Last updated by Hein Wijma, March 12th 2026
#
#
# Need to define on commandline the MacroTarget AvoidedResidue and the AvoidDistance
# E.g. yasara -txt "MacroTarget = '1NWW'" "AvoidResidue = 'HPN'" "AvoidDistance = 5"
# This will result in a table file with all residues more than 5 Angstrom from the residue HPN in 1NWW

cifTarget = FILESIZE (MacroTarget).cif
if cifTarget 
  LOADCIF (MacroTarget)
else 
  LOADPDB (MacroTarget)
TABULATE
LISTRES RES protein with distance > (AvoidDistance) from res (AvoidResidue), FORMAT=RESNAME1 MOLNAME RESNUM
SAVETAB 1, (MacroTarget)_MoreThan(AvoidDistance)AngstromFrom(AvoidResidue)_forDisulfides, columns=1,NumFormat=-1.2f
TABULATE 'END'
SAVETAB 1, (MacroTarget)_MoreThan(AvoidDistance)AngstromFrom(AvoidResidue), columns=1,NumFormat=-1.2f

exit
