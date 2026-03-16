#
#  Name Script: CombineMutations.mcr
#               ---------------------------
# 
# 
#  Purpose:     This scripts predicts the 3D-structure of a protein variant that has undergone one or multiple mutations. 
#
#  User modify: See instructions between two lines of pound (#) signs below for how to set the desired mutations and the name of the resulting PDB file.
#
#  Author:      Hein J. Wijma, University of Groningen, The Netherlands, H.J.Wijma@rug.nl
#
#  Requirement: Yasara-Structure 
#
#  Other info:  This script applies the desired side-chain mutations and allows the surrounding residues to adapt to the mutations. The script has been tested succesfully on 
#               several haloalkane dehalogenase mutant variants of which the structure had been solved by crystallography (3FWH, 3FBW,1HDE, 1BEE). A template PDB file of the 
#               unaltered protein should be provided (The MacroTarget). Before using this script, hydrogens should have been added to the PDB file, and the structure should 
#               have been inspected to see if the predicted protonation states are indeed reasonable and in agreement with the catalytic mechanism. If this is not the case 
#               for certain sidechains,the protonation state of those residues should be corrected manually to obtain maximal accuracy .  
# 
#  Usage:       Under any  operating system: the script can be run from the Yasara menu (Options -> Macro & Movie -> Play Macro). Then the template PDB file should be set as 
#               the MacroTarget in the usual way while one also needs to set the table file via the command line below with for example: "TableFile = 'Hhec2360MutList'"
#
#        


#  TableFile:   The Table file should have first a single line with comments that will not be read by yasara follwed by Original amino acid [SPACES]  position number  [SPACES]  new amino acid
#  =========----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
#Example Content TableFile (This would convert Ile 150 to a Gln, etc.  
#  I  150   Q
#  T  230   V
#  A  287   F



# This removes anything left in memory that can perturb the current calculations.
CLEAR

##########################################################################################################################################################################################
################################# Adapt this part by hand if needed ######################################################################################################################

# By default, a name of the mutant protein is provided based on the name of the TableFile. 
# Alternatively, it is possible to have the script provide a systematic name automatical (of the type: _D271A_C398S). Then set the following option to 'YES' ('Yes' does not work)

AutomaticNameExtension = 'YES'


# Decide if you want to save a Yasara scene file with a superposition of the predicted structure on the template ('YES').

KeepSuperPosition = 'YES'



##########################################################################################################################################################################################
####################### The part below does not need to be changed #######################################################################################################################

# Load the table and convert it to a list of residues. 
LOADTAB (TableFile)
ListEverythingInTable() = TAB 1
DELTAB 1

# convert the ListEverythingInTable into 3 letter codes
for i = 1 to count ListEverythingInTable
  if (ListEverythingInTable(i) == 'A')
    ListEverythingInTable(i) = 'ALA' 
  if (ListEverythingInTable(i) == 'R')
    ListEverythingInTable(i) = 'ARG' 
  if (ListEverythingInTable(i) == 'N')
    ListEverythingInTable(i) = 'ASN' 
  if (ListEverythingInTable(i) == 'D')
    ListEverythingInTable(i) = 'ASP' 
  if (ListEverythingInTable(i) == 'C')
    ListEverythingInTable(i) = 'CYS' 
  if (ListEverythingInTable(i) == 'Q')
    ListEverythingInTable(i) = 'GLN' 
  if (ListEverythingInTable(i) == 'E')
    ListEverythingInTable(i) = 'GLU' 
  if (ListEverythingInTable(i) == 'G')
    ListEverythingInTable(i) = 'GLY' 
  if (ListEverythingInTable(i) == 'H')
    ListEverythingInTable(i) = 'HIS' 
  if (ListEverythingInTable(i) == 'I')
    ListEverythingInTable(i) = 'ILE' 
  if (ListEverythingInTable(i) == 'L')
    ListEverythingInTable(i) = 'LEU' 
  if (ListEverythingInTable(i) == 'K')
    ListEverythingInTable(i) = 'LYS' 
  if (ListEverythingInTable(i) == 'M')
    ListEverythingInTable(i) = 'MET' 
  if (ListEverythingInTable(i) == 'F')
    ListEverythingInTable(i) = 'PHE' 
  if (ListEverythingInTable(i) == 'P')
    ListEverythingInTable(i) = 'PRO' 
  if (ListEverythingInTable(i) == 'S')
    ListEverythingInTable(i) = 'SER' 
  if (ListEverythingInTable(i) == 'T')
    ListEverythingInTable(i) = 'THR' 
  if (ListEverythingInTable(i) == 'W')
    ListEverythingInTable(i) = 'TRP' 
  if (ListEverythingInTable(i) == 'Y')
    ListEverythingInTable(i) = 'TYR' 
  if (ListEverythingInTable(i) == 'V')
    ListEverythingInTable(i) = 'VAL' 
 


j = 1
for i = 1 to count ListEveryThingInTable step 3
  ListOriginalResidues(j) = '(ListEverythingInTable((i)))'
  ListResidueNumbers(j)   = '(0+(ListEverythingInTable((i)+ 1)))'
  ListReplacingResidues(j)= '(ListEverythingInTable((i) + 2))'
  # For Debugging
  # print (ListOriginalResidues(j)) (0+(ListResidueNumbers(j))) (ListReplacingResidues(j))
  j = j+ 1
  
  
# Convert to 3 letter codes for the residues


# Provide a name for the current variant.The resulting PDB file can be recognized by this extension of the name of the template. 
CurrentName = '(TableFile)'


# This loads the template PDB file
PdbFilePresent = FILESIZE (MacroTarget).pdb
if (PdbFilePresent)
  LOADPDB (MacroTarget)
else 
  LOADYOB (MacroTarget) 

# This sets parameters for forcefield and energy calculations. 
ForceField YASARA2, SETPAR = YES
SURFPAR Resolution = 3 


# This does the desired mutagenesis and provides a clear presentation of what is happening to the protein. 
STYLE TUBE
for i = 1 to count ListOriginalResidues
  # verify that the original residue is correct, otherwise exit
  CorrectIsOne = COUNTRES (ListOriginalResidues(i)) (0+(ListResidueNumbers(i)))
  if (CorrectIsOne)  
    # This does the mutagenesis
    SWAPRES (ListOriginalResidues(i)) (0+(ListResidueNumbers(i))), (ListReplacingResidues(i))
    # This is to prevent clashes with existing water molecules in the template file, initial testing suggested 2.3 Angstrom was a reasonable distance.
    DELRES HOH with distance < 2.3 from (ListReplacingResidues(i)) (0+(ListResidueNumbers(i)))
    # This makes the altered residues stand out.
    COLORATOM res (ListReplacingResidues(i)) (0+(ListResidueNumbers(i))) element C, red
    SHOWRES  (ListReplacingResidues(i)) (0+(ListResidueNumbers(i)))
    STICKRES (ListReplacingResidues(i)) (0+(ListResidueNumbers(i)))
  else
    RAISEERROR Wrong original residue in list at (ListOriginalResidues(i)) (0+(ListResidueNumbers(i)))
    exit

# This makes a string of the positions that are mutated, needed later in this script
MutatedResiduesNumbers =''
for i = 1 to count ListResidueNumbers
  MutatedResiduesNumbers = '(MutatedResiduesNumbers) (0+(ListResidueNumbers(i)))'


# Uncomment below to show show the residues nearby the mutated residues
#SHOWRES protein with distance < 6 from protein res (MutatedResiduesNumbers)

# uncomment this to let the user inspect the structure and positions of the mutations before minimization
WAIT ContinueButton

# These 6 cycles optimize the new mutations with stepwise DEE optimization based on rotamers, followed by a local energy minimization in water. 
# The volume that is energy optimized starts at 7 Angstrom from the mutated residues and increases with 1 Angstrom at every cycle. 
# At the end there is an energy minimization of the entire protein. 
for i = 1 to 6
  NICEORIOBJ 1
  SHOWMESSAGE Round (i) in optimizing
  CELL Auto, Extension = 10,shape=cuboid
  OPTIMIZERES protein res (MutatedResiduesNumbers), method = SCWALL
  if (i == 6)
    NICEORIOBJ 1
    CELL Auto, Extension=7.5 ,shape=cuboid, protein
    BOUNDARY Periodic
  else  
    CELL Auto, Extension=(6.0 + (i)) ,shape=cuboid,  protein res (MutatedResiduesNumbers)
    BOUNDARY Wall
    LONGRANGE None
  FillCellWater Density=0.997,Probe=1.4,BumpSum=1.0,DisMax=0
  HIDERES HOH
  Experiment Minimization
  Experiment On
  STICKRES HOH
  COLORRES OBJ 1 res HOH, magenta
  Wait ExpEnd
  DELRES OBJ 3
  ForceField YASARA2, SETPAR = YES


# The following procedure creates the name of the mutations in the format _D16C_C150S_A201C
# This name will only be used if above the automatic naming is switched on (AutomaticNameExtension = 'YES')
MutationsList =''
for i = 1 to count ListOriginalResidues
  # Add an underscore '_'
  MutationsList = '(MutationsList)_'
  # Add the name of the original residue
  Mut(i) = '(ListOriginalResidues(i))'
  if (Mut(i) == 'ALA')
    Mut(i) = 'A' 
  if (Mut(i) == 'ARG')
    Mut(i) = 'R' 
  if (Mut(i) == 'ASN')
    Mut(i) = 'N' 
  if (Mut(i) == 'ASP')
    Mut(i) = 'D' 
  if (Mut(i) == 'CYS')
    Mut(i) = 'C' 
  if (Mut(i) == 'GLN')
    Mut(i) = 'Q' 
  if (Mut(i) == 'GLU')
    Mut(i) = 'E' 
  if (Mut(i) == 'GLY')
    Mut(i) = 'G' 
  if (Mut(i) == 'HIS')
    Mut(i) = 'H' 
  if (Mut(i) == 'ILE')
    Mut(i) = 'I' 
  if (Mut(i) == 'LEU')
    Mut(i) = 'L' 
  if (Mut(i) == 'LYS')
    Mut(i) = 'K' 
  if (Mut(i) == 'MET')
    Mut(i) = 'M' 
  if (Mut(i) == 'PHE')
    Mut(i) = 'F' 
  if (Mut(i) == 'PRO')
    Mut(i) = 'P' 
  if (Mut(i) == 'SER')
    Mut(i) = 'S' 
  if (Mut(i) == 'THR')
    Mut(i) = 'T' 
  if (Mut(i) == 'TRP')
    Mut(i) = 'W' 
  if (Mut(i) == 'TYR')
    Mut(i) = 'Y' 
  if (Mut(i) == 'VAL')
    Mut(i) = 'V' 
  MutationsList = '(MutationsList)(Mut(i))'
  # Add the number
  MutationsList = '(MutationsList)(0+ (ListResidueNumbers(i)))'
  # Add the name of the new residue
  Mut(i) = '(ListReplacingResidues(i))'
  if (Mut(i) == 'ALA')
    Mut(i) = 'A' 
  if (Mut(i) == 'ARG')
    Mut(i) = 'R' 
  if (Mut(i) == 'ASN')
    Mut(i) = 'N' 
  if (Mut(i) == 'ASP')
    Mut(i) = 'D' 
  if (Mut(i) == 'CYS')
    Mut(i) = 'C' 
  if (Mut(i) == 'GLN')
    Mut(i) = 'Q' 
  if (Mut(i) == 'GLU')
    Mut(i) = 'E' 
  if (Mut(i) == 'GLY')
    Mut(i) = 'G' 
  if (Mut(i) == 'HIS')
    Mut(i) = 'H' 
  if (Mut(i) == 'ILE')
    Mut(i) = 'I' 
  if (Mut(i) == 'LEU')
    Mut(i) = 'L' 
  if (Mut(i) == 'LYS')
    Mut(i) = 'K' 
  if (Mut(i) == 'MET')
    Mut(i) = 'M' 
  if (Mut(i) == 'PHE')
    Mut(i) = 'F' 
  if (Mut(i) == 'PRO')
    Mut(i) = 'P' 
  if (Mut(i) == 'SER')
    Mut(i) = 'S' 
  if (Mut(i) == 'THR')
    Mut(i) = 'T' 
  if (Mut(i) == 'TRP')
    Mut(i) = 'W' 
  if (Mut(i) == 'TYR')
    Mut(i) = 'Y' 
  if (Mut(i) == 'VAL')
    Mut(i) = 'V' 
  MutationsList = '(MutationsList)(Mut(i))'  



# Save with the new name
if (AutomaticNameExtension == 'YES')
  Tag = '(MutationsList)'
else
  Tag = '(CurrentName)'
if (PdbFilePresent)
  SAVEPDB OBJ 1, '(MacroTarget)_(Tag)'
else
  SAVEYOB OBJ 1, '(MacroTarget)_(Tag)'


# now compare visually with the original    
PdbFilePresent = FILESIZE (MacroTarget).pdb
if (PdbFilePresent)
  LOADPDB (MacroTarget)
else 
  LOADYOB (MacroTarget) 
SUPATOM OBJ 3 atom CA, OBJ 1 atom CA, match = yes


# Further visualization
showres all with distance < 1 from visible 
HIDEATOM element H with bond to atom element C


# Save the superposition if desired
if (KeepSuperPosition == 'YES')
  # Also save the superposition
  SAVESCE '(MacroTarget)_(Tag)_SuperPositionOriginal'
  
SHOWMESSAGE Finished with modelling, shown is superposition with original
RENAMEOBJ 1, Model
RENAMEOBJ 3, Original
