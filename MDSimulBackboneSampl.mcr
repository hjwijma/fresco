#
#  Name Script: 	MDASimulBackboneSampl.mcr
#               	------------------
#
#  Author:      	Hein J. Wijma, University of Groningen, The Netherlands, Biochemical Laboratory, Biotransformation and Biocatalysis group
#                       H.J.Wijma@rug.nl
#                       Parts of this script were inspired by Elmar Krieger his macros. 
# 
#  Purpose:             - Can be used to create a series of differently initialized trajectories of which the relevant angles and distance are sampled on the fly, during 
#                         the MD simulation.
#                       - script written such that it can be easily adapted for a new enzyme
#
#  Requirements: 	- YASARA-Structure 
#                	- yob file with correct protonation states of the entire protein and of the ligand, this is the MacroTarget. 
#                          
#  Results:     	- MD trajectory(ies) for the (differently seeded) simulations(s) that can be inspected with Elmar Krieger's md_play macro
#               	- Several Table Files with Near Attack Conformation percentages and other statistical data, such as for which NAC criteria the requirements are met
#                       - overview flexibility of the protein versus flexibility of the X-ray structure (_RMSF.tab)
#                       - PDB file with average structure from the MD simulation
#                       - table files with energy and RMSD from starting structure versus simulation time
#  
#  User modify: 	- if no definitions are available for your enzyme, you will have to insert suitable definition for the Near Attack Conformation (NAC) in your enzyme.
#                       - suitable definitions consist of (see below for more details)
#                            	- name of the criterion (12 characters, e.g. '  AngletoNH1', spaces are only possible at the start)
#                               - how to measure the criterion
#                               - minimal requirement, should always be defined even if the smaller the better (e.g. define then -0.01 Angstrom for the minimal distance)
#                               - maximally allowed, should also always be defined, even if the larger the better (e.g. define an angle of 1000 degrees as maximum).
# 
#  Important:   	- The produced average PDB files standard have RMSF instead of b-factors. This is done since with RMSF, with its unit in Angstrom, the flexibility is 
#                         more easy to grasp than with b-factors. This can be changed back to b-factors if desired as described below under UseRMSFInPdbYOBFile
# 	
#  Further:     	Under linux this script is easiest printed with the following or an equivalent command: mpage -l2 -W180 -H -m50 -r DMDAnalysisNacs.mcr | lpr
#

# ========================================================================================================================================================================================
# =========================  Part A of the Script: The settings and definitions ==========================================================================================================
# ========================================================================================================================================================================================

# ==== WarmUpRegime ====
# How long to take to warm up from 5 K to 298 K (or the assay temperature )
# 'F' =  3000 fs
# 'N' = 30000 fs
# Given here two options that determine also how the name of the snapshots are called, this to prevent and MD simulation to restart with the wrong snapshots(e.g., otherwise a 
# 30 ns warmup trajectory run might accidentely continue with the snapshots of a 3 ns warmup trajectory. Options are 'F' for Fast and 'N' for Normal

# ==== EquilibrationTime ====
# How long to equilibrate after warming up? The unit here is fs. Please NOTICE that the WarmUpRegime and Equilibration time should be dividable by the SnapshotInterval below. 

# ==== ProductionTime ====
# How long should the production run then take? The unit here is fs.

# ==== MS ====
# How many individual trajectories to run. These trajectories will start with different initial velocities 

# ==== IntervalAnalysis ====
# ===== TrajectorySubAnalysisIntervals() =5000,50000
# TRUE 

# MoreLOGFile
# TRUE, or anything else

# ==== UseRMSFInPdbYOBFile ====
# anything else than 'TRUE' will result in a pdb file with the bfactors from the MD simulation


# ==== SnapShotInterval ====      
# How often to take snapshot to record (most desirable is once every 2500 fs but once every 5000 fs is OK for these calculations to prevent cluttering the harddisk). 
# THIS VALUE SHOULD BE DIVIDABLE BY THE NACSamplingInterval, SEE BELOW

# ==== NACSamplingInterval ====
# How often to sample NACs on the fly? (desirable is every 20 or 100 femtosecond (unit is not picosecond). This value needs to be dividable through 2.5 and 4 to be compatible with the 
# possible timesteps (2.5 for normal, 4.0 for LINCS/SETTLE algorithms switched on).  

# ==== NACShotEveryXSnapshots ====
# Per how many snapshots to record a tab file with the NAC information stored. It is desirable to do this less often than snapshots (less files) but often enough not to waste time 
# during a restart. BEST TO KEEP IT AT 1 unless the run will produced many snapshots.DO NOT CHANGE THIS NUMBER HALF WAY DURING A SIMULATION.

# ==== AssayTemperature ====
# Temperature during the run (standard 298 K, higher temperatures would need a shorter timestep in simulations which needs to be modified at line 270).  

# ==== EmployedForcefield ====
# which forcefield to use, will not be visible in final file names

# ==== LincsSettle ====
# Use LINCS and SETTLE algorithms to make the simulation faster, constrains angles of water molecules and bond lengths that involve hydrogen atoms
# Has to be either 'On' or 'Off', Nothing else. Visible in all files as LSOff and LSOn. If on the simulation will be much faster but possibly less accurate. 

# ==== CalculateAverageTrajectory ====
# calculate the average trajectory for RMSD (= 'On' or anything else to set it to Off).



##########################################################################################################################################################################################
##################### The four settings below have to altered by the users regulary ######################################################################################################

# Define whether to automatically exit at end, YES if you are running macros under -txt mode
# Set here whether we want Yasara to exit after finishing the script by choosing for 'YES' or 'NO'
AutomaticExitAtEnd = 'YES'

# Set the number of processors that Yasara is allowed to use. 
PROCESSORS 4

# Define Current Target of the script, the script needs to know which enzyme your are working with
CT='NoNACMeasurements'


# Set the mode (see below for what these settings of ForDebuggingDefinitions, MultiShort (good sampling of conformations at low cpu cost), and SingleLong (allow long MD simulation
# to allow the protein to find alternative conformations. You can also define your own settings
CurrentSetting = 'MDforDisulfideBonds'


##########################################################################################################################################################################################
####################### The part below does not need to be adapted unless a new enzyme (CT) is adopted or unless a new MD protocol is desired ############################################


# Parameters 1: Simulation time periods for the phases of warming, equilibration, and production 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if (CurrentSetting == 'ForDebuggingDefinitions')
  WarmUpRegime           = 'F'
  EquilibrationTime      = 1000
  ProductionTime         = 2000
  MS                     = 2
  InterValAnalysis       = 'Off'
  UseRMSFInPdbYOBFile    = 'TRUE'
  MoreLOGFile            = 'TRUE'
  SnapShotInterval       = 1000
  NACSamplingInterval    = 20
  NACShotEveryXSnapshots = 1
  AssayTemperature       = 298
  EmployedForcefield     =  'Yamber3'
  LincsSettle            = 'Off'
  CalculateAverageTrajectory = 'Off'
  SkipDifferentMinimization = 'TRUE'


if (CurrentSetting == 'MDforDisulfideBonds')
  WarmUpRegime           = 'N'
  EquilibrationTime      = 470000
  ProductionTime         = 2000000
  MS                     = 1
  InterValAnalysis       = 'Off'   
  TrajectorySubAnalysisIntervals() =5000,50000
  UseRMSFInPdbYOBFile    = 'TRUE'
  MoreLOGFile            = 'UnTRUE'
  SnapShotInterval       = 25000
  NACSamplingInterval    = 100
  NACShotEveryXSnapshots = 1
  AssayTemperature       = 298
  EmployedForcefield     =  'Yamber3' 
  LincsSettle            = 'On'
  CalculateAverageTrajectory = 'Off'
  SkipDifferentMinimization = 'TRUE'
  MakeSnapshotPdbs  = 'Yes'
  SnapShotPdbsTimeps= 500,750, 1000,1250,1500,1750,2000,2250,2500 

if (CurrentSetting == 'OwnSettings')
  # here the user needs to define his own settings as above, see instructions above
  # especially be sure to have the ratios for NAC sampling interval with the snapshots right


# the settings for how the active site looks like
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
if CT=='NoNACMeasurements'
  pH = 7
  ActiveSiteResidues = 'Sub'
  StayAwayDistance = 5
  print No NAC definitions, therefore in all files 100 % NAC will be reported
  


# This verifies that some needed files are present and no setting are wrong.Some other files and settings are checked elsewhere
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 

if MacroTarget == ''
  RaiseError The Enzyme Substrate complex Target is missing
if CT == ''
  RaiseError The Current Target, which enzyme whe are working with, is missing
if WarmUpRegime != 'F'
  if WarmUpRegime != 'N'
    RaiseError Illegal Warm up regime defined in script, (WarmUpRegime) should be either 'N' or 'F' for normal or fast warming


# ========================================================================================================================================================================================
# =========================  Part B of the Script: Neutralization with salt and energy minimisation ======================================================================================
# ========================================================================================================================================================================================
# This removes anything left in memory that can hinder the current calculations.
CLEAR
#  Assumes the protonation states of the His/Glu/Asp residues have been set manually. 
#  What the procedure does is that it uses the existing Yasara Neutralization procedure to place salt ions. Then it reverts to the original structure
#  (with the manually set ionization states), and places waters with the command FillCellWater (does not alter ionization states) and puts all salt ions that are 
#  acceptable inside. After that the cell is unlikely to be neutral, and it positions additional salt ions to neutralize the cell again. 
#

if MoreLOGFile == 'TRUE'
  CONSOLE On
else
  CONSOLE Off


UnMinimizedStructureExists = FILESIZE (MacroTarget)_SaltWaterBoxNotMinimized.sce
if !(UnMinimizedStructureExists)
  FileThere = FILESIZE (MacroTarget).yob
  if (FileThere) 
    LOADYOB (MacroTarget)
  else 
    FileThere = FILESIZE (MacroTarget).pdb
    if (FileThere)
      LOADPDB (MacroTarget)
    else
      RAISEERROR No MacroTarget Found
  # some settings
  NICEORIOBJ 1
  FORCEFIELD (EmployedForcefield), SETPAR = YES

  # 2018-07-23 replaced because something changed in Yasara  
  #CELL AUTO, extension = 7.5, Scale = Yes
  CELL Auto, Extension=7.5,Shape=Cuboid

  # make a copy of the protein and remove it temporarily
  DUPLICATEOBJ 1
  RENUMBEROBJ 3,4
  REMOVEOBJ 4
  
  # Get salt ions vio the normal procedure
  FORCEFIELD (EmployedForcefield)
  EXPERIMENT Neutralization
    WaterDensity 0.997
    pH (pH)
    NaCl 0.5
    pKaFile (MacroTarget).pka
    Speed Fast
  FIXATOM water 
  EXPERIMENT On
  WAIT ExpEnd
  
  # now keep the salt ions, but delete the protein since the ionizations states have altered now. The macro requires the user to provide a pdb file with the correct ionization states.
  DELOBJ 1
  DELRES element O
  ADDOBJ 4
  RENUMBEROBJ 4,1
  RENUMBEROBJ 3,4
  REMOVEOBJ 4
  
  # use command FillCellWater since with original neutralization too many waters placed inside the protein. Then a protocol based on that used by Elmar Krieger
  # -----------------------------------------------------------------------------------------------------------------------------------------------------------  
  LONGRANGE none
  FILLCELLWATER Density=0.997,Probe=1.4,BumpSum=1.0,DisMax=0
  # Fix all the heavy atoms of object 1, the original
  FIXATOM OBJ 1
  FREERES HOH
  Style Stick
  
  # some settings of pressure and fast electrostatics
  PRESSURECTRL Off 
  SIMSPEED Fast
  TEMPERATURE (AssayTemperature)
  # switch electrostatics off to prevent local minima
  INTERACTIONS Bond,Angle,Dihedral,Planarity,VdW
  TIMESTEP 2,2.00
  # now 25 timesteps of steepes descent minimization of OBj 3 water and hydrogen atoms
  TempCtrl SteepDes
  FIXBOND water, water
  SIM ON
  WAIT 25
  # switch electrostatics back on and do a short simulated annealing, too long will create vacuum bubbles
  INTERACTIONS Bond,Angle,Dihedral,Planarity,Coulomb,VdW
  TEMPCTRL Anneal
  WAIT 100
  TEMPERATURE (AssayTemperature)  
  # Followed by 100 steps normal dynamics to ensure the water is OK.   
  Wait 100
  FREERES ALL
  Timestep 2,1.25
  SIM OFF
  
  
  # Now add the previously added salt ions and remove waters that are too close to the salt or salt that is too close to the active site. 
  ADDOBJ 4
  RENUMBEROBJ 3,5
  JOINOBJ 5,4
  RENUMBEROBJ 4,3
  DELRES OBJ 3 res element Na Cl with distance < (StayAwayDistance) from res (ActiveSiteResidues)
  DELRES OBJ 3 res element O with distance < 1.75 from res OBJ 3 atom element Cl
  DELRES OBJ 3 res element O with distance < 1.02 from res OBJ 3 atom element Na
  
  
  # some settings of pressure and fast electrostatics
  PRESSURECTRL Type=Combined,Pressure=1.000,Name=HOH,Density=0.997,Axis=XYZ
  TIMESTEP 2,1.25
  LONGRANGE Coulomb
  BOUNDARY periodic
  ENERGYUNIT kJ/mol
  
  # now need an algorithms to check every water from OBj 3 that lives more than 4 angrstrom from the protein, 
  # find the one at the most positive/negative place place, change it, relist all waters, find again the most pos/negative, untill neutral. 
  # present them by Ballres 
  
  # Check what the charge of the system is before neutralization
  CurrentChargeList() = CHARGEOBJ all 
  CurrentCharge = sum CurrentChargeList
  CurrentCharge = 0 + (CurrentCharge)
  
  # determine what kind of ions, and how many, to introduce for neutralization
  ReplaceWithSodiumIons = 0
  ReplaceWithChlorideIons = 0
  if (CurrentCharge < 0.5)
    ReplaceWithSodiumIons = - (CurrentCharge)
  if (CurrentCharge > 0.5)  
    ReplaceWithChlorideIons = (CurrentCharge)
  
  # now do the replacements, start with the ion for neutralization only after those are finished add the ones for salt.
  If ((ReplaceWithSodiumIons) > (ReplaceWithChlorideIons))
    FirstIon = 'SOD'
    SecondIon = 'CHL'
    FirstChargeCount  = 0 + (ReplaceWithSodiumIons)
    SecondChargeCount =  0 + (ReplaceWithChlorideIons)
  else
    FirstIon = 'CHL'
    SecondIon = 'SOD'  
    FirstChargeCount  = 0 + (ReplaceWithChlorideIons)
    SecondChargeCount =  0 + (ReplaceWithSodiumIons)
  
  if (FirstChargeCount)
    HIDERES HOH
    REMOVEENVRES HOH
    for CurrentStage in 'FirstRound', 'SecondRound'
    if CurrentStage == 'FirstRound'
      CurrentIon = 'FirstIon'
      ChargeCount = (FirstChargeCount)
    if CurrentStage == 'SecondRound'
      CurrentIon = 'SecondIon'
      ChargeCount = (SecondChargeCount)
    for i = 1 to ChargeCount
      # reset the lists 
      EligibleH2OList() = 0
      ESPSURFLIST()     = 0  
      # get the list of waters that are part of Object 3 and more than 6 angstrom away from the , the list consists of atom numbers
      SIM on
      SIM pause
      EligibleH2OList() = LISTRES OBJ 3 res HOH with distance > 6 from res SOD CHL CIM CIP
      SIM Off
      # determine the Electrostatic potential at their surface (the list consists of surface and potential energy assuming the atom charge = +1)
      ESPSURFLIST() = SURFESPRES OBJ 3 res HOH with distance > 6 from res SOD CHL, type = accessible, method = PME, unit = RES
      # seems bug in yasara in which the first value for ESPATOM is always nan
      ESPSURFLIST(2) = 0
      for j = 1 to count EligibleH2OList
        CurrentCoulomb = ESPSURFLIST(j*2)
        if j==1
          MostNegativeSpot  = (CurrentCoulomb)
          MostPositiveSpot  = (CurrentCoulomb)
          IndexMostNegative = 1
          IndexMostPositive = 1
        if CurrentCoulomb<MostNegativeSpot
          CriticalDistance = GroupDistance res (ActiveSiteResidues),  atom (EligibleH2OList(j))
          if (CriticalDistance > StayAwayDistance)
            MostNegativeSpot  = (CurrentCoulomb)
            IndexMostNegative = (j)
        if CurrentCoulomb>MostPositiveSpot
          CriticalDistance = GroupDistance res (ActiveSiteResidues),  atom (EligibleH2OList(j))
          if (CriticalDistance > StayAwayDistance)
            MostPositiveSpot  = (CurrentCoulomb)
            IndexMostPositive = (j)
      # do the replacements, do not forget to delete the hydrogen atoms that belong to it (be carefull the atom numbers might shift)
      if (ReplaceWithSodiumIons) > 0.5
        if ((CurrentIon) == 'SOD')
          RENAMERES atom (EligibleH2OList(IndexMostNegative)), SOD
          ReplaceWithSodiumIons = (ReplaceWithSodiumIons) - 1
      if (ReplaceWithChlorideIons) > 0.5
        if ((CurrentIon) == 'CHL')
          RENAMERES atom (EligibleH2OList(IndexMostPositive)), CHL
          ReplaceWithChlorideIons = (ReplaceWithChlorideIons) - 1
      DELATOM RES CHL SOD element H
      SWAPATOM RES SOD, Na
      SWAPATOM RES CHL, Cl
      BALLRES RES SOD CHL
      ADDENVRES SOD CHL
      SIM init

  SIM OFF
  DELRES OBJ 3 Res HOH with distance < 2.8 from Res Sub
  
  SAVESCE (MacroTarget)_SaltWaterBoxNotMinimized



# also check the conformations in the Unminimized structure, if that has not yet been done
# Magic number, i.e. NACs, and all criteria with their NAC
MagicNumber = 1 + (2*(count Label_))

LOADSCE (MacroTarget)_SaltWaterBoxNotMinimized
FORCEFIELD (EmployedForcefield), SETPAR = YES
SIM On
for c = 1 to count Label_
  (Label_(c))_Q_       = ((Label_(c))_Measurement) 
# determine whether NAC is achieved 
SIM Off
NAC = 1
for c = 1 to count Label_
  (Label_(c))_NAC = 0.00000
  if (  ((Label_(c))_Q_)  > ((Label_(c))_Q_Min))
    if (((Label_(c))_Q_)  < ((Label_(c))_Q_Max))
      (Label_(c))_NAC = 1.00000
  NAC = (NAC) * ((Label_(c))_NAC)  
TABULATE (100.000 * (NAC))
for c = 1 to count Label_
  TABULATE ((Label_(c))_Q_)
  TABULATE ((Label_(c))_NAC)
  
LabelCollection = ' ________NAC?'
for c = 1 to count Label_
  LabelCollection = '(LabelCollection) (Label_(c)) CriteriaMet?'
SaveTab 1,(MacroTarget)_NAC_Results_UnMinimized,Format=Text,Columns=(MagicNumber),NumFormat=%12.2f,(LabelCollection)
DELTAB 1
  
# ========================================================================================================================================================================================
# =========================  Part C of the Script: The MD simulation(s) ==================================================================================================================
# ========================================================================================================================================================================================


#  **********************************  HERE STARTS THE LOOP THAT DOES THE ENTIRE MD SIMULATION   *****************************************************************************************
#  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 

# This removes anything left in memory that can hinder the current calculations.
CLEAR

MaximumSeeds = 00 + (MS)
for CSN = 01 to (MaximumSeeds)
  # for multi trajectories
  SeedingNumber = (0.001234567*(CSN))
  
  # The warmuptime is 
  #  do not go under 3000 fs since otherwise it goes so fast that the simulation is out of control by the time room temperature is reached
  if WarmUpRegime == 'F'
    WarmUpTime = 3000
  if WarmUpRegime == 'N'
    WarmUpTime = 30000
  

  # Settings of the forcefield and temperature pressure control, saving snapshots, etcetera
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  
  # Use the above defined ForceField
  FORCEFIELD (EmployedForcefield ),SETPAR = Yes
  
  # Ensure no fixed atoms present
  FREE

  
  # use periodic boundary conditions
  BOUNDARY Periodic

  # Use the earlier made saltwater box not minimized structure, and the minimized _water.sce if it exists, else do an independently seeded minimization
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  
  # ensure some files exists, to prevent problems later
  UnMinimizedStructureExists = FILESIZE (MacroTarget)_SaltWaterBoxNotMinimized.sce
  if (UnMinimizedStructureExists)
    MinimizedStructureExists = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)_water.sce
    if (MinimizedStructureExists)
      LOADSCE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)_water.sce
    else
      if (SkipDifferentMinimization != 'TRUE')
        LOADSCE (MacroTarget)_SaltWaterBoxNotMinimized
        # now a trick to get differently initilized minimizations, by putting a water in front
        JOINRES OBJ 3
        for i = 1 to (CSN) 
          AllWaters()    = LISTRES OBJ 3 res HOH element O
          NumberOfWaters = COUNTRES OBJ 3 res HOH element O
          SPLITRES OBJ 3 res (AllWaters(NumberOfWaters)) 
          SPLITOBJ 3
          JOINOBJ 1,4
          RENUMBEROBJ 4,1
        EXPERIMENT Minimization
        EXPERIMENT On
        WAIT ExpEnd
        # restore the pre-existing state
        JOINRES HOH
        SPLITOBJ 1
        TotalObjects = COUNTOBJ All
        JOINOBJ 1,3
        for i = 4 to ((CSN)+2) 
        RENUMBEROBJ 4,1
        for i = 5 to TotalObjects
          JOINOBJ (i), 1
        RENAMEOBJ 1, protein
        RENAMEOBJ 3, solution
        SAVESCE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)water
      else
        MinimizedStructureExists = FILESIZE (MacroTarget)_(WarmUpRegime)_01_LS(LincsSettle)water.sce 
        if (MinimizedStructureExists)
          shell cp (MacroTarget)_(WarmUpRegime)_01_LS(LincsSettle)water.sce (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)water.sce
        else 
          LOADSCE (MacroTarget)_SaltWaterBoxNotMinimized
          EXPERIMENT Minimization
          EXPERIMENT On
          WAIT ExpEnd
          SAVESCE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)_water
     
  else
    RAISEEROR Could not find (MacroTarget)_NotMinimized.sce, needed later in procedure
  
  
  
  # determine whether NAC is achieved, and tabulate it, after that reset it again 
  LOADSCE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)_water.sce
  NAC_Seed(CSN)_M = 1
  FORCEFIELD (EmployedForcefield), SETPAR = YES
  SIM On
  for c = 1 to count Label_
    (Label_(c))_Q_Seed(CSN)_M = ((Label_(c))_Measurement) 
  SIM Off
  for c = 1 to count Label_
    (Label_(c))_NAC_Seed(CSN)_M = 0.00000
    if (  ((Label_(c))_Q_Seed(CSN)_M)  > ((Label_(c))_Q_Min))
      if (((Label_(c))_Q_Seed(CSN)_M)  < ((Label_(c))_Q_Max))
        (Label_(c))_NAC_Seed(CSN)_M = 1.00000
    NAC_Seed(CSN)_M = ((NAC_Seed(CSN)_M) * ((Label_(c))_NAC_Seed(CSN)_M  ) )
  TABULATE (100.000*(NAC_Seed(CSN)_M))
  for c = 1 to count Label_
    TABULATE ((Label_(c))_Q_Seed(CSN)_M)
    TABULATE ((Label_(c))_NAC_Seed(CSN)_M)
  LabelCollection = ' ________NAC?'
  for c = 1 to count Label_
    LabelCollection = '(LabelCollection) (Label_(c)) CriteriaMet?'
  SAVETAB 1,(MacroTarget)_NAC_Results_Seed(CSN)_Minimized,Format=Text,Columns=(MagicNumber),NumFormat=%12.2f,(LabelCollection)
  DELTAB 1
  
  # use a cutoff for the non-bonded of 7.86 Angstrom
  CUTOFF 7.86
  
  # calculate further away than 7.86 Angstrom with Particle Mesh Ewald (PME) algorithm
  LONGRANGE coulomb
  
  # use the 4th degree splines for PME
  SIMSPEED normal
  
  # Correct for the diffusion of the protein
  CORRECTDRIFT On
  
  # settings for timesteps 
  if LincsSettle == 'Off'
    # use a timestep of 1.25 fs, update the non-bonded interactions every 2 timesteps
    CalculationTimeStep = 1.250
    UpdateCycle         = 2
    TIMESTEP (UpdateCycle),(CalculationTimeStep)
    # remove LINCS and SETTLE constraints
    FREEBOND All, All
    FREEANGLE All, All, All
  else
    if  LincsSettle == 'On' 
      # NEW FIXHYDANGLE COMMAND
      FIXHYDANGLE All
      # constrain hydrogen atom bond distances 
      FIXBOND All, element H
      # use a timestep of 2 fs, update the non-bonded interactions every 2 timesteps (as in new MD_simul.mcr example
      CalculationTimeStep = 2
      UpdateCycle         = 2
      TIMESTEP (UpdateCycle),(CalculationTimeStep)
    else
      RAISERROR Illegal lincs settle type definitions defined in script, (LincsSettle) should be either 'On' or 'Off'
         
  # use standard SI units rather than kcal/mol or pound per square inch or stones or other medieval units that are only used in three countries in the world.
  ENERGYUNIT kj/mol
  
  # exit if a warning occurs
  WARNISERROR On
  
  # Start the simulation 
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  SIM On

  # Keep the pressure constant with the solvent density
  PRESSURECTRL Type=SolventProbe,Name=HOH,Density=0.997,Axis=XYZ
  
  # Control the temperature by rescaling the velocities
  TEMPCTRL Rescale

  # save snapshots with an appropriate name at appropriate intervals
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  IntervalsSnapShots = (SnapShotInterval) / ( (CalculationTimeStep)*(UpdateCycle ))
  i = 00000
  FileName = '(MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(i).sim'
  SAVEsim (FileName),(IntervalsSnapShots)

  # Check if the simulation was already done?
  FinalTime = 0+ (( (WarmUpTime)+(EquilibrationTime)+(ProductionTime)) / 1000)

  # *********************** do the WARMING UP or load a snapshot of it after it is done ***************************************************************************************************
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # check if the snapshot at which time the simulation is done is there
  WarmUpDone = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)WarmUpDone.sim
  if (WarmUpDone)
    LOADSIM (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)WarmUpDone.sim
  else
    # the 5 and 10 come from the do while loop below
    NumberOfTemperatureJumps = 0 + (( (AssayTemperature)  - 5) / 10)
    TimeForEquilibrationInTimeSteps = (WarmUpTime) / ( (CalculationTimeStep)*(UpdateCycle ))
    CurrentWaitPeriod = (TimeForEquilibrationInTimeSteps) / (NumberOfTemperatureJumps)

    # start at 5 K, not 0 K, to conserve the existing motions, their speeds are rescaled only, seems like a poor idea if the original speed was 0 K. 
    t = 5

    # The following section sets the temperature with different random seeds. 
    # reliable than the following method. 
    TEMP ((SeedingNumber)+(t))
    # This do while loop increases the temperature, in steps of 10 K to allow the temperature to equilibrate before the next step 
    do
      TEMP (t), REASSIGN=NO
      wait (CurrentWaitPeriod)
      t = (t) + 10
    while (t < (AssayTemperature))
    # save the snapshot of the warmup being done
    SAVESIM  (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)WarmUpDone
    # ensure still saved continuously (unclear if this is necessary)
    SAVESIM (FileName),(IntervalsSnapShots)

    
  # ensure the temperature and time is set correctly now
  TEMP (AssayTemperature), REASSIGN=NO  
  Time (WarmUpTime)
  
  # ********************** do the equilibrationphase or load a snapshot of it after it was done *******************************************************************************************
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  SnapShotEquilibrationDone = 00000 + ( ((WarmUpTime) + (EquilibrationTime) )/(SnapShotInterval) )
  SnapShotWarmUpDone = 00000 + 1+ ( ((WarmUpTime)  )/(SnapShotInterval) )
  EquilibrationDone = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(SnapShotEquilibrationDone).sim
  if (EquilibrationDone)
    LOADSIM (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(SnapShotEquilibrationDone).sim
  else
    LastSnapShot = 0
    CorrectionForSavedSnapShots = 0
    # check which snapshot is there, continue with the last saved snapshot during equilbration
    for j = SnapShotWarmUpDone to SnapShotEquilibrationDone
      EquilibrationSnapshotHalfway = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(j).sim
      if (EquilibrationSnapshotHalfway)
        LastSnapShot = (j)
    if (LastSnapShot)
      LOADSIM (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(LastSnapShot).sim    
      # correct for the remaining waiting time
      CorrectionForSavedSnapshots =  (SnapShotInterval) * ((LastSnapShot) - (SnapShotWarmUpDone))
    EquilibrationTime = (EquilibrationTime) - (CorrectionForSavedSnapShots)
    WaitingTime = (EquilibrationTime) / ( (CalculationTimeStep)*(UpdateCycle ))
    Wait (WaitingTime)
  
  
  # **************************** This is the PRODUCTION run part, sample NACs on the fly **************************************************************************************************
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  # some time tracking calculations and checking whether part of the production has already been done before running the production phase loop
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # How long does a NAC take expressed in intervals rather than fs?
  SamplingNACInterval = 0 + (NACSamplingInterval) / ( (CalculationTimeStep)*(UpdateCycle ))
  
  # How often to sample a NAC?
  TotalNACIntervals = 0 + ((ProductionTime) / (NACSamplingInterval))
    
  # some initialization
  LastSnapShot = (SnapShotEquilibrationDone)
  LastTimeCounter = 0
  
  # Magic number, i.e. time, NACs, and all criteria with their NAC
  MagicNumber = 2 + (2*(count Label_))
  
  # Are there already snapshots and periodically saved tables available for part of the production run??
  # check which snapshot is there, continue with the last saved snapshot during equilbration
  SnapShotProductionDone = 00000 + ( ((WarmUpTime) + (EquilibrationTime)+ (ProductionTime))/(SnapShotInterval) )
  for j = SnapShotEquilibrationDone to SnapShotProductionDone
    EquilibrationSnapshotHalfway = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(j).sim
    if (EquilibrationSnapshotHalfway)
      TimeCounterfs = 0 +  ((j) - (SnapShotEquilibrationDone)) * (SnapShotInterval)
      AlsoATableFile = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(j).tab
      if (AlsoATableFile)
        LastSnapShot = (j)
        LastTimeCounter = TimeCounterfs
    else
      break      
  if (LastSnapShot > SnapShotEquilibrationDone)
    LOADSIM (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(LastSnapShot).sim
 
  # This keeps track of how far we are into production time
  TimeCounterfs = 0 + (LastTimeCounter)
  StartInterval = 000001 + ((TimeCounterfs)/(NACSamplingInterval))
  # the counter is to coincide with the snapshots, the snapshots for the table are saved independently
  AtTimeSnapShot = (SnapShotInterval) * (NACShotEveryXSnapshots)/ (NACSamplingInterval)
  AlreadyObtainedSnapShotsShouldBe = (StartInterval) / (AtTimeSnapShot)
  CounterToTimeSnapShot = 0 + ((StartInterval) - ( (AlreadyObtainedSnapShotsShouldBe) * (AtTimeSnapShot) )) - 1
  # do the actual sampling, this loop runs during the production time, when it is finished it saves a final table but during the snapshots tables are saved as well
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  for j =  StartInterval to TotalNACIntervals 
    PRINT (j) interval out of (TotalNACIntervals)
    if (CounterToTimeSnapShot == AtTimeSnapShot)
      CurrentSnapShot = (SnapShotEquilibrationDone) + ((TimeCounterfs)/(SnapShotInterval))
      CounterToTimeSnapShot = 0
      # also save a table with only the new NACs 
      for q = (000000+ ((j) - (AtTimeSnapShot))) to ((j) - 1)
        # TABULATE The results
        TABULATE (0.001 *(q)*(NACSamplingInterval) )
        TABULATE (NAC(q))
        for c = 1 to count Label_
          TABULATE ((Label_(c))_Q_(q))
          TABULATE ((Label_(c))_NAC(q))
      SAVETAB 1, (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(CurrentSnapShot),Format=Text,Columns=(MagicNumber),NumFormat=%12.6f, Mixed_Data 
      DELTAB 1
    CounterToTimeSnapShot = (CounterToTimeSnapShot) + 1
    wait (SamplingNACInterval)
    TimeCounterfs  = (TimeCounterfs) + (NACSamplingInterval)
    # do the measurements
    for c = 1 to count Label_
      (Label_(c))_Q_(j)    = ((Label_(c))_Measurement) 
    # determine whether NAC is achieved 
    NAC(j) = 1
    for c = 1 to count Label_
      (Label_(c))_NAC(j) = 0.00000
      if (  ((Label_(c))_Q_(j))  > ((Label_(c))_Q_Min))
        if (((Label_(c))_Q_(j))  < ((Label_(c))_Q_Max))
          (Label_(c))_NAC(j) = 1.00000
      NAC(j) = (NAC(j)) * ((Label_(c))_NAC(j))  
     
  # stop the simulations and save a final scene and intermediate tab file
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
  # ensure last snapshots saved first for the tab file
  CurrentSnapShot = (SnapShotEquilibrationDone) + ((TimeCounterfs)/(SnapShotInterval))
  # Only save it if it is not already there (!) 
  LastTabFileAlreadyExists = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(CurrentSnapShot).tab
  if (LastTabFileAlreadyExists == 0)
    for q = ((j) - (AtTimeSnapShot) + 1) to ((j) )
      # TABULATE The results
      TABULATE (0.001 *(q)*(NACSamplingInterval) )
      TABULATE (NAC(q))
      # First all the distances
      for c = 1 to count Label_
        TABULATE ((Label_(c))_Q_(q))
        TABULATE ((Label_(c))_NAC(q))
    SAVETAB 1, (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(CurrentSnapShot),Format=Text,Columns=(MagicNumber),NumFormat=%12.6f, Mixed_Data 
    DELTAB 1
  

  # and here for the sim snapshot
  wait 2
  
  SIM Off
  # How long has it taken in total)
  FinalTime = 0+ (( (WarmUpTime)+(EquilibrationTime)+(ProductionTime)) / 1000)
  


# ========================================================================================================================================================================================
# =========================  Part D of the Script: The Analysis of the MD simulation(s)===================================================================================================
# ========================================================================================================================================================================================


# ***************************** This Loop analyses the NAC data for every seed individually *********************************************************************************************
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
for CSN = 01 to (MaximumSeeds)  
  # save two tables, one with the original sampled data, one with the averaged data and their standard deviations
  # ------------------------------------------------
  
  # first load the original data from the intermediate tables
  CounterToTimeSnapShot = (AtTimeSnapShot)
  PreviousSnapShot      = (SnapShotEquilibrationDone)
  PreviousTimeCounter   = 0
  q = 000001
  for j =  000001 to TotalNACIntervals
    if (CounterToTimeSnapShot == AtTimeSnapShot)
      CounterToTimeSnapShot   = 0
      q                       = 000001
      LastSnapShot            = (PreviousSnapShot) + (NACShotEveryXSnapshots) 
      LastTimeCounter         = (PreviousTimeCounter) + (NACShotEveryXSnapshots) * (SnapShotInterval)
      LOADTAB (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(LastSnapShot).tab  
      CounterToTimeSnapShot = 0
      PreviousSnapShot      = (LastSnapShot)    
      PreviousTimeCounter   = (LastTimeCounter) 
      RawData() = TAB 1
      DELTAB 1
    # now convert this loaded table into the lists 
    
    TotalNacIntervalsAlreadyAnalyzed = 00000 + (count RawData)/ (MagicNumber)
    NACSamplingIntervalList(j) = 1000.000 *(RawData( ((MagicNumber)*((q)-1)) +1))
    NAC(j)                     =           (RawData( ((MagicNumber)*((q)-1)) +2))
    for c = 1 to count Label_
      (Label_(c))_Q_(j)          =           (RawData( ((MagicNumber)*((q)-1)) +((2*(c))+1) ))
      (Label_(c))_NAC(j)         =           (RawData( ((MagicNumber)*((q)-1)) +((2*(c))+2) ))
    # to see how far it is in the log file, if it is not simply stock
    if MoreLOGFile == 'TRUE'
      print reading in (q) out of (TotalNacIntervalsAlreadyAnalyzed) in this tab file
      print reading in (j) out of (TotalNacIntervals) in total for this MD run
    # update counter
    CounterToTimeSnapShot = (CounterToTimeSnapShot) + 1   
    q                     = (q) + 1
  
  # Save a table with all the data. 
  for j = 000001 to TotalNACIntervals
    # TABULATE The results
    TABULATE (0.001 *(j)*(NACSamplingInterval) )
    TABULATE (NAC(j))
    for c = 1 to count Label_
      TABULATE ((Label_(c))_Q_(j))
      TABULATE ((Label_(c))_NAC(j))
  LabelCollection = '____Interval ________NAC?'
  for c = 1 to count Label_
    LabelCollection = '(LabelCollection) (Label_(c)) CriteriaMet?'
  SAVETAB 1,(MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(EquilibrationTime)fs_(ProductionTime)fs_NAC_DATA_ProductionRun,Format=Text,Columns=(MagicNumber),NumFormat=%12.2f,(LabelCollection)
  DELTAB 1

  # TABULATE The interval results, see if even distribution of NACs over the intervals or not. 
  if InterValAnalysis == 'On'
    for q = 1 to count TrajectorySubAnalysisIntervals
      CountTrajectoryInterVal = 1
      CurrentTimeInterValEnd = 0
      EndTrajectoryCount =  (TrajectorySubAnalysisIntervals(q))/ (NACSamplingInterval) 
      for j = 000001 to TotalNACIntervals
        if (CountTrajectoryInterval ==  EndTrajectoryCount)
          CurrentTimeInterValEnd = (CurrentTimeInterValEnd) + (TrajectorySubAnalysisIntervals(q))  
          TABULATE (CurrentTimeInterValEnd)
          CurrentTrajectoryAverage  = (100.0000000000000 * (mean TemporaryIntervalNAC))
          TABULATE (CurrentTrajectoryAverage)
          # print Current number of intervals is (count TemporaryIntervalNAC)
          CountTrajectoryInterVal = 1
        else 
          CountTrajectoryInterVal = (CountTrajectoryInterVal)+ 1
        TemporaryIntervalNAC(CountTrajectoryInterVal) = (NAC(j))
      SAVETAB 1,(MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(EquilibrationTime)fs_(ProductionTime)fs_NAC_Interval(TrajectorySubAnalysisIntervals(q))fsProductionRun,Format=Text,Columns=2,NumFormat=%12.2f,____Interval _____NACPerc  
      DELTAB 1
    
  TABULATE 'Total_Number'
  TABULATE (TotalNACIntervals)
  TABULATE 'NA'
  TABULATE 'NA'
  
  TABULATE 'Total_NACs'
  TABULATE (sum NAC)
  TABULATE 'NA'
  TABULATE 'NA'
  
  TABULATE 'TotalTime_ps'
  TABULATE ((ProductionTime)/1000)
  TABULATE 'NA'
  TABULATE 'NA'
  
  TABULATE 'NACPercentag'
  TABULATE (100.000*(mean NAC))
  TABULATE 'NA'
  TABULATE 'NA'
  
  for c = 1 to count Label_
    TABULATE '(Label_(c))'
    TABULATE (mean (Label_(c))_Q_)
    TABULATE (stddev (Label_(c))_Q_)
    TABULATE (100.000*(mean (Label_(c))_NAC))
  SAVETAB 1,(MacroTarget)_LS(LincsSettle)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_NAC_Results_Seed(CSN)_Summary,Format=Text,Columns=4,NumFormat=%12.3f,___Parameter Average_Value Standard_Dev PassCritPerc
  DELTAB 1


# This calculates the averages and standard deviations of the ensemble rather than that of the individual runs
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
for CSN = 01 to (MaximumSeeds)
  # loat the table and convert it to a temporary list
  LOADTAB (MacroTarget)_LS(LincsSettle)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_NAC_Results_Seed(CSN)_Summary
  IntermediateList() = TAB 1
  DELTAB 1
  # convert the variables
  PRINT round (CSN) of CSN = 01 to (MaximumSeeds)
  TotalNumbers(CSN)    = (IntermediateList(2))
  TotalNacs(CSN)       = (IntermediateList(6))
  TotalTimes(CSN)      = (IntermediateList(10))
  NACPercentages(CSN)  = (IntermediateList(14))
  # this part has not recently (2026-03-12) been tested for functionality
  for c = 1 to count Label_
    (Label_(c))s(CSN)  = (IntermediateList(12 + (c*4)))
    (Label_(c))SD(CSN) = (IntermediateList(12 + (c*4)+1))
    (Label_(c))Pass(CSN) = (IntermediateList(12 + (c*4)+2))    
  PRINT round (CSN) of CSN = 01 to (MaximumSeeds)




# Make a table with all the numbers, both mean and standard deviation

TABULATE 'Total_Number'
TABULATE (mean TotalNumbers)
if MS == 1
  TABULATE 'NA'
else 
  TABULATE (stddev TotalNumbers)
for i = 1 to 4 
  TABULATE 'NA'

TABULATE 'Total_NACs'
TABULATE (mean   TotalNacs)
if MS == 1
  TABULATE 'NA'
else 
  TABULATE (stddev TotalNacs)
for i = 1 to 4 
  TABULATE 'NA'

TABULATE 'TotalTime_ps'
TABULATE (mean     TotalTimes)
if MS == 1
  TABULATE 'NA'
else 
  TABULATE (stddev   TotalTimes)
for i = 1 to 4 
  TABULATE 'NA'

TABULATE 'NACPercentag'
TABULATE (mean       NACPercentages)  
if MS == 1
  TABULATE 'NA'
else 
  TABULATE (stddev     NACPercentages)  
for i = 1 to 4 
  TABULATE 'NA'

for c = 1 to count Label_
  TABULATE '(Label_(c))'
  TABULATE (mean     (Label_(c))s)    
  if MS == 1
    TABULATE 'NA'
  else 
    TABULATE (stddev   (Label_(c))s)   
  TABULATE (mean     (Label_(c))SD) 
  if MS == 1
    TABULATE 'NA'
  else 
    TABULATE (stddev   (Label_(c))SD) 
  TABULATE (mean (Label_(c))Pass)   
  if MS == 1
    TABULATE 'NA'
  else 
    TABULATE (stddev   (Label_(c))Pass)  

SAVETAB 1, (MacroTarget)_LS(LincsSettle)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_NAC_Results_All(MS)Seeds_Summary,Format=Text,Columns=7,NumFormat=%12.3f,___Parameter Average_Value SD_from average  mean SD Standard_Dev PercCritPass Standard_Dev
DELTAB 1



# *************** Now also do the analysis of the snapshots: record energy and RMSD and prepare an averaged structure of the complex with Bfactors (during the equilibration time)****** 
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


LOADSCE (MacroTarget)_SaltWaterBoxNotMinimized.sce
# Analyze the Bfactors of the crystal structure (or at least the earlier structure)
ResidueListProtein() = LISTRES protein, FORMAT=RESNUM
BfactorsOriginal()   = BFACTORATOM Res Protein atom CA
for i = 1 to count BfactorsOriginal
  RMSFOriginal(i) = SQRT (BfactorsOriginal(i)*0.037995443)
# Make a copy to use later for RMSD calculations
OriginalStructure = DUPLICATEOBJ 1
REMOVEOBJ (OriginalStructure)

ENERGYUNIT kj/mol

# In case we want to obtain snapshots
if MakeSnapshotPdbs == 'Yes'
  for CSN = 01 to (MaximumSeeds)
    for i = 1 to count SnapShotPdbsTimeps 
      TimeInFs = (1000* (SnapShotPdbsTimeps(i)))
      LOADSIM (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(00000+(TimeInFs)/(SnapShotInterval))
      SAVEPDB OBJ 1, (MacroTarget)_(CSN)_(0000+(SnapShotPdbsTimeps(i)))ps

  
for CSN = 01 to (MaximumSeeds)  
  i = 00000
  NextSnapShotExists = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(i).sim
  while (NextSnapShotExists)
    LOADSIM (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(i).sim
    CurrentTime = TIME
    TimeSnapShot(i) =(CurrentTime) / 1000
    EnergySnapShot(i) = ENERGY
    SIM OFF
    
    
    # Now analyse the rmsd from the crystal structure or designed structure
    ADDOBJ (OriginalStructure)
    RMSDCASnapShot(i)       = SUPATOM OBJ 1 atom CA, OBJ (OriginalStructure) atom CA
    RMSDBackBoneSnapShot(i) = SUPATOM OBJ 1 atom backbone, OBJ (OriginalStructure) atom backbone
    RMSDAllHeavySnapShot(i) = SUPATOM OBJ 1 protein atom element !H, OBJ (OriginalStructure) protein atom element !H
    REMOVEOBJ (OriginalStructure)
    # if during production phase
    if ((TimeSnapShot(i)*1000) >= ((WarmUpTime)+(EquilibrationTime)))
      ADDPOSATOM OBJ 1
      # if also want to predict binding energy
    i = i + 1
    NextSnapShotExists = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(i).sim
  # save a table with the energy and RMSD versus time
  for i = 00000 to ((count TimeSnapShot) - 1)
    TABULATE (TimeSnapShot(i))
    TABULATE (EnergySnapShot(i))
    TABULATE (RMSDCASnapShot(i))
    TABULATE (RMSDBackBoneSnapShot(i))
    TABULATE (RMSDAllHeavySnapShot(i))
  # store the Final RMSD of CA snapshot
  EndRMSDCA(CSN) = (RMSDCASnapShot(00000 + ((count TimeSnapShot) - 1)))
  # ensure to reset them in case some snapshots are missing
  for i = 00000 to ((count TimeSnapShot) - 1)
    TimeSnapShot(i)         = 0
    EnergySnapShot(i)       = 0
    RMSDCASnapShot(i)       = 0
    RMSDBackBoneSnapShot(i) = 0
    RMSDAllHeavySnapShot(i) = 0
  SAVETAB 1,(MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_(CSN)_LS(LincsSettle)_EnergyRMSDTime,Format=Text,Columns=5,NumFormat=%12.3f,____Timeinps _Energykjmol ______RMSDCA ____RMDSBack ___RMSDHeavy 
  DELTAB 1

  

# Take the average of the positions to save a PDB and YOB File of it
AVERAGEPOSATOM OBJ 1

# make a list of the CA RMSD values
LongListRMSFMDRun()   = RMSFATOM OBj 1,UNIT=A
LongListAtoms()       = LISTATOM OBJ 1
j = 1
k = 1
for i = 1 to count LongListRMSFMDRun
  CurrentAtomName = NAMEATOM Atom (LongListAtoms(i))
  if (CurrentAtomName == 'CA')
    ShortListRMSFMDRun(j) = (LongListRMSFMDRun(i))
    j = j + 1
  CurrentElement = ElementAtom (LongListAtoms(i))
  if (CurrentElement > 1.5)
    HeavyAtomsRMSFMDRun(k) = (LongListRMSFMDRun(i))
    k = k + 1
        

# Take the average of the positions to save a PDB and YOB File of it
if UseRMSFInPdbYOBFile== 'TRUE'
  for i = 1 to count LongListRMSFMDRun
    BFACTORATOM (i), (LongListRMSFMDRun(i))
  SAVEPDB OBJ 1, (MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_All(MS)Seeds_LS(LincsSettle)_AvgRMSF
else
  RMSFATOM OBj 1,UNIT=bfactor
  for i = 1 to count LongListRMSFMDRun
    Value10timesTooHigh =  BFACTORATOM (i)
    if Value10timesTooHigh < 9999
      BFACTORATOM (i), ((Value10timesTooHigh)/10)
    else
      BFACTORATOM (i),999.9
  SAVEPDB OBJ 1, (MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_All(MS)Seeds_LS(LincsSettle)_Avg




# save a table suitable for a plot of Bfactor versus residue number
for i = 1 to count ResidueListProtein
  TABULATE (ResidueListProtein(i))
  TABULATE (RMSFOriginal(i))
  TABULATE (ShortListRMSFMDRun(i))
SAVETAB 1,(MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_All(MS)Seeds_LS(LincsSettle)_RMSF,Format=Text,Columns=3,NumFormat=%12.3f,__ResnNumber RMSFOriginal ___RMSFMDrun 
DELTAB 1


# ******************** if desired, get the deviation from the original structure of the average structure versus time, average of all the seeds *****************************************
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if CalculateAverageTrajectory == 'On'
  for i = 00000 to ((count TimeSnapShot) - 1)
    CLEAR
    GreenFlag = 1
    LOADSCE (MacroTarget)_SaltWaterBox.sce
    OriginalStructure = DUPLICATEOBJ 1
    REMOVEOBJ (OriginalStructure)
    for CSN = 01 to (MaximumSeeds)
      NextSnaphotExists = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(i).sim
      if (NextSnaphotExists)
        LOADSIM (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(i).sim
        CurrentTime = TIME
        TimeSnapShot(i) =(CurrentTime) / 1000
        SIM Off
        ADDOBJ (OriginalStructure)
        SUPATOM OBJ 1 protein atom element !H, OBJ 4 protein atom element !H
        ADDPOSATOM OBJ 1
        REMOVEOBJ (OriginalStructure)
      else 
        GreenFlag = 0
    if (GreenFlag)  
      AVERAGEPOSATOM OBJ 1
      ADDOBJ (OriginalStructure)
      RMSDCASnapShot(i)       = SUPATOM OBJ 1 atom CA, OBJ 4 atom CA
      RMSDBackBoneSnapShot(i) = SUPATOM OBJ 1 atom backbone, OBJ 4 atom backbone
      RMSDAllHeavySnapShot(i) = SUPATOM OBJ 1 protein atom element !H, OBJ 4 protein atom element !H
    else
      TimeSnapShot(i)         = 0
      RMSDCASnapShot(i)       = 0
      RMSDBackBoneSnapShot(i) = 0
      RMSDAllHeavySnapShot(i) = 0
  for i = 00000 to ((count TimeSnapShot) - 1)
    TABULATE (TimeSnapShot(i))
    TABULATE (RMSDCASnapShot(i))
    TABULATE (RMSDBackBoneSnapShot(i))
    TABULATE (RMSDAllHeavySnapShot(i))
  SAVETAB 1,(MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_All(MS)Seeds_LS(LincsSettle)_RMSDversusTime,Format=Text,Columns=4,NumFormat=%12.3f,____Timeinps ______RMSDCA ____RMDSBack ___RMSDHeavy 
  DELTAB 1



# ******************************* The end of the script, leave automatically if set ******************************************************************************************************
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 


if AutomaticExitAtEnd == 'YES'
  PRINT Finished the script succesfully. 
  EXIT
else
  SHOWMESSAGE Finished with analysis trajectories 
  CONSOLE ON





