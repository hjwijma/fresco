#
#
#   DisulfideDiscovery.mcr
#
#  This script is used to design disulfide bonds in proteins
# It features
# - looking for suitable positions 
# - testing all 9 possible geometries of a disulfide bond
# - verification that both the energy and the geometry is OK
# - provide multiple conformations of a single SS bond if multiple conformations are similarly stable
#
#  example use: yasara -txt DisulfideDiscovery.mcr "MacroTarget = '1NWW'" "TableFile = '1NWW_AllResidues.tab'" >LOG&



# USER settings, could be modified by user 
# ========================================

# Only calculate if at least this distance in between, no natural reason to use (most disulfide bonds in proteins bridge very few residues, could be set to 0, instead of 15)
MinimalInBetweenDistance  = 0

# Make list of residues with sidechain < X from disulfide bonds (suggested value = 7)
MinimalSidechainDistance  = 7

# if SkipEliminations is set to On,the it would just build all disulfide bonds. Only useful if previously  reported disulfide bonds are not reproduced by the software. 
# Preferentially use with very limited TableFile (otherwise pile ofPDB files saved).  
SkipGeometricCriteria     = 'Off' 
SkipEnergyCriteria        = 'Off' 

# standard setting of energy of 10 kJ/mol for dihedrals and angles of SS atoms seems to work well to discern between feasible and unfeasible bonds, together with geometric criteria.
MaximalEnergySSBOND       = 10.00

# number of processors to use
PROCESSORS 4

# General Settings, should not be touched
# ==========================================

# These number are from Pellequer, Chen, Proteins, 2006,65,192-202
Min_CB_Distance           = (3.180 )
Max_CB_Distance           = (4.780 )
Min_CA_Distance           = (3.720 )
Max_CA_Distance           = (6.770 )

# These numbers are derived from plots by Petersen, 1999
MinimalSSDistance         = 1.970
MaximalSSDistance         = 2.090
LeftTurnMinimalCADistance = 4.200
MinimalCSSangle           =  99.0
MaximalCSSangle           = 112.0
MinimalCCSangle           = 107.0
MaximalCCSangle           = 121.0
LeftTurnMinimalChi3       = -120.0
LeftTurnMaximalChi3       = -60.0
RighTurnMinimalChi3       = 70.0
RighTurnMaximalChi3       = 130.0
Chi1WrongZoneLow1         = -140.0
Chi1WrongZoneHigh1        = -110.0
Chi1WrongZoneLow2         = -30.0
Chi1WrongZoneHigh2        = 30.0
Chi1WrongZoneLow3         = 100.0
Chi1WrongZoneHigh3        = 160.0
Chi2WrongZoneLow          = -20.0
Chi2WrongZoneHigh         = +40.0



# Here the calculation is set up
# ==========================================

# Load the protein, ensure to maintain original coordinates
# --------------------------------------------
TestFilePresent =FILESIZE (MacroTarget).pdb
if !TestFilePresent
  RAISEERROR (MacroTarget).pdb not found!
LOADPDB (MacroTarget),Center=NO
STYLE Stick

# load the table file and verify that the loaded residues exist
# -----------------------
LOADTAB (TableFile)
TemporaryList() = TAB 1
DELTAB 1
TotalRows = count TemporaryList
TotalRows = (0 + (TotalRows)) / 3
for i = 1 to TotalRows
  ResType(i)   = '(TemporaryList( (((i)-1) *3)+1))'
  ResChain(i)  = '(TemporaryList( (((i)-1) *3)+2))'
  ResNumber(i) = (0+(TemporaryList( (((i)-1) *3)+3)))

# convert the Restypes into 3 letter codes
for i = 1 to count ResType
  if (ResType(i) == 'A')
    ResType(i) = 'ALA' 
  if (ResType(i) == 'R')
    ResType(i) = 'ARG' 
  if (ResType(i) == 'N')
    ResType(i) = 'ASN' 
  if (ResType(i) == 'D')
    ResType(i) = 'ASP' 
  if (ResType(i) == 'C')
    ResType(i) = 'CYS' 
  if (ResType(i) == 'Q')
    ResType(i) = 'GLN' 
  if (ResType(i) == 'E')
    ResType(i) = 'GLU' 
  if (ResType(i) == 'G')
    ResType(i) = 'GLY' 
  if (ResType(i) == 'H')
    ResType(i) = 'HIS' 
  if (ResType(i) == 'I')
    ResType(i) = 'ILE' 
  if (ResType(i) == 'L')
    ResType(i) = 'LEU' 
  if (ResType(i) == 'K')
    ResType(i) = 'LYS' 
  if (ResType(i) == 'M')
    ResType(i) = 'MET' 
  if (ResType(i) == 'F')
    ResType(i) = 'PHE' 
  if (ResType(i) == 'P')
    ResType(i) = 'PRO' 
  if (ResType(i) == 'S')
    ResType(i) = 'SER' 
  if (ResType(i) == 'T')
    ResType(i) = 'THR' 
  if (ResType(i) == 'W')
    ResType(i) = 'TRP' 
  if (ResType(i) == 'Y')
    ResType(i) = 'TYR' 
  if (ResType(i) == 'V')
    ResType(i) = 'VAL' 

# check if they exist
for i = 1 to (TotalRows)
  TestCondition = countres (ResType(i)) (ResChain(i)) (ResNumber(i))
  if !(TestCondition)
    RAISEERROR (ResType(i)) (ResChain(i)) (ResNumber(i)) does not exist
  
# Make a PolyGly copy, will be used repeatedly
DUPLICATEOBJ 1
REMOVEOBJ 1
DELRES !protein
SWAPRES All, Gly
CLEANALL
REMOVEOBJ 2

# set up energy calculations
ENERGYUNIT kj/mol
FORCEFIELD AMBER03

# set up standard table of donor dihedrals to try
DD1 = -60
DD2 = -60
DD3 = -60
DD4 = 60
DD5 = 60
DD6 = 60
DD7 = 180
DD8 = 180
DD9 = 180

# set up standard table of acceptro dihedrals to try.
AD1 = -60
AD2 =  60
AD3 = 180
AD4 = -60
AD5 =  60
AD6 = 180
AD7 = -60
AD8 =  60
AD9 = 180



# here the real calculation starts
# =================================

# switching off the log file speeds up the calculation
Console Off

# This keeps track of total stored conformations
StoredConformationsTotal = 0

for i = 1 to count ResNumber
  # make a fresh copy of the polyglycine scaffold 
  ADDOBJ 2
  DUPLICATEOBJ 2
  REMOVEOBJ 2
  # get a list of residues (complete with chain IDs) within Carbon Alpha distance and test the other four criteria one by one for them. 
  # reinit the lists to empty
  NearbyResList()   = 0
  NearbyChainList() = 0
  NearbyResList()  = Listres res atom CA with distance < (Max_CA_Distance) from res (Resnumber(i)) atom CA, FORMAT = RESNUM
  NearbyChainList()= Listres res atom CA with distance < (Max_CA_Distance) from res (Resnumber(i)) atom CA, FORMAT = MOLNAME
  # mutate them to ALA
  SWAPRES protein with distance < (Max_CA_Distance) from res (Resnumber(i)) atom CA, Ala
  SWAPRES res (Resnumber(i)), Ala
  for j = 1 to count NearbyResList
    #initialize that it does not pass the criteria
    CriterionPassed(j) = 0
    # check that part of the list of residues load (otherwise also other residues that are not allowed to mutate would be tried) 
    for k = 1 to count ResNumber
      if (ResNumber(k)) == (NearbyResList(j))
        if ('(ResChain(k))' == '(NearbyChainList(j))')
          # check if already had it before, we do not want to make twice exactly the same SSbonds
          AtomNumber1 = LISTATOM MOL (ResChain(i))        res (ResNumber(i)) atom CA
          AtomNumber2 = LISTATOM MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CA
          Criterion = (AtomNumber2) - (AtomNumber1)
          if (Criterion > 1)
            Criterion = Distance MOL (ResChain(i))        res (ResNumber(i)) atom CA,  MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CA
            if Criterion > Min_CA_Distance
              Criterion = Distance MOL (ResChain(i))        res (ResNumber(i)) atom CB,  MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CB
              if Criterion > Min_CB_Distance
                if Criterion < Max_CB_Distance
                  if ('(ResChain(i))' != '(NearbyChainList(j))')
                    CriterionPassed(j) = 1
                  if ((NearbyResList(j)) - (ResNumber(i))) > (MinimalInBetweenDistance)
                    CriterionPassed(j) = 1
  for j = 1 to count NearbyResList  
    if (CriterionPassed(j))
      # print to the log file
      print The combination (ResChain(i)) (ResNumber(i)) with (NearbyChainList(j)) (NearbyResList(j)) passes the CA en CB distance criteria and the minimal number of residues in between, now try mutagenesis with energy minimization
      # mutate them back to GLY
      SWAPRES All, Gly
      # make the mutation
      SWAPRES MOL (ResChain(i)) protein res (ResNumber(i)), Cys
      SWAPRES MOL (NearbyChainList(j)) protein res (NearbyResList(j)), Cys
      # since will be used 9 times, make a backup
      DUPLICATEOBJ 3
      REMOVEOBJ 3
      
      # ThiskeepsTrackOf the Number of stored conformations for this round
      StoredConformationsCurrentResidueCombination = 0
      # calculate for 9 different start conformations of a disulfide bond whether it passes energetic and geometric criteria 
      # --------------------------------------------------------------------------------------------------------------------
      for k = 1 to 9
        # prevent that the backbone is modified in position
        FIXATOM backbone
        # set the donor and acceptor dihedrals, see list above
        Dihedral MOL (ResChain(i)) res (ResNumber(i)) atom C, MOL (ResChain(i)) res (ResNumber(i)) atom CA,MOL (ResChain(i)) res (ResNumber(i)) atom CB,MOL (ResChain(i)) res (ResNumber(i)) atom SG, SET = (DD(k))
        Dihedral MOL (NearbyChainList(j)) res (NearbyResList(j)) atom C, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CA, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CB, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom SG, SET=(AD(k))
        # connect the residues
        ADDBOND MOL (ResChain(i)) res (ResNumber(i)) atom SG,  MOL (NearbyChainList(j)) res (NearbyResList(j)) atom SG, update = yes
        AtomNumber1 = LISTATOM MOL (ResChain(i))        res (ResNumber(i)) atom CA
        AtomNumber2 = LISTATOM MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CA
        #CELL auto, extension = 2.0, shape=cuboid, res atom (AtomNumber1) (AtomNumber2)
        CELL auto, extension = 5.0, shape=cuboid, res atom (AtomNumber1) (AtomNumber2)
        Boundary Wall
        Experiment Minimization
        EXPERIMENT ON
        WAIT EXPEND
        # store the energy results of the calculation
        INTERACTIONS Bond,Angle,Dihedral,Planarity
        EnergyList() = ENERGYATOM SG
        INTERACTIONS All
        CurrentEnergy(k) = SUM EnergyList
        # if energy criterium passed, continue
        # ------------------------------------
        if ( (CurrentEnergy(k)) < (MaximalEnergySSBOND) or (SkipEnergyCriteria   == 'On')  )
          # write to LOG file
          print Did pass the criterium with (0.00 + (CurrentEnergy(k))) SS bond energy, from residue (ResNumber(i)) to residue (NearbyResList(j))
          
          # store the results (energy, dihedrals, type of SS bond (classified per dihedral class)
          CurrentChi1Donor        = DIHEDRAL MOL (ResChain(i))        res (ResNumber(i))     atom N,  MOL (ResChain(i))        res (ResNumber(i))     atom CA, MOL (ResChain(i))        res (ResNumber(i))     atom CB, MOL (ResChain(i))        res (ResNumber(i))     atom SG
          CurrentChi2Donor        = DIHEDRAL MOL (ResChain(i))        res (ResNumber(i))     atom CA, MOL (ResChain(i))        res (ResNumber(i))     atom CB, MOL (ResChain(i))        res (ResNumber(i))     atom SG, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom SG
          DihedralSS              = DIHEDRAL MOL (ResChain(i))        res (ResNumber(i))     atom CB, MOL (ResChain(i))        res (ResNumber(i))     atom SG, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom SG, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CB
          CurrentChi2Acceptor     = DIHEDRAL MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CA, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CB, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom SG, MOL (ResChain(i))        res (ResNumber(i))     atom SG 
          CurrentChi1Acceptor     = DIHEDRAL MOL (NearbyChainList(j)) res (NearbyResList(j)) atom N,  MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CA, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CB, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom SG
          
          DistanceSS              = DISTANCE MOL (ResChain(i))        res (ResNumber(i))     atom SG, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom SG 
          DistanceCA              = DISTANCE MOL (ResChain(i))        res (ResNumber(i))     atom CA, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CA 
          DistanceCB              = DISTANCE MOL (ResChain(i))        res (ResNumber(i))     atom CB, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CB 
          
          CurrentAngleCCSDonor    = ANGLE MOL (ResChain(i))           res (ResNumber(i))     atom CA, MOL (ResChain(i))        res (ResNumber(i))     atom CB, MOL (ResChain(i))        res (ResNumber(i))     atom SG
          CurrentAngleCCSAcceptor = ANGLE MOL (NearbyChainList(j))    res (NearbyResList(j)) atom CA, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CB, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom SG
          CurrentAngleCSSDonor    = ANGLE MOL (ResChain(i))           res (ResNumber(i))     atom CB, MOL (ResChain(i))        res (ResNumber(i))     atom SG, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom SG
          CurrentAngleCSSAcceptor = ANGLE MOL (NearbyChainList(j))    res (NearbyResList(j)) atom CB, MOL (NearbyChainList(j)) res (NearbyResList(j)) atom SG, MOL (ResChain(i))        res (ResNumber(i))     atom SG
          
          CurrentChi1Donor        = 0.000 + (CurrentChi1Donor)
          CurrentChi2Donor        = 0.000 + (CurrentChi2Donor)
          CurrentChi2Acceptor     = 0.000 + (CurrentChi2Acceptor)
          CurrentChi1Acceptor     = 0.000 + (CurrentChi1Acceptor)
          
          CurrentAngleCCSDonor    = 0.000 + (CurrentAngleCCSDonor   )
          CurrentAngleCCSAcceptor = 0.000 + (CurrentAngleCCSAcceptor)
          CurrentAngleCSSDonor    = 0.000 + (CurrentAngleCSSDonor   )
          CurrentAngleCSSAcceptor = 0.000 + (CurrentAngleCSSAcceptor)
          
          DistanceSS = 0.000 + (DistanceSS)
          DihedralSS = 0.000 + (DihedralSS)
          DistanceCA = 0.000 + (DistanceCA)       
          
          
          # Now determine if it passes the geometric criteria, partly different for left and right turning dihedrals (Chi3, SS  bond)
          # -------------------------------------------------------------------------------------------------------------------------
          
          # initalize that it passes to 1, set to zero if exceed or meets a criterion
          PassesGeometricCriteria = 1
            
          # length disulfide bond
          if DistanceSS < MinimalSSDistance
            PassesGeometricCriteria = 0
          if DistanceSS > MaximalSSDistance
            PassesGeometricCriteria = 0
          
          # stricter distances criterium for left hand SS bridges
          if DihedralSS < 0
            if DistanceCA < LeftTurnMinimalCADistance
              PassesGeometricCriteria = 0
          
          # angles in CCS and CSS
          if CurrentAngleCSSDonor < MinimalCSSangle
            PassesGeometricCriteria = 0
          if CurrentAngleCSSAcceptor < MinimalCSSangle
            PassesGeometricCriteria = 0
          if CurrentAngleCSSDonor > MaximalCSSangle
            PassesGeometricCriteria = 0
          if CurrentAngleCSSAcceptor > MaximalCSSangle
            PassesGeometricCriteria = 0               
          if CurrentAngleCCSDonor < MinimalCCSangle
            PassesGeometricCriteria = 0
          if CurrentAngleCCSAcceptor < MinimalCCSangle
            PassesGeometricCriteria = 0
          if CurrentAngleCCSDonor > MaximalCCSangle
            PassesGeometricCriteria = 0
          if CurrentAngleCCSAcceptor > MaximalCCSangle
            PassesGeometricCriteria = 0               
          
          # criteria for left and right handed disulfide bridges
          if DihedralSS < 0
            if DihedralSS < LeftTurnMinimalChi3
              PassesGeometricCriteria = 0
            if DihedralSS > LeftTurnMaximalChi3
              PassesGeometricCriteria = 0
          if DihedralSS > 0
            if DihedralSS < RighTurnMinimalChi3 
              PassesGeometricCriteria = 0
            if DihedralSS > RighTurnMaximalChi3
              PassesGeometricCriteria = 0
          
          # for chi1 three different zones that are illegal
          if CurrentChi1Donor > Chi1WrongZoneLow1
            if CurrentChi1Donor < Chi1WrongZoneHigh1
              PassesGeometricCriteria = 0
          if CurrentChi1Donor > Chi1WrongZoneLow2
            if CurrentChi1Donor < Chi1WrongZoneHigh2
              PassesGeometricCriteria = 0
          if CurrentChi1Donor > Chi1WrongZoneLow3 
            if CurrentChi1Donor < Chi1WrongZoneHigh3
              PassesGeometricCriteria = 0
          if CurrentChi1Acceptor > Chi1WrongZoneLow1 
            if CurrentChi1Acceptor < Chi1WrongZoneHigh1
              PassesGeometricCriteria = 0
          if CurrentChi1Acceptor > Chi1WrongZoneLow2 
            if CurrentChi1Acceptor < Chi1WrongZoneHigh2
              PassesGeometricCriteria = 0
          if CurrentChi1Acceptor > Chi1WrongZoneLow3
            if CurrentChi1Acceptor < Chi1WrongZoneHigh3
              PassesGeometricCriteria = 0
          
          # for chi2 lots more or less permitted, with the exception of:
          if CurrentChi2Donor > Chi2WrongZoneLow 
            if CurrentChi2Donor < Chi2WrongZoneHigh
              PassesGeometricCriteria = 0
          if CurrentChi2Acceptor > Chi2WrongZoneLow 
            if CurrentChi2Acceptor < Chi2WrongZoneHigh
              PassesGeometricCriteria = 0
          
          if (SkipGeometricCriteria == 'On')
            PassesGeometricCriteria = 1
            
          if PassesGeometricCriteria
            # print to LOG File
            print from residue (ResNumber(i)) to residue (NearbyResList(j)) measured are Chi1Donor: (CurrentChi1Donor), Chi2Donor: (CurrentChi2Donor), Chi3SSbond:(DihedralSS), Chi2Acceptor: (CurrentChi2Acceptor), Chi1Acceptor: (CurrentChi1Acceptor), distance: (DistanceSS), and energy:(0.000+( CurrentEnergy(k)))
            print from residue (ResNumber(i)) to residue (NearbyResList(j)) measured are CCS_angle_Donor: (CurrentAngleCCSDonor), CSS_angle_Donor: (CurrentAngleCSSDonor), CSS_angle_Acceptor:(CurrentAngleCSSAcceptor), CCS_angle_Acceptor: (CurrentAngleCCSAcceptor), and Calpha distance: (DistanceCA)
            
            
            # Classify its dihedral angles, not store same conformation twice
            # Chi1 donor, -60,+60,+180
            Chi1DonorStringCurrent = 'Undecided'
            if CurrentChi1Donor > -120
              if CurrentChi1Donor < 0
                Chi1DonorStringCurrent = 'm60'
            if CurrentChi1Donor >= 0
              if CurrentChi1Donor < 120
                Chi1DonorStringCurrent = 'p60'
            if Chi1DonorStringCurrent == 'Undecided' 
              Chi1DonorStringCurrent = 'p180'      
            # Chi1 acceptor, -60, +60,+180
            Chi1AcceptorStringCurrent = 'Undecided'
            if CurrentChi1Acceptor > -120
              if CurrentChi1Acceptor < 0
                Chi1AcceptorStringCurrent = 'm60'
            if CurrentChi1Acceptor >= 0
              if CurrentChi1Acceptor < 120
                Chi1AcceptorStringCurrent = 'p60'
            if Chi1AcceptorStringCurrent == 'Undecided' 
              Chi1AcceptorStringCurrent = 'p180'      
            # Chi2 donor, -90, +90
            if  CurrentChi2Donor < 0
              Chi2DonorStringCurrent = 'm90'
            else
              Chi2DonorStringCurrent = 'p90'  
            # Chi2 acceptor, -90,+90 
            if  CurrentChi2Acceptor < 0
              Chi2AcceptorStringCurrent = 'm90'
            else
              Chi2AcceptorStringCurrent = 'p90'  
            # chi3, -90, or +90
            if DihedralSS < 0
              Chi3StringCurrent = 'm90'  
            else
              Chi3StringCurrent = 'p90'
            SummaryString = '(Chi1DonorStringCurrent)(Chi2DonorStringCurrent)(Chi3StringCurrent)(Chi2AcceptorStringCurrent)(Chi1AcceptorStringCurrent)'
            print from residue (ResNumber(i)) to residue (NearbyResList(j)) measures are summarized (SummaryString) with energy (0.000 +(CurrentEnergy(k)))  
            
            # Check if already stored an disulfide bonded, for this round, if so compare, otherwise make an extra object
            StoreCurrentSSbond = 'No'
            
            if StoredConformationsCurrentResidueCombination >= 0
              # check for two conditions, Found and equal that is better (if no, act is if no SS bound found for this combination yet)
              # and check for equal but that one is worse than the current one, then act to delete that structure and have it replaced
              FoundEqual = 'No'
              FoundEqualButWorse  = 'No'
              
              # loop through all stored combinations to see if one is identical and if so if the current one is lower in energy
              for l = (1 + (StoredConformationsTotal) - (StoredConformationsCurrentResidueCombination)) to (StoredConformationsTotal)
                print stored combinations should be (StoredConformationsTotal) , for this combination (StoredConformationsCurrentResidueCombination), and l is (l)
                if '(SummaryString)' == '(CentralListResiduesSummaryString(l))'  
                  FoundEqual = 'Yes'  
                  if  CentralListEnergy(l) > (CurrentEnergy(k))
                    print Decided to replace conformation with (CentralListEnergy(l)) energy by conformation with (CurrentEnergy(k)) energy
                    FoundEqualButWorse ='Yes'
                    # assign object and numbers
                    CurrentProteinObject = (l) + 5
                    CurrentConformationNumber = (l)
               
              if FoundEqual == 'No'
                StoreCurrentSSbond = 'Yes'
                StoredConformationsCurrentResidueCombination = (StoredConformationsCurrentResidueCombination) + 1
                StoredConformationsTotal = (StoredConformationsTotal) + 1
                # Assign the object number and the current conformation Number
                CurrentProteinObject      = (StoredConformationsTotal) + 5
                CurrentConformationNumber = (StoredConformationsTotal)
              
              if FoundEqualButWorse == 'Yes'
                StoreCurrentSSbond = 'Yes'
                
                # Assign the object number and the current conformation Number, happened really above !, here just for clarity
                CurrentProteinObject      = (CurrentProteinObject      )
                CurrentConformationNumber = (CurrentConformationNumber )
                
                # Remove the current object
                DELOBJ (CurrentProteinObject      )
              
            if StoredConformationsCurrentResidueCombination == 0
              # keep track of the numbers
              StoredConformationsCurrentResidueCombination = (StoredConformationsCurrentResidueCombination) + 1
              StoredConformationsTotal = (StoredConformationsTotal) + 1
              
              # Assign the object number and the current conformation Number
              CurrentProteinObject = (StoredConformationsTotal) + 5
              CurrentConformationNumber = (StoredConformationsTotal)
              
              # This signals that the Current SS bond needs to be stored
              StoreCurrentSSbond = 'Yes'
            
            if StoreCurrentSSbond == 'Yes'
              # make a copy of the disulfide bond, inside the normal protein, not in the glycine scaffold. 
              # ------------------------------------------------------------------------------------------
              
              # make a copy of the original protein
              ADDOBJ 1
              DUPLICATEOBJ 1
              REMOVEOBJ 1  
              # mutate the residues
              SWAPRES OBJ (CurrentProteinObject) MOL (ResChain(i))        protein res (ResNumber(i)),     Cys
              SWAPRES OBJ (CurrentProteinObject) MOL (NearbyChainList(j)) protein res (NearbyResList(j)), Cys
              
              # set their conformation and make the disulfide bond
              ANGLE    OBJ (CurrentProteinObject) MOL (ResChain(i))        res (ResNumber(i))     atom CA, OBJ (CurrentProteinObject) MOL (ResChain(i))        res (ResNumber(i))     atom CB, OBJ (CurrentProteinObject) MOL (ResChain(i))        res (ResNumber(i))     atom SG       , Set =( CurrentAngleCCSDonor     ) 
              ANGLE    OBJ (CurrentProteinObject) MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CA, OBJ (CurrentProteinObject) MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CB, OBJ (CurrentProteinObject) MOL (NearbyChainList(j)) res (NearbyResList(j)) atom SG       , Set =( CurrentAngleCCSAcceptor  ) 
              
              DIHEDRAL OBJ (CurrentProteinObject) MOL (ResChain(i))        res (ResNumber(i))     atom N,  OBJ (CurrentProteinObject) MOL (ResChain(i))        res (ResNumber(i))     atom CA, OBJ (CurrentProteinObject) MOL (ResChain(i))        res (ResNumber(i))     atom CB, OBJ (CurrentProteinObject) MOL (ResChain(i))        res (ResNumber(i))     atom SG, Set =(CurrentChi1Donor   )  
              DIHEDRAL OBJ (CurrentProteinObject) MOL (NearbyChainList(j)) res (NearbyResList(j)) atom N,  OBJ (CurrentProteinObject) MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CA, OBJ (CurrentProteinObject) MOL (NearbyChainList(j)) res (NearbyResList(j)) atom CB, OBJ (CurrentProteinObject) MOL (NearbyChainList(j)) res (NearbyResList(j)) atom SG, Set =(CurrentChi1Acceptor)  
              
              ADDBOND  OBJ (CurrentProteinObject) MOL (ResChain(i))        res (ResNumber(i))     atom SG, OBJ (CurrentProteinObject) MOL (NearbyChainList(j)) res (NearbyResList(j)) atom SG, update = yes
               
              # Also record the original residue at donor and acceptor positions 
              ADDOBJ 1
              OriginalDonor    = LISTRES OBJ 1 MOL (ResChain(i))        res (ResNumber(i))     , FORMAT=RESNAME1
              OriginalAcceptor = LISTRES OBJ 1 MOL (NearbyChainList(j)) res (NearbyResList(j)) , FORMAT=RESNAME1
              REMOVEOBJ 1 
              
              # store all of their data in a table
              CentralListChi1Donor(CurrentConformationNumber)                  = (CurrentChi1Donor)
              CentralListChi2Donor(CurrentConformationNumber)                  = (CurrentChi2Donor)
              CentralListChi2Acceptor(CurrentConformationNumber)               = (CurrentChi2Acceptor)
              CentralListChi1Acceptor(CurrentConformationNumber)               = (CurrentChi1Acceptor)
              
              CentralListAngleCCSDonor(CurrentConformationNumber)              = (CurrentAngleCCSDonor   )
              CentralListAngleCCSAcceptor(CurrentConformationNumber)           = (CurrentAngleCCSAcceptor)
              CentralListAngleCSSDonor(CurrentConformationNumber)              = (CurrentAngleCSSDonor   )
              CentralListAngleCSSAcceptor(CurrentConformationNumber)           = (CurrentAngleCSSAcceptor)
              
              CentralListDistanceSS(CurrentConformationNumber)                 = (DistanceSS)
              CentralListDihedralSS(CurrentConformationNumber)                 = (DihedralSS)
              CentralListDistanceCA(CurrentConformationNumber)                 = (DistanceCA)
              CentralListDistanceCB(CurrentConformationNumber)                 = (DistanceCB)
              CentralListResiduesSummaryString(CurrentConformationNumber)      = '(SummaryString)'
              CentralListResiduesInformation(CurrentConformationNumber)        = '(ResChain(i))_(OriginalDonor)(ResNumber(i))C_(NearbyChainList(j))_(OriginalAcceptor)(NearbyResList(j))C_(SummaryString)' 
              CentralListResiduesInformationShort(CurrentConformationNumber)   = '(ResChain(i))(OriginalDonor)(ResNumber(i))C(NearbyChainList(j))(OriginalAcceptor)(NearbyResList(j))C' 
                            
              CentralListEnergy(CurrentConformationNumber)                     = (0.000 + (CurrentEnergy(k)))
              
              CentralListNamePdbFile(CurrentConformationNumber)                = '(MacroTarget)_(CentralListResiduesInformation(CurrentConformationNumber))'
              
              CentralListMoleculeDonor(CurrentConformationNumber)              = '(ResChain(i))'
              CentralListMoleculeAcceptor(CurrentConformationNumber)           = '(NearbyChainList(j))' 
              CentralListResDonor(CurrentConformationNumber)                   = (ResNumber(i))
              CentralListResAcceptor(CurrentConformationNumber)                = (NearbyResList(j))
              
              
              # And remove the object to prevent problems
              REMOVEOBJ  (CurrentProteinObject)     
                             
          else
            print from residue (ResNumber(i)) to residue (NearbyResList(j)) did not pass all geometric criteria
        else 
          # write to LOG file
          print Did NOT pass the criterium with (0.00 + (CurrentEnergy(k))) SS bond energy from residue (ResNumber(i)) to residue (NearbyResList(j))
        DELOBJ 4
        ADDOBJ 3
        if ( k != 9)
          DUPLICATEOBJ 3
          REMOVEOBJ 3
  DELOBJ 3 


SHELL mkdir disulfideBonds_(MacroTarget)__(TableFile)

# Currently Object 6 till .. hold all the pdb files with a disulfide bond. Now 
CONSOLE Off
MAKETAB 1
for i = 1 to StoredConformationsTotal
  # Right ObjectNumber
  CurrentObjectNumber = 5 + (i)
  REMOVEOBJ All
  ADDOBJ (CurrentObjectNumber)
  Style Ribbon
  
  # save the pdb file

  SAVEPDB OBJ (CurrentObjectNumber), disulfideBonds_(MacroTarget)__(TableFile)/(CentralListNamePdbFile(i)),transform = no 
    
  # Save all the data in a table
  SELECTTAB 1
  TABULATE (CentralListDistanceCA(i))
  TABULATE (CentralListDistanceCB(i))
  TABULATE (CentralListDistanceSS(i))
  
  TABULATE (CentralListAngleCCSDonor(i))    
  TABULATE (CentralListAngleCCSAcceptor(i)) 
  TABULATE (CentralListAngleCSSDonor(i))    
  TABULATE (CentralListAngleCSSAcceptor(i)) 
  
  TABULATE (CentralListChi1Donor(i))      
  TABULATE (CentralListChi1Acceptor(i))   
  TABULATE (CentralListChi2Donor(i))      
  TABULATE (CentralListChi2Acceptor(i))   
  
  TABULATE (CentralListDihedralSS(i))

  TABULATE (CentralListEnergy(i))

  TABULATE '(CentralListNamePdbFile(i)).pdb'
  
# Save the resulting table file
SELECTTAB 1
TABULATE 'd = distance ad = angle donor, aa = angle acceptor, hd = dihedral donor, ha = dihedral acceptor, hhSS = dihedral of SS bond'
SAVETAB 1, disulfideBonds_(MacroTarget)__(TableFile)/(MacroTarget)_ConformationsDisulfideBonds, columns=14,NumFormat=10.3f,'   d_CA_CA    d_CB_CB    d_SG_SG adCA_CB_SG aaCA_CB_SG adCB_SG_SG aaCB_SG_SG    hd_Chi1    ha_Chi1    hd_Chi2    ha_Chi2      hh_SS  Energy_SS Name_Pdb_File  '
DELTAB 1

DELOBJ 5
ADDOBJ All
for i = 1 to StoredConformationsTotal
  # Right ObjectNumber
  CurrentObjectNumber = 5 + (i)
  RENAMEOBJ (CurrentObjectNumber), (CentralListResiduesInformationShort(i))

ZOOMRES All  
CENTERRES All
RENAMEOBJ 1, WT
RENAMEOBJ 2, PolyGly
REMOVEOBJ 1 2
SAVESCE disulfideBonds_(MacroTarget)__(TableFile)/(MacroTarget)_AllDisulfides

CONSOLE Off
print Succesfully reached end of DisulfideDiscovery.mcr
print '============================================'
exit
