#
#   ConversionSelectedDesigns2NamedAndHydratedPdbs.mcr
#
#   Hein J. Wijma, University of Groningen, h.j.wijma@rug.nl
#  
#   This collects all the selected variants in named and hydrated PDB files
#   Please verify by inspection that all the cofactors were put back correctly


########### To be set by user ####################################################################################################################################################
##################################################################################################################################################################################

# the name of the list with selected mutations
DesiredList         = 'list_SelectedMutations'

# the first part of the name of all the pdb files to be read
BasicName           = '1NWW_cleaned'
 
# the total number of subdirectories
TotalSubdirectories = 3

# The pdb file where the water molecules can be found
HydrationShellFrom = '1NWW_cleaned'

########### Parts below only need to be modified if there are errors in the hydrated designs that prevent proper MD simulations ##################################################
##################################################################################################################################################################################


# make a new directory
shell mkdir NamedPdbFiles

# immediately save the template file in there
shell mkdir NamedPdbFiles/Subdir_template
shell cp (HydrationShellFrom).pdb NamedPdbFiles/Subdir_template/

# This loads the list of desired mutations, the list is in the format mutation, energy change, error from different calculations
LOADTAB (DesiredList)
DesiredListRaw() = TAB 1
DELTAB 1
TotalDesired = (count DesiredListRaw)
TotalDesired = (TotalDesired) / 3
for i = 1 to (TotalDesired)
  print (i)
  DesiredList(i) = '(DesiredListRaw((3 * ((i)- 1))+1))'
  print (DesiredList(i))

# This speeds up the calculation, no LOG file output
CONSOLE Off

#This loop goes through the directories
for h = 1 to (TotalSubdirectories)
  CurrentSubdirectory = 'Subdirectory(h)'
  # This generates a list of mutations of the subdirectory in the right format for YASARA to read
  shell echo ' ' > ListYasaraMutations.tab 
  shell cat (CurrentSubdirectory)/List_Mutations_readable.txt >> ListYasaraMutations.tab
  LOADTAB ListYasaraMutations.tab
  #shell rm ListYasaraMutations.tab
  AllMutationsInDirectory()= TAB 1
  DELTAB 1
  
  
  # determine which mutations can be found in current directory
  NumberInDirectory = (count AllMutationsInDirectory)
  NumberInDirectory = (NumberInDirectory) / 4
  for i = 1 to (NumberInDirectory)
    CurrentNumberInDirectory(i)  = (0+(AllMutationsInDirectory((4 * ((i)- 1))+1)))
    CurrentMutationInDirectory(i)= '(AllMutationsInDirectory((4 * ((i)- 1))+2))(0+(AllMutationsInDirectory((4 * ((i)- 1))+3)))(AllMutationsInDirectory((4 * ((i)- 1))+4))'
  
  # now check if something identical found
  for i = 1 to TotalDesired
    for j = 1 to NumberInDirectory
      if '(DesiredList(i))' == '(CurrentMutationInDirectory(j))'
        # Load the pdb files of the target and the pdb file from which the hydration shell is coming.  
        LOADPDB (CurrentSubdirectory)/(BasicName)_(CurrentNumberInDirectory(j))_0.pdb
        CLEANALL
        OPTHYDALL
        PRINT carried out optimize h-bonding network.
        LOADPDB (HydrationShellFrom)
        
        # This alignment method is needed since otherwise renumbering gives problems even if identical sequences
        ALIGNOBJ 2,1, method = globalseq
        
        # The following gets water from one side to the other, and saves the results and report how many waters were deleted
        DELRES OBJ 2 RES !HOH
        ObjectNum = CountObj All
        if (ObjectNum) > 1
          JOINOBJ 2,1
        A = COUNTRES HOH
        DELRES HOH atom element O with distance < 1.5 from res protein atom element O C N P S
        B = COUNTRES HOH
        PRINT Deleted (A - B) waters from original (A) waters of structure since protein was within 1.5 Angstrom now
        
        # This ensures metals are bound to their ligands
        ListMetals() = LISTATOM element Fe Co Zn Cu Ni 
        for j = 1 to count ListMetals
          ListLigands() = Listatom atom element N S O with distance < 3.0 from (ListMetals(j))
          for k= 1 to count ListLigands
            ShowArrow Start=AtAtom,(ListLigands(k)),End=AtAtom,(ListMetals(j)), heads = 0 
        
        
        # WOW means with original Water
        MakeNewDirectory = FILESIZE NamedPdbFiles/Subdir_(DesiredList(i))/(BasicName)_(DesiredList(i))_WOW.pdb
        if !(MakeNewDirectory)
          SHELL mkdir NamedPdbFiles/Subdir_(DesiredList(i))
        OPTHYDALL
        SAVEPDB OBJ 1, NamedPdbFiles/Subdir_(DesiredList(i))/(BasicName)_(DesiredList(i))_WOW
        DELOBJ All

exit
