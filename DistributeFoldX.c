/***********************************************************************************************************************************************************
*                                                                                                                                                          *
*                                                                                                                                                          *
*                DistributeFoldX.c                                                                                                                         *
*                =================                                                                                                                         *
*                                                                                                                                                          *
*  Author: Hein Wijma, University of Groningen, last updated 12th of March 2026, h.j.wijma@rug.nl                                                          *
*                                                                                                                                                          *
*  Can be used to automate the calculation of a very large number of mutations by FoldX                                                                    *
*                                                                                                                                                          *
*  FoldX can be obtained at the website foldx.crg.es                                                                                                       *
*                                                                                                                                                          *
*                                                                                                                                                          *
*                                                                                                                                                          *
*                                                                                                                                                          *
***********************************************************************************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "DistributeFoldX.h"

int main(int argc, char *argv[])
{
  /* first some important declarations */
  char NamePDBFile[200], NameTableFile[200];           /* names of the PDB file to move to the separate directories and of the table file with mutations to open */
  char strand1[4], strand2[4], strand3[4], strand4[4], strand5[4], strand6[4], strand7[4], strand8[4];  /* program for maximally an octameric protein */
  int NumberOfStrands, NumberOfMutationsPerDirectory;  /* obvious */
  char WhichStageAreWeIn[20];                          /* check if it is Phase1 or Phase2 */
  float EnergyCutOff = 0;                              /* in Phase2, make also a list of all the solutions with an energy better than this in kJ mol -1 */ 
  char FoldXApplication[500];                          /* Which FoldX version to use */
  
  /* provide the user with info about what the program does and examples of its use */
  short_info();
  example_use();
  
  
  if (DEBUG) /* list input */ 
  {  int i;
     for (i = 0; i < argc; i++)
     {  
        printf("%s\n",argv[i]);
     }
  }
  
  /* first check, number of arguments should be > 6 */
  if (argc < 7)
     error_message("Too few arguments, see above");
  
  
  /* ==========================================assign all the char in the command line to what they should be ================================================*/
  /* =========================================================================================================================================================*/
  strcpy (WhichStageAreWeIn, argv[1]);
  if (DEBUG)
    printf("The name of the stage we are in as read in is %s\n", WhichStageAreWeIn);
  
  
  strcpy(NamePDBFile, argv[2]);
  if (DEBUG)
    printf("The name of the PDB file as read in is %s\n", NamePDBFile);
 
 
  NumberOfStrands = (int) atof(argv[3]);
  if (DEBUG)
    printf("The number of strands read in is %d\n", NumberOfStrands);
  
  
  if (NumberOfStrands == 1)
  {
    strcpy(strand1,argv[4]);
    if (DEBUG)
      printf("The letter of this strand is %s\n", strand1);
  }  
  
  
  if (NumberOfStrands == 2)
  {
    strcpy(strand1,argv[4]);
    strcpy(strand2,argv[5]);
    if (DEBUG)
      printf("The letters of these strand are %s, %s\n", strand1, strand2);
  }  
  
  
  if (NumberOfStrands == 3)
  {
    strcpy(strand1,argv[4]);
    strcpy(strand2,argv[5]);
    strcpy(strand3,argv[6]);
    if (DEBUG)
      printf("The letters of these strand are %s, %s, %s\n", strand1, strand2, strand3);
  }  
  
  if (NumberOfStrands == 4)
  {
    strcpy(strand1,argv[4]);
    strcpy(strand2,argv[5]);
    strcpy(strand3,argv[6]);
    strcpy(strand4,argv[7]);
    if (DEBUG)
      printf("The letters of these strand are %s, %s, %s, %s\n", strand1, strand2, strand3, strand4);
  }  
   
  if (NumberOfStrands == 5)
  {
    strcpy(strand1,argv[4]);
    strcpy(strand2,argv[5]);
    strcpy(strand3,argv[6]);
    strcpy(strand4,argv[7]);
    strcpy(strand5,argv[8]);
    if (DEBUG)
      printf("The letters of these strand are %s, %s, %s, %s, %s\n", strand1, strand2, strand3, strand4, strand5);
  } 

  if (NumberOfStrands == 6)
  {
    strcpy(strand1,argv[4]);
    strcpy(strand2,argv[5]);
    strcpy(strand3,argv[6]);
    strcpy(strand4,argv[7]);
    strcpy(strand5,argv[8]);
    strcpy(strand6,argv[9]);
    if (DEBUG)
      printf("The letters of these strand are %s, %s, %s, %s, %s, %s\n", strand1, strand2, strand3, strand4, strand5, strand6);
  }


  if (NumberOfStrands == 7)
  {
    strcpy(strand1,argv[4]);
    strcpy(strand2,argv[5]);
    strcpy(strand3,argv[6]);
    strcpy(strand4,argv[7]);
    strcpy(strand5,argv[8]);
    strcpy(strand6,argv[9]);
    strcpy(strand7,argv[10]);
    if (DEBUG)
      printf("The letters of these strand are %s, %s, %s, %s, %s, %s, %s\n", strand1, strand2, strand3, strand4, strand5, strand6, strand7);
  }


  if (NumberOfStrands == 8)
  {
    strcpy(strand1,argv[4]);
    strcpy(strand2,argv[5]);
    strcpy(strand3,argv[6]);
    strcpy(strand4,argv[7]);
    strcpy(strand5,argv[8]);
    strcpy(strand6,argv[9]);
    strcpy(strand7,argv[10]);
    strcpy(strand8,argv[11]);
    if (DEBUG)
      printf("The letters of these strand are %s, %s, %s, %s, %s, %s, %s, %s\n", strand1, strand2, strand3, strand4, strand5, strand6, strand7, strand8);
  }









  strcpy(NameTableFile, argv[4 + NumberOfStrands]);    
  if (DEBUG)
    printf("The name of the table file is %s\n", NameTableFile);
  NumberOfMutationsPerDirectory = (int) atof(argv[5 + NumberOfStrands]);
  if (DEBUG)
    printf("The number of mutations per directory is %d\n", NumberOfMutationsPerDirectory);
  
  
  if (strcmp(WhichStageAreWeIn,"Phase1")== 0)
  {
    strcpy(FoldXApplication, argv[6 + NumberOfStrands]);
    printf("The name and location of FoldX is %s\n", FoldXApplication);
  }
  
  if (strcmp(WhichStageAreWeIn,"Phase2")== 0)
  { EnergyCutOff = (float) atof(argv[6 + NumberOfStrands]);
     printf("The energy cutoff is %.6f", EnergyCutOff);
  }   
  /* ===========================  read in the table of residues and if it is a multimer, ensure the residues are found in every strand defined above =============================*/
  /* =============================================================================================================================================================================*/
    
  
  /* load the sequence */
  struct SequenceProtein * RawSequence = LoadYasaraTable(NameTableFile);
   
  /*print out the table of residues to ensure it was properly read (only if debug)*/
  if (DEBUG)
    PrintSequenceProtein(RawSequence);
  
  /* purge the sequence of residues that are not present in all chains, just copy the raw list if just one strand*/ 
  struct SequenceProtein * PurgedSequence;
  
  if (NumberOfStrands == 1)
  {
     PurgedSequence = RawSequence;  
  }
  if (NumberOfStrands == 2)
  {
    PurgedSequence = GetPurgedListResidues(RawSequence, strand1);
  }
    if (NumberOfStrands == 3)
  {
    PurgedSequence = GetPurgedListResidues(RawSequence, strand1);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand2);
  }
  if (NumberOfStrands == 4)
  {
    PurgedSequence = GetPurgedListResidues(RawSequence, strand1);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand2);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand3);
  }
  if (NumberOfStrands == 5)
  {
    PurgedSequence = GetPurgedListResidues(RawSequence, strand1);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand2);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand3);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand4);
  }
  if (NumberOfStrands == 6)
  {
    PurgedSequence = GetPurgedListResidues(RawSequence, strand1);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand2);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand3);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand4);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand5);
  }
  if (NumberOfStrands == 7)
  {
    PurgedSequence = GetPurgedListResidues(RawSequence, strand1);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand2);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand3);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand4);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand5);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand6);
  }
  if (NumberOfStrands == 8)
  {
    PurgedSequence = GetPurgedListResidues(RawSequence, strand1);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand2);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand3);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand4);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand5);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand6);
    PurgedSequence = GetPurgedListResidues(PurgedSequence, strand7);
  }

  if (DEBUG)
    PrintSequenceProtein(PurgedSequence); 

  /*report to LOG file */
  printf("\nFinished eliminating residues that are not present in all of the reported %d strands\n", NumberOfStrands);


  /* ===========================  make the list of mutations and write them to a series of files that can be used by FoldX =======================================================*/
  /* =============================================================================================================================================================================*/
  
  struct MutationProtein * MutatedProteinList = SaturationScanning(PurgedSequence);
  
  if (DEBUG)
    PrintMutationsProtein(MutatedProteinList);
  
  int NumberOfMutations = CountNumberMutations(MutatedProteinList);
  
  /*report to LOG file */
  printf("\nThe number of mutations is %d\n", NumberOfMutations);
  
  
  /* decide on how many subdirectories to make */
  int NumberOfSubdirectories = ((NumberOfMutations) / (NumberOfMutationsPerDirectory));
  if ( ( (NumberOfMutations) % (NumberOfMutationsPerDirectory)) != 0)
    NumberOfSubdirectories++;

  /*report to LOG file */
  printf("\nWith maximum %d mutations per subdirectory this requires %d directories\n", NumberOfMutationsPerDirectory, NumberOfSubdirectories);

  if (strcmp(WhichStageAreWeIn,"Phase1")== 0)
  { /* --------------------------------------------------------here starts a huge if loop for step 1---------------------------------------- */
    /* things beloware not done as a sub function since otherwise strand1, strand2, strand3, and strand 4 needed to be passed in some dynamic way that is too much time for me */
    
    /* make the subdirectories and write the mutations lists to them */
    MakeSubdirectories(NumberOfSubdirectories);
     
    char inbetweenstring[1000];
  
    /*write the sequence mutations to every one of them */
    FILE *mutationfile;
    FILE *readablemutationfile;
    FILE *listPDBfile;
    FILE *todolist;
    char CurrentDirectoryName[200];
    char PdbListFileName[400];
    char MutationFileName[400];
    char ReadableMutationFileName[400];
    struct MutationProtein * p = MutatedProteinList; 
    int h, i;
  
    /* initialize the file todolist, if already present it will be overwritten */
    todolist = fopen("todolist", "w");
    
    /* the following writes mutation files to the different directories and the additional files */
    for (h = 0; h < NumberOfSubdirectories; h++)
    {
      /*first get all the file names used later */
      
      /* the name of the current subdirectory */
      strcpy(CurrentDirectoryName, "Subdirectory");
      titoa((h + 1), inbetweenstring);
      strcat(CurrentDirectoryName, inbetweenstring);
      strcat(CurrentDirectoryName, "/");
      
      /* the file names */     
      strcpy(MutationFileName,CurrentDirectoryName);
      strcpy(ReadableMutationFileName,CurrentDirectoryName);
      strcpy(PdbListFileName, CurrentDirectoryName);
      
      strcat(MutationFileName,"individual_list.txt"); /*FoldX does accept this name as the input list for mutations, for example mutationlist is not accepted */
      strcat(ReadableMutationFileName,"List_Mutations_readable.txt");
      strcat(PdbListFileName, "list.txt");
      
      
        /* fill first the mutation files */
        mutationfile = fopen(MutationFileName, "w");
        readablemutationfile = fopen(ReadableMutationFileName, "w");
        for (i = 0; i < NumberOfMutationsPerDirectory && p != NULL; i++, p = p->next)
          {  
            fprintf(mutationfile,"%s%s%d%c",p->CurrentAminoAcid,strand1,p->residueNumber,p->NewAminoAcid);
          fprintf(readablemutationfile, "%d  %s  %d  %c\n", (i + 1),p->CurrentAminoAcid,p->residueNumber,p->NewAminoAcid);
          if (NumberOfStrands >= 2)
            fprintf(mutationfile,",%s%s%d%c",p->CurrentAminoAcid,strand2,p->residueNumber,p->NewAminoAcid);
          if (NumberOfStrands >= 3)
          fprintf(mutationfile,",%s%s%d%c",p->CurrentAminoAcid,strand3,p->residueNumber,p->NewAminoAcid);
          if (NumberOfStrands >= 4)
            fprintf(mutationfile,",%s%s%d%c",p->CurrentAminoAcid,strand4,p->residueNumber,p->NewAminoAcid);
          if (NumberOfStrands >= 5)
            fprintf(mutationfile,",%s%s%d%c",p->CurrentAminoAcid,strand5,p->residueNumber,p->NewAminoAcid);
          if (NumberOfStrands >= 6)
            fprintf(mutationfile,",%s%s%d%c",p->CurrentAminoAcid,strand6,p->residueNumber,p->NewAminoAcid);
          if (NumberOfStrands >= 7)
            fprintf(mutationfile,",%s%s%d%c",p->CurrentAminoAcid,strand7,p->residueNumber,p->NewAminoAcid);
          if (NumberOfStrands >= 8)
            fprintf(mutationfile,",%s%s%d%c",p->CurrentAminoAcid,strand8,p->residueNumber,p->NewAminoAcid);
          fprintf(mutationfile, ";\n"); 
        }
      fclose(mutationfile);
      fclose(readablemutationfile);
      
      /* and then the list file which contains the PDB file */
      
      listPDBfile = fopen(PdbListFileName, "w");
      fprintf(listPDBfile,"%s.pdb\n\n",NamePDBFile); 
      fclose (listPDBfile);
      
      /* and also copy necessary files from the starting directory to the subdirectory */
     
      /*strcpy(inbetweenstring, "cp rotabase.txt ");
      strcat(inbetweenstring, CurrentDirectoryName);
      system(inbetweenstring);*/

      strcpy(inbetweenstring, "cp ");
      strcat(inbetweenstring, NamePDBFile);
      strcat(inbetweenstring, ".pdb ");
      strcat(inbetweenstring, CurrentDirectoryName);
      system(inbetweenstring);
      
      /* add to the file todolist the following 3 lines */
      fprintf(todolist,"cd %s\n%s --command=BuildModel --pdb=%s.pdb  --mutant-file=individual_list.txt --numberOfRuns=5 > LOG&\ncd ..\n",CurrentDirectoryName, FoldXApplication,NamePDBFile); 
    }
    /* close the file todolist and make it executable */
    fclose(todolist);
    system("chmod +x todolist");
   
   
  /* here end a large loop for if this is step 1 */
  }
   
  /* ============================================ Here collect the results of the calculations, again not in seperate function due to the variable number of strands======================= */
  if (strcmp(WhichStageAreWeIn,"Phase2")== 0)
    {
  printf("\n========= This software assume the user has waited for the FOLDX calculations to finish, e.g. type ps aux | grep FoldX ============\n");
   
  /* Collect the data, analyse it and write it to the right files, write it in kJ mol instead of the kcal mol -1*/
  /* -----------------------------------------------------------------------------------------------------------*/
  
  
  /* for each of the subdirectories, read in the file Average_BuildModel_?????????.fxout, store the results */
  FILE *readablemutationfile;
  FILE *foldxenergyfile;
  char CurrentDirectoryName[200], inbetweenstring[400];
  char ReadableMutationFileName[400];
  char FoldXEnergyFileName[400];
  char c; /* check for end of file EOF */
  char TargetNameToLookFor[400];
  int i, h, interestingpart;

  /* initialize the structure of the mutation energy list*/
  struct MutationProteinEnergy *Before = NULL;
  struct MutationProteinEnergy * MutationsEnergyList; 
  
  
  /* read in the energies from the different subdirectories */ 
  for (h = 0; h < NumberOfSubdirectories; h++)
  {
    /*first get all the file names used later */
    
    /* the name of the current subdirectory */ 
    strcpy(CurrentDirectoryName, "Subdirectory");
    titoa((h + 1), inbetweenstring);
    strcat(CurrentDirectoryName, inbetweenstring);
    strcat(CurrentDirectoryName, "/");
    
    /* the file names */      
    strcpy(ReadableMutationFileName,CurrentDirectoryName);
    strcpy(FoldXEnergyFileName, CurrentDirectoryName);
    strcat(ReadableMutationFileName,"List_Mutations_readable.txt");
    strcat(FoldXEnergyFileName, "Average_");
    strcat(FoldXEnergyFileName, NamePDBFile);
    strcat(FoldXEnergyFileName, ".fxout");
    
    /* test whether the FoldX energy file exists */
    printf("Currently looking for %s....\n", FoldXEnergyFileName);
    foldxenergyfile = fopen(FoldXEnergyFileName, "r");
    if (foldxenergyfile == '\0')
      error_message("can't find it or empty file "); 
    else
      printf("Found it, start analyzing now\n");
    
    /* idem for the file with the mutations */
    readablemutationfile = fopen(ReadableMutationFileName, "r");
    if  (readablemutationfile == '\0')
      error_message("can't find the list with mutations  "); 
    
    /* decide whether there is a useful part to read in or not, basiccaly the foldxfiles start with a huge index, followed by now and then by a piece of filename, SD, averageenergy */
    interestingpart = FALSE;
    i = 1;
    strcpy(TargetNameToLookFor, NamePDBFile);
    strcat(TargetNameToLookFor, "_");
    titoa((i), inbetweenstring);
    strcat(TargetNameToLookFor, inbetweenstring);
    fscanf(foldxenergyfile," %s", inbetweenstring);
    do
    { 
      /*check whether OK to read */
      if (strcmp(TargetNameToLookFor,inbetweenstring)==0)
        interestingpart = TRUE;
      if (interestingpart == TRUE)
      {  /* get the data */
         MutationsEnergyList = malloc(sizeof(struct MutationProteinEnergy));
         fscanf(foldxenergyfile, "%f %f",&MutationsEnergyList->EnergyChangeSD, &MutationsEnergyList->EnergyChange);/* here the standard deviation and the energy value */
         fscanf(readablemutationfile, "%d %s %d %c",&MutationsEnergyList->IndexInDirectory, MutationsEnergyList->CurrentAminoAcid, &MutationsEnergyList->residueNumber, &MutationsEnergyList->NewAminoAcid);/* here the mutations */ 
          
         /* convert the energies to kJ mol -1 */
         MutationsEnergyList->EnergyChangeSD = (float) 4.1840000 * MutationsEnergyList->EnergyChangeSD;
         MutationsEnergyList->EnergyChange   = (float) 4.1840000 * MutationsEnergyList->EnergyChange;
         
         /* give an error message if lost the frame somehow */
         if (MutationsEnergyList->IndexInDirectory != i)
           error_message("index problem");
         /* prepare for the next part of the list */
         MutationsEnergyList->next  =  Before;   
         Before                     = MutationsEnergyList;      
         
         /* prepare for next part of the cycle*/
         i++;
         strcpy(TargetNameToLookFor, NamePDBFile);
         strcat(TargetNameToLookFor, "_");
         titoa((i), inbetweenstring);
         strcat(TargetNameToLookFor, inbetweenstring);

         interestingpart = FALSE;
      }     
      
      fscanf(foldxenergyfile," %s", inbetweenstring);
      
    } while((c = fgetc(foldxenergyfile)) != EOF);
    
    fclose(foldxenergyfile);
    fclose(readablemutationfile);
    
  }    
     
  if (DEBUG)
    PrintMutationsProteinEnergy(MutationsEnergyList); 
  
  
  /* get the list of mutations below cuttoff energy */
  struct MutationProteinEnergy *MutationsEnergyListBelowCutOff = CuttOffMutationsEnergy(MutationsEnergyList, EnergyCutOff);
  
  /* get the list of best mutations per position */
  struct MutationProteinEnergy *MutationsEnergyListBestPerPosition =KeepBestMutationsEnergy(MutationsEnergyList);
  
  /* get the list of best mutations per position below cutoff */
  struct MutationProteinEnergy *MutationsEnergyListBestPerPositionBelowCutOff = CuttOffMutationsEnergy(MutationsEnergyListBestPerPosition, EnergyCutOff);

   /* count the lists and report that to the LOG*/
  int NumberOfMutationsCollected                       = CountMutationsProteinEnergy(MutationsEnergyList);
  int NumberOfMutationsCollectedBelowCutOff            = CountMutationsProteinEnergy(MutationsEnergyListBelowCutOff);
  int NumberOfMutationsCollectedPerPosition            = CountMutationsProteinEnergy(MutationsEnergyListBestPerPosition);
  int NumberOfMutationsCollectedPerPositionBelowCutOff = CountMutationsProteinEnergy(MutationsEnergyListBestPerPositionBelowCutOff);
  printf("\n");
  printf("The total number of mutations collected from FoldX is ................................................... %5d\n", NumberOfMutationsCollected);
  printf("Of those the following number of mutations passes the energy cut off..................................... %5d\n" ,NumberOfMutationsCollectedBelowCutOff );
  printf("Of those the following number of mutations was the best at that position................................. %5d\n" ,NumberOfMutationsCollectedPerPosition );
  printf("Of those the following number of mutations was the best at that position and passed the energy cut off... %5d\n" ,NumberOfMutationsCollectedPerPositionBelowCutOff );
  printf("The used energy cut off was ............................................................................. %9.3f kJ mol-1\n", EnergyCutOff);
  /* write them to files */
  WriteMutationsEnergyToFile(MutationsEnergyList, "MutationsEnergies_CompleteList.tab");
  WriteMutationsEnergyToFile(MutationsEnergyListBelowCutOff, "MutationsEnergies_BelowCutOff.tab");
  WriteMutationsEnergyToFile(MutationsEnergyListBestPerPosition, "MutationsEnergies_BestPerPosition.tab");
  WriteMutationsEnergyToFile(MutationsEnergyListBestPerPositionBelowCutOff, "MutationsEnergies_BestPerPositionBelowCutOff.tab");
  
  /*Sorting of files is carried out by linux */
  system(" cat MutationsEnergies_CompleteList.tab | sed 's/Below are/Below -1000000 are/g' | sort -n -k2 | sed 's/Below -1000000 are/Below are/g'>MutationsEnergies_CompleteList_SortedByEnergy.tab");
  system(" cat MutationsEnergies_BelowCutOff.tab | sed 's/Below are/Below -1000000 are/g' | sort -n -k2 | sed 's/Below -1000000 are/Below are/g'>MutationsEnergies_BelowCutOff_SortedByEnergy.tab");
  system(" cat MutationsEnergies_BestPerPosition.tab | sed 's/Below are/Below -1000000 are/g' | sort -n -k2 | sed 's/Below -1000000 are/Below are/g'>MutationsEnergies_BestPerPosition_SortedByEnergy.tab");
  system(" cat MutationsEnergies_BestPerPositionBelowCutOff.tab | sed 's/Below are/Below -1000000 are/g' | sort -n -k2 | sed 's/Below -1000000 are/Below are/g'>MutationsEnergies_BestPerPositionBelowCutOff_SortedByEnergy.tab");
  

  /*  give the file names explicitly in the LOG file */
  printf("\nThe results were written to files called:\n - MutationsEnergies_CompleteList.tab\n - MutationsEnergies_BelowCutOff.tab\n - MutationsEnergies_BestPerPosition.tab\n - MutationsEnergies_BestPerPositionBelowCutOff.tab\n");
  printf("\nThe same results sorted by best energy were written to files called:\n - MutationsEnergies_CompleteList_SortedByEnergy.tab\n - MutationsEnergies_BelowCutOff_SortedByEnergy.tab\n - MutationsEnergies_BestPerPosition_SortedByEnergy.tab\n - MutationsEnergies_BestPerPositionBelowCutOff_SortedByEnergy.tab\n");
  
 
  /* check if agrees with sequence list according to test list, if not give an error message after cleaning up the memory */
  MutationsEnergyList                           =  RemoveMutationProteinEnergy(MutationsEnergyList);
  MutationsEnergyListBelowCutOff                =  RemoveMutationProteinEnergy(MutationsEnergyListBelowCutOff);
  MutationsEnergyListBestPerPosition            =  RemoveMutationProteinEnergy(MutationsEnergyListBestPerPosition);
  MutationsEnergyListBestPerPositionBelowCutOff =  RemoveMutationProteinEnergy(MutationsEnergyListBestPerPositionBelowCutOff);
  }
  

  /* prevent leaving memory occupied until the computer is restarted*/
  MutatedProteinList = RemoveMutationProtein(MutatedProteinList); 
  RawSequence = RemoveSequenceProtein(RawSequence);
  if (NumberOfStrands != 1)
    PurgedSequence = RemoveSequenceProtein(PurgedSequence); 
    
  printf("\nFinished succesfully stage %s, exiting....\n\n", WhichStageAreWeIn);
  exit(EXIT_SUCCESS);
}
