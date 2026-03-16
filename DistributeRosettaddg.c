/***********************************************************************************************************************************************************
*                                                                                                                                                          *
*                                                                                                                                                          *
*                DistributeRosettaddg.c                                                                                                                    *
*                                                                                                                                                          *
*                                                                                                                                                          *
*  Author: Hein Wijma, University of Groningen, January 2011, h.j.wijma@rug.nl                                                                             *
*                                                                                                                                                          *
*  Can be used to automate the calculation by Rosetta of the change of folding energy for a large number of mutations                                      *
*                                                                                                                                                          *
*                                                                                                                                                          *
*   Rosetta itself can be obtained from www.rosettacommons.org                                                                                             *
*                                                                                                                                                          *
*                                                                                                                                                          *
*                                                                                                                                                          *
*                                                                                                                                                          *
***********************************************************************************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "DistributeRosettaddg.h"

int main(int argc, char *argv[])
{
  /* first some important declarations */
  char NamePDBFile[200], NameTableFile[200];                 /* names of the PDB file to move to the separate directories and of the table file with mutations to open */
  char strand1[4], strand2[4], strand3[4], strand4[4], strand5[4], strand6[4], strand7[4], strand8[4];      
                                                             /* program for maximally an octameric protein */
  int StartStrand1, StartStrand2, StartStrand3, StartStrand4, StartStrand5, StartStrand6, StartStrand7, StartStrand8;
                                                             /*  "       "    "       "    "         "     */
  int NumberOfStrands, NumberOfMutationsPerDirectory;        /* obvious */
  char WhichStageAreWeIn[20];                                /* check if it is Phase1 or Phase2 */
  float EnergyCutOff = 0;                                    /* in Phase2, make also a list of all the solutions with an energy better than this in kJ mol -1 */ 
  char FLAGFile[200], RosettaApplication[500];               /* If Phase1, need to know which FLAG file to use and which RosettaApplication */

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
  
  /* first check, number of arguments should be at least 8*/
  if (argc < 9)
     error_message("Too few arguments, see above");
  
  
  /* ==========================================assign all the char in the command line to what they should be ================================================*/
  /* =========================================================================================================================================================*/
  strcpy (WhichStageAreWeIn, argv[1]);
  printf("The name of the stage we are in as read in is %s\n", WhichStageAreWeIn);
  
  strcpy(NameTableFile, argv[2]);
  printf("The name of the table file is %s\n", NameTableFile);
 
  NumberOfStrands = (int) atof(argv[3]);
  printf("The number of strands read in is %d\n", NumberOfStrands);
  
  if (NumberOfStrands == 1)
  {
    strcpy(strand1,argv[4]);
    StartStrand1 = (int) atof(argv[5]);
    printf("The letter of this strand is %s and it starts at residue %d \n", strand1, StartStrand1);
  }  
  
  if (NumberOfStrands == 2)
  {
    strcpy(strand1,argv[4]);
    StartStrand1 = (int) atof(argv[5]);
    strcpy(strand2,argv[6]);
    StartStrand2 = (int) atof(argv[7]);
    printf("The letters of these strand are %s, %s, which start at residue %d, %d\n", strand1, strand2, StartStrand1, StartStrand2);
  }  
  
  if (NumberOfStrands == 3)
  {
    strcpy(strand1,argv[4]);
    StartStrand1 = (int) atof(argv[5]);
    strcpy(strand2,argv[6]);
    StartStrand2 = (int) atof(argv[7]);
    strcpy(strand3,argv[8]);
    StartStrand3 = (int) atof(argv[9]);
    printf("The letters of these strand are %s, %s, %s which start at residue %d, %d, %d\n", strand1, strand2, strand3, StartStrand1, StartStrand2, StartStrand3);
  }  
  
  if (NumberOfStrands == 4)
  {
    strcpy(strand1,argv[4]);
    StartStrand1 = (int) atof(argv[5]);
    strcpy(strand2,argv[6]);
    StartStrand2 = (int) atof(argv[7]);
    strcpy(strand3,argv[8]);
    StartStrand3 = (int) atof(argv[9]);
    strcpy(strand4,argv[10]);
    StartStrand4 = (int) atof(argv[11]);
    printf("The letters of these strand are %s, %s, %s, %s which start at residue %d, %d, %d, %d\n", strand1, strand2, strand3, strand4, StartStrand1, StartStrand2, StartStrand3, StartStrand4);
  }  
    
  if (NumberOfStrands == 5)
  {
    strcpy(strand1,argv[4]);
    StartStrand1 = (int) atof(argv[5]);
    strcpy(strand2,argv[6]);
    StartStrand2 = (int) atof(argv[7]);
    strcpy(strand3,argv[8]);
    StartStrand3 = (int) atof(argv[9]);
    strcpy(strand4,argv[10]);
    StartStrand4 = (int) atof(argv[11]);
    strcpy(strand5,argv[12]);
    StartStrand5 = (int) atof(argv[13]);
    printf("The letters of these strand are %s, %s, %s, %s, %s which start at residue %d, %d, %d, %d, %d\n", strand1, strand2, strand3, strand4, strand5,StartStrand1, StartStrand2, StartStrand3, StartStrand4, StartStrand5);
  }  

  if (NumberOfStrands == 6)
  {
    strcpy(strand1,argv[4]);
    StartStrand1 = (int) atof(argv[5]);
    strcpy(strand2,argv[6]);
    StartStrand2 = (int) atof(argv[7]);
    strcpy(strand3,argv[8]);
    StartStrand3 = (int) atof(argv[9]);
    strcpy(strand4,argv[10]);
    StartStrand4 = (int) atof(argv[11]);
    strcpy(strand5,argv[12]);
    StartStrand5 = (int) atof(argv[13]);
    strcpy(strand6,argv[14]);
    StartStrand6 = (int) atof(argv[15]);
    printf("The letters of these strand are %s, %s, %s, %s, %s, %s which start at residue %d, %d, %d, %d, %d, %d\n", strand1, strand2, strand3, strand4, strand5,strand6, StartStrand1, StartStrand2, StartStrand3, StartStrand4, StartStrand5, StartStrand6);
  }

  if (NumberOfStrands == 7)
  {
    strcpy(strand1,argv[4]);
    StartStrand1 = (int) atof(argv[5]);
    strcpy(strand2,argv[6]);
    StartStrand2 = (int) atof(argv[7]);
    strcpy(strand3,argv[8]);
    StartStrand3 = (int) atof(argv[9]);
    strcpy(strand4,argv[10]);
    StartStrand4 = (int) atof(argv[11]);
    strcpy(strand5,argv[12]);
    StartStrand5 = (int) atof(argv[13]);
    strcpy(strand6,argv[14]);
    StartStrand6 = (int) atof(argv[15]);
    strcpy(strand7,argv[16]);
    StartStrand7 = (int) atof(argv[17]);
    printf("The letters of these strand are %s, %s, %s, %s, %s, %s, %s which start at residue %d, %d, %d, %d, %d, %d, %d\n", strand1, strand2, strand3, strand4, strand5,strand6,strand7, StartStrand1, StartStrand2, StartStrand3, StartStrand4, StartStrand5, StartStrand6, StartStrand7);
  }

  if (NumberOfStrands == 8)
  {
    strcpy(strand1,argv[4]);
    StartStrand1 = (int) atof(argv[5]);
    strcpy(strand2,argv[6]); 
    StartStrand2 = (int) atof(argv[7]);
    strcpy(strand3,argv[8]);
    StartStrand3 = (int) atof(argv[9]);
    strcpy(strand4,argv[10]);
    StartStrand4 = (int) atof(argv[11]);
    strcpy(strand5,argv[12]);
    StartStrand5 = (int) atof(argv[13]);
    strcpy(strand6,argv[14]);
    StartStrand6 = (int) atof(argv[15]);
    strcpy(strand7,argv[16]);
    StartStrand7 = (int) atof(argv[17]);
    strcpy(strand8,argv[18]);
    StartStrand8 = (int) atof(argv[19]);
    printf("The letters of these strand are %s, %s, %s, %s, %s, %s, %s, %s which start at residue %d, %d, %d, %d, %d, %d, %d, %d\n", strand1, strand2, strand3, strand4, strand5,strand6,strand7,strand8, StartStrand1, StartStrand2, StartStrand3, StartStrand4, StartStrand5, StartStrand6, StartStrand7, StartStrand8);
  }

  strcpy(NamePDBFile, argv[4 + (2*NumberOfStrands)]);    
  printf("The name of the PDB file as read from the command line is %s\n", NamePDBFile);

  NumberOfMutationsPerDirectory = (int) atof(argv[5 + (2*NumberOfStrands)]);
  printf("The number of mutations per directory is %d\n", NumberOfMutationsPerDirectory);
  
  if (strcmp(WhichStageAreWeIn,"Phase1")== 0)
  { 
    strcpy(FLAGFile, argv[6 + (2*NumberOfStrands)]);
    printf("The desired Flag file is %s\n", FLAGFile);
  }   
 
  if (strcmp(WhichStageAreWeIn,"Phase1")== 0)
  { 
    strcpy(RosettaApplication, argv[7 + (2*NumberOfStrands)]);
    printf("The selected RosettaApplication is %s\n", RosettaApplication);
  }   
 
  if (strcmp(WhichStageAreWeIn,"Phase2")== 0)
  { 
    EnergyCutOff = (float) atof(argv[6 + (2*NumberOfStrands)]);
    printf("The energy cutoff is %.6f\n", EnergyCutOff);
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


  /* ===========================  make the list of mutations and write them to a series of files that can be used by Rosetta =====================================================*/
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
    /* calculations below are not done as a sub function since otherwise strand1, strand2, strand3, and strand 4 needed to be passed in some dynamic way that takes too much time for me to program*/
    
    /* make the subdirectories and write the mutations lists to them */
    MakeSubdirectories(NumberOfSubdirectories);
    
    char inbetweenstring[1000];
  
    /*write the sequence mutations to every one of them */
    FILE *mutationfile;
    FILE *readablemutationfile;
    FILE *todolist;
    char CurrentDirectoryName[200];
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
      
      strcat(MutationFileName,"RosettaFormatMutations.mut"); /*Rosetta does accept this name as the input list for mutations, for example mutationlist is not accepted */
      strcat(ReadableMutationFileName,"List_Mutations_readable.txt");
      
      /* fill first the readable mutation file mutation files */
      mutationfile = fopen(MutationFileName, "w");
      readablemutationfile = fopen(ReadableMutationFileName, "w");
      /*first print the total number of mutations, simple if not the last directory, otherwise calculate number of mutations in last directory  */
      if (h < (NumberOfSubdirectories - 1) )
        fprintf(mutationfile,"total %d\n", (NumberOfMutationsPerDirectory*NumberOfStrands));
      else
        fprintf(mutationfile,"total %d\n", ( (NumberOfMutations - (NumberOfMutationsPerDirectory*(NumberOfSubdirectories - 1) ) ) *NumberOfStrands)); 
      /* now really fill the mutation files */
      for (i = 0; i < NumberOfMutationsPerDirectory && p != NULL; i++, p = p->next)
      {  
        fprintf(mutationfile,"%d\n", NumberOfStrands);
        if (NumberOfStrands == 1)
        { 
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
        }
        if (NumberOfStrands == 2)
        { 
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand2 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand2 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
        }
        if (NumberOfStrands == 3)
        { 
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand2 - StartStrand1 - StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand3 - StartStrand1 - StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand2 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber +  StartStrand3 - StartStrand1 - StartStrand1 + 1),p->NewAminoAcid);
        }
        if (NumberOfStrands == 4)
        { 
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand2 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand3 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand4 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand2 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand3  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand4  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
        }
        if (NumberOfStrands == 5)
        { 
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand2 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand3 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand4 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand5 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand2 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand3  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand4  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand5  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
        }
        if (NumberOfStrands == 6)
        {
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand2 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand3 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand4 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand5 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand6 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand2 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand3  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand4  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand5  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand6  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
        }
        if (NumberOfStrands == 7)
        {
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand2 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand3 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand4 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand5 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand6 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand7 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);         
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand2 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand3  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand4  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand5  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand6  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand7  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
        }
        if (NumberOfStrands == 8)
        {
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand2 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand3 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand4 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand5 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand6 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand7 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(mutationfile,"%s %d %c\n",p->CurrentAminoAcid,(p->residueNumber + StartStrand8 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber -StartStrand1  + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand2 - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand3  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand4  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand5  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand6  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand7  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
          fprintf(readablemutationfile,"%s%d%c",p->CurrentAminoAcid,(p->residueNumber + StartStrand8  - StartStrand1 -StartStrand1 + 1),p->NewAminoAcid);
        }

        fprintf(readablemutationfile, " is %s %d %c\n", p->CurrentAminoAcid,(p->residueNumber),p->NewAminoAcid);
      }
      fclose(mutationfile);
      fclose(readablemutationfile);
      
      /* and also copy necessary files from the starting directory to the subdirectory */
      strcpy(inbetweenstring, "cp ");
      strcat(inbetweenstring, FLAGFile);
      strcat(inbetweenstring, " ");
      strcat(inbetweenstring, CurrentDirectoryName);
      system(inbetweenstring);
      
      strcpy(inbetweenstring, "cp ");
      strcat(inbetweenstring, NamePDBFile);
      strcat(inbetweenstring, " ");
      strcat(inbetweenstring, CurrentDirectoryName);
      system(inbetweenstring);
      
      /* add to the file todolist the following 3 lines */
      fprintf(todolist,"cd %s\n%s @%s -in:file:s %s -ddg::mut_file RosettaFormatMutations.mut >LOG&\ncd ..\n",CurrentDirectoryName,RosettaApplication,FLAGFile ,NamePDBFile); 
    }
    /* close the file todolist and make it executable */
    fclose(todolist);
    system("chmod +x todolist");
   
  /* here end a large loop for if this is step 1 */
  }
   
  /* ============================================ Here collect the results of the calculations, again not in separate function due to the variable number of strands======================= */
  if (strcmp(WhichStageAreWeIn,"Phase2")== 0)
    {
  printf("\n========= This software assume the user has waited for the Rosetta calculations to finish, e.g. type ps aux | grep Rosetta ============\n");
   
  /* Collect the data, analyse it and write it to the right files, write it in kJ mol instead of the kcal mol -1*/
  /* -----------------------------------------------------------------------------------------------------------*/
  
  
  /* for each of the subdirectories, read in the file Average_BuildModel_?????????.fxout, store the results */
  FILE *readablemutationfile;
  FILE *Rosettaenergyfile;
  char CurrentDirectoryName[200], inbetweenstring[400];
  char ReadableMutationFileName[400];
  char RosettaEnergyFileName[400];
  char c; /* check for end of file EOF */
  char TargetNameToLookFor[400];
  int h;

  /* initialize the structure of the mutation energy list*/
  struct MutationProteinEnergy *Before = NULL;
  struct MutationProteinEnergy * MutationsEnergyList; 
  
  printf("starting to read the energies from all the subdirectories\n");
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
    strcpy(RosettaEnergyFileName, CurrentDirectoryName);
    strcat(ReadableMutationFileName,"List_Mutations_readable.txt");
    strcat(RosettaEnergyFileName, "ddg_predictions.out");
    
    /* test whether the Rosetta energy file exists */
    printf("Currently looking for %s....\n", RosettaEnergyFileName);
    Rosettaenergyfile = fopen(RosettaEnergyFileName, "r");
    if (Rosettaenergyfile == '\0')
      error_message("can't find it or empty file "); 
    else
      printf("Found it, start analyzing now\n");
    
    /* idem for the file with the mutations */
    readablemutationfile = fopen(ReadableMutationFileName, "r");
    if  (readablemutationfile == '\0')
      error_message("can't find the list with mutations  "); 
    
    /* decide whether there is a useful part to read in or not, basiccaly the Rosettafiles start with a huge index, followed by now and then the mutations followed by the calculated total energy */
    /*the mutations need to come from the list of readable mutations in ListMutationsReadbale*/
    fscanf(readablemutationfile, "%s %s",TargetNameToLookFor, inbetweenstring);
    fscanf(Rosettaenergyfile," %s", inbetweenstring);
    do
    { 
      /*check whether OK to read */
      if (strcmp(TargetNameToLookFor,inbetweenstring)==0)
      {  /* get the data */
         MutationsEnergyList = malloc(sizeof(struct MutationProteinEnergy));
         fscanf(Rosettaenergyfile, " %f", &MutationsEnergyList->EnergyChange);                                                                                   /* here the energy is read in */
         fscanf(readablemutationfile, "%s %d %c",MutationsEnergyList->CurrentAminoAcid, &MutationsEnergyList->residueNumber, &MutationsEnergyList->NewAminoAcid);/* here the mutations */ 
         if (DEBUG)
           printf("Just read in %f after finding %s which has to be mutation %s %d %c \n",MutationsEnergyList->EnergyChange,TargetNameToLookFor, MutationsEnergyList->CurrentAminoAcid, MutationsEnergyList->residueNumber, MutationsEnergyList->NewAminoAcid);
          
         /* convert the energies to kJ mol -1 */
         MutationsEnergyList->EnergyChangeSD = (float) 4.1840000 * 0;                                  /*not calculated here */
         MutationsEnergyList->EnergyChange   = (float) 4.1840000 * MutationsEnergyList->EnergyChange;
         
         /* prepare for the next part of the list */
         MutationsEnergyList->next  =  Before;   
         Before                     = MutationsEnergyList;      
        
         /* prepare for next part of the cycle*/
         fscanf(readablemutationfile, "%s %s",TargetNameToLookFor, inbetweenstring);
      }     
      fscanf(Rosettaenergyfile," %s", inbetweenstring);

    } while((c = fgetc(Rosettaenergyfile)) != EOF);
    
    fclose(Rosettaenergyfile);
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
  printf("The total number of mutations collected from Rosetta is ................................................. %5d\n", NumberOfMutationsCollected);
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
