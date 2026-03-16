/***********************************************************************************************************************************************************
*                                                                                                                                                          *
*                                                                                                                                                          *
*                DistributeRosettaddg.h                                                                                                                    *
*                                                                                                                                                          *
*                                                                                                                                                          *
*  Author: Hein Wijma, University of Groningen, January 2011, h.j.wijma@rug.nl                                                                             *
*                                                                                                                                                          *
*                                                                                                                                                          *
***********************************************************************************************************************************************************/

/* if set to 1, much more information is provided which makes it easier to find out what is is going (wr)on(g) */
#define DEBUG 0


#define TRUE 1
#define FALSE 0


/*==============================================================The structure definitions================================================================*/
/*=======================================================================================================================================================*/

/*This is the structure definition for the sequences */
struct SequenceProtein
{
  char CurrentAminoAcid[4]; 
  char Chain[4];
  int  residueNumber;
  struct SequenceProtein *next;
};


/*This is the structure definition for the mutations */
struct MutationProtein
{
  char CurrentAminoAcid[4]; 
  char Chain[4];
  int  residueNumber;
  char NewAminoAcid;
  struct MutationProtein *next;
};


/*This is the structure definition for the mutations with energy and standard deviation of energy */
struct MutationProteinEnergy
{
  int IndexInDirectory;
  char CurrentAminoAcid[4]; 
  int  residueNumber;
  char NewAminoAcid;
  float EnergyChange;
  float EnergyChangeSD;
  struct MutationProteinEnergy *next;
};


/*=====================================The declarations of the functions ===================================================================================*/
/*==========================================================================================================================================================*/

/* This can print an error message */                                          
void error_message(char *Message);

/* this prints who made the program, etcetera */                                                           
void short_info(); 

/* this prints how to use this program */                                          
void example_use();

/* this reads the table file from yasara in a format M A 1\n Q A 2\n N A 3\n etcetera */   
struct SequenceProtein * LoadYasaraTable(char *NameTableFile);

/* this prints the sequence of a protein */
void PrintSequenceProtein(struct SequenceProtein *);

/* this prints the list of mutations of protein */
void PrintMutationsProtein(struct MutationProtein *);

/* this produces a list of the residues in strand Y, that are present also in strand X, this can be used to purge residues that are not present in both strands */ 
struct SequenceProtein * GetPurgedListResidues(struct SequenceProtein *, char *NameStrandX);

/* this produces a list of all the mutations in case of casting with two or maximally three residues at the same time */
struct MutationProtein * CastingMutagenesis(struct SequenceProtein *);

/* This writes a mutation with energies to a file */
void WriteMutationsEnergyToFile(struct MutationProteinEnergy *, char *NameFile);

/* This removes all mutations with a value above the cuttoff value.  */
struct MutationProteinEnergy * CuttOffMutationsEnergy(struct MutationProteinEnergy *, float);

/* This counts the number of recorded mutations with energy change  */
int CountMutationsProteinEnergy(struct MutationProteinEnergy *);

/* count number of mutations */
int CountNumberMutations(struct MutationProtein *);

/* This removes the stored list of sequences out of memory */
struct SequenceProtein * RemoveSequenceProtein(struct SequenceProtein *);

/* This removes the stored list of mutations out of memory */
struct MutationProtein * RemoveMutationProtein( struct MutationProtein *); 

/* This removes the stored list of mutations with energy out of memory */
struct MutationProteinEnergy * RemoveMutationProteinEnergy( struct MutationProteinEnergy *); 

/* This makes the subdirectories to do the Rosetta calculations in */
void MakeSubdirectories( int NumberOfSubdirectories);

/* these standard functions (itoa/reverse) are absent from some compilers  */
void treverse(char s[]);
void titoa(int n, char s[]);

/*============================================================the actual functions=========================================================================*/
/*========================================================================================================================================================*/ 

/* 
the following should
- LOAD a YASARA table file
- skip the first line
- store the result in a protein sequence structure
- give an error message if the table file does not exist
*/
struct SequenceProtein * LoadYasaraTable(char *NameTableFile)
{
  /* initalize the first bit of protein sequence */
  struct SequenceProtein *Before = NULL;
  
  FILE *OpenedTableFile;
  
  OpenedTableFile = fopen(NameTableFile,"r");
 
  if (OpenedTableFile == '\0') /* check that the file really exists*/
    {  
      char CompositeMessage[200] = "\0";
      strcat(CompositeMessage,"Table file ");
      strcat(CompositeMessage,NameTableFile);
      strcat(CompositeMessage, " does not exist"); 
      error_message(CompositeMessage);
    }
   
  /* put the data from the file into the SequenceProtein Structure ,assume all the lines are useful to read unless at says END */

  /* skip the first line */
  char FirstLineDiscarded[2000];
  fgets(FirstLineDiscarded,1999,OpenedTableFile);
  if (DEBUG)
    printf("This line was discarded while reading the table file: %s", FirstLineDiscarded);
  
  char AnAminoAcid[4];
  fscanf(OpenedTableFile, "%3s",AnAminoAcid);
  while (strcmp(AnAminoAcid,"END") != 0)
  {
    /* add to the list */
    struct SequenceProtein * NewAminoAcid       = malloc(sizeof(struct SequenceProtein));
    
    /* fill the following by simple assingment */
    strcpy(NewAminoAcid->CurrentAminoAcid, AnAminoAcid);
    
    /* add the other data straight from fscanf */
    fscanf(OpenedTableFile,"%3s %d",NewAminoAcid->Chain, &NewAminoAcid->residueNumber);
   
    /* if debug, report what is scanned */
    if (DEBUG)
      printf("JustScanned %s %s %d\n", NewAminoAcid->CurrentAminoAcid, NewAminoAcid->Chain,NewAminoAcid->residueNumber );
    /* now connect it to the previous atom and create a new Atom in the list */
    NewAminoAcid->next = Before;
    Before        = NewAminoAcid;
    
    /* scan the next line to see if it is not END */
    fscanf(OpenedTableFile, "%3s",AnAminoAcid);
  }

  fclose(OpenedTableFile);

     
  return Before;  
}


/* this prints the sequence of a protein in inverse order, just for debugging*/
void PrintSequenceProtein(struct SequenceProtein *NameSequenceProtein)
{
  if (NameSequenceProtein == NULL)
    printf("Called sequence protein which does not (longer) exist\n");
  
  struct SequenceProtein *p;
  
  for (p = NameSequenceProtein; p != NULL; p = p->next)
    printf("Stored aminoacid %s on strand %s of residue number %d\n" ,p->CurrentAminoAcid ,p->Chain ,p->residueNumber);  
}


/* this prints the list of mutations of protein which happen to be in the correct order */
void PrintMutationsProtein(struct MutationProtein *NameMutationsList)
{
  if (NameMutationsList == NULL)
    printf("Called mutationslist which is empty or does not longer exist\n");
  
  struct MutationProtein *p;
  
  for (p = NameMutationsList; p != NULL; p = p->next)
    printf("Planned the mutation of aminoacid %s on strand %s of residue number %d to residue %c\n" ,p->CurrentAminoAcid ,p->Chain ,p->residueNumber, p->NewAminoAcid);  
}


/* this prints the list of mutations of protein with predicted energy which happen to be in the inverse order */
void PrintMutationsProteinEnergy(struct MutationProteinEnergy *NameMutationsEnergyList)
{
  if (NameMutationsEnergyList == NULL)
    printf("Called mutationslist which is empty or does not longer exist\n");
  
  struct MutationProteinEnergy *p;
  
  for (p = NameMutationsEnergyList; p != NULL; p = p->next)
    printf("The mutation of aminoacid %s on of residue number %d to residue %c with local index %d is predicted to have a change of %.3f +/- %.3f kJ mol -1\n" ,p->CurrentAminoAcid,p->residueNumber, p->NewAminoAcid, p->IndexInDirectory, p->EnergyChange, p->EnergyChangeSD);  
}


void WriteMutationsEnergyToFile(struct MutationProteinEnergy * MutationProteinEnergyList, char *NameFile)
{
  /* since the list is currently in inverse order, make a copy that puts it in the right order, before doing anything else*/
  struct MutationProteinEnergy *p;
  struct MutationProteinEnergy *q;
  
  struct MutationProteinEnergy *Before =NULL;
  
  for (q = MutationProteinEnergyList; q != NULL; q = q->next)
  {   /* claim memory */
      p = malloc(sizeof(struct MutationProteinEnergy));
     
      /* assigning values */
      strcpy(p->CurrentAminoAcid,q-> CurrentAminoAcid);
      p->NewAminoAcid = q-> NewAminoAcid;
      p-> IndexInDirectory= q-> IndexInDirectory;
      p-> residueNumber= q-> residueNumber;
      p-> EnergyChange= q->EnergyChange ;
      p-> EnergyChangeSD= q->EnergyChangeSD ;
     
      /* and prepare for the next residue */
      p-> next = Before;
      Before   = p;
  }  
  
  /* start making the file */
  FILE *MutationEnergyFile = fopen(NameFile, "w");
  fprintf(MutationEnergyFile, "Below are the mutations and the change in stability and SD from 5 calculations, all in kJ mol -1\n");
  for (q = p; q != NULL; q = q->next)
    fprintf(MutationEnergyFile,"%s%d%c %8.3f %8.3f\n" ,q->CurrentAminoAcid, q->residueNumber, q->NewAminoAcid, q->EnergyChange, q->EnergyChangeSD);
  fclose(MutationEnergyFile);

  /*remove from memory the intermediate table */
  p = RemoveMutationProteinEnergy(p);
}


/* this produces a list of the residues in strands other than X, that are present also in strand X, this can be used to purge residues that are not present in both strands */ 
struct SequenceProtein * GetPurgedListResidues(struct SequenceProtein *NameSequenceProtein, char *NameStrandX)
{
  /* first check it still exists */
  if (NameSequenceProtein == NULL)
    printf("Called sequence protein which does not (longer) exist at the stage of ==== Purging ======\n");
  
  /* now make a copy of all residues that are present in strand X, start by initalization */
  struct SequenceProtein *ReferenceSet;
  struct SequenceProtein *Before  = NULL;
  struct SequenceProtein *p;
  for (p = NameSequenceProtein; p != NULL; p = p->next)
  {
    if (strcmp(NameStrandX, p->Chain) == 0)
    {
      /* claim memory */
      ReferenceSet = malloc(sizeof(struct SequenceProtein));
 
      /* assigning values */
      strcpy(ReferenceSet->CurrentAminoAcid,p->CurrentAminoAcid);
      strcpy(ReferenceSet->Chain,p->Chain);
      ReferenceSet->residueNumber = p->residueNumber;
 
      /* and prepare for the next residue */
      ReferenceSet-> next = Before;
      Before              = ReferenceSet;
    }
  }

  /* now make a copy of all residues that are present in the other strands and have the same amino acid type and number */
  struct SequenceProtein *PurgedSet;
  Before  = NULL;
  struct SequenceProtein *q;  
  int FoundOnReferenceStrand = FALSE;
  for (p = NameSequenceProtein; p!= NULL;p= p->next)
  {
    if (strcmp(NameStrandX, p->Chain) != 0)
    {
      /* we know now that this residue is not the one on the reference strand, now check if the same residue is present on the reference strand */
      for (q = ReferenceSet; q != NULL; q = q->next)
      {
        if (p->residueNumber == q->residueNumber)
        {
        if (strcmp(p->CurrentAminoAcid,q->CurrentAminoAcid)==0) 
          FoundOnReferenceStrand = TRUE;
        }
      }
      
      if (FoundOnReferenceStrand == TRUE)
      {
        /* claim memory */
        PurgedSet = malloc(sizeof(struct SequenceProtein));
 
        /* assigning values */
        strcpy(PurgedSet->CurrentAminoAcid,p->CurrentAminoAcid);
        strcpy(PurgedSet->Chain,p->Chain);
        PurgedSet->residueNumber = p->residueNumber;
 
        /* and prepare for the next residue */
        PurgedSet-> next   = Before;
        Before             =PurgedSet;
      }
      /* reset for next residue*/
      FoundOnReferenceStrand = FALSE;
    }
  }
  
  /* if DEBUG, report the current result */
  if (DEBUG)
  {
    printf("Reached the Stage of ==== Purging ==== resulting in \n");
    PrintSequenceProtein(PurgedSet);
  }

  /* put them in the normal order again, i.e. the inverse sequence*/
  struct SequenceProtein *InversePurgedSet;
  Before  = NULL;
  for (p = PurgedSet; p!= NULL;p= p->next)
  {
    /* claim memory */
    InversePurgedSet = malloc(sizeof(struct SequenceProtein));
 
    /* assigning values */
    strcpy(InversePurgedSet->CurrentAminoAcid,p->CurrentAminoAcid);
    strcpy(InversePurgedSet->Chain,p->Chain);
    InversePurgedSet->residueNumber = p->residueNumber;
 
    /* and prepare for the next residue */
    InversePurgedSet-> next   = Before;
    Before             =InversePurgedSet;
  }
  
  /* remove the intermediate table from memory */
  ReferenceSet = RemoveSequenceProtein(ReferenceSet);
  PurgedSet    = RemoveSequenceProtein(PurgedSet); 
  /* return a pointer to the memory address of the new protein sequence */
  return(InversePurgedSet);
}


/* this produces a list of all the mutations in the case of saturation scanning mutagenesis. It does 19 aminoacids, no Cysteine since no code for disulfide bridges and multimerization down the line */
struct MutationProtein * SaturationScanning(struct SequenceProtein *NameSequenceProtein)
{
  /* check if input is OK */
  if (NameSequenceProtein == NULL)
    error_message("Called sequence protein which does not (longer) exist at the stage of making mutations");
  
  /* initalize the mutant list */
  struct MutationProtein *ListMutations;
  struct MutationProtein *Before = NULL;
  struct SequenceProtein *p;
  int i;
   
  /* simply go to every residue of the sequence and make all the mutations, i.e. 19, not 20 except cystein, the resulting list will be in normal order again */
  char ListUsedAminoAcids[30];

  /*strcpy(ListUsedAminoAcids, "ADEFGHIKLMNPQRSTVWY");*/
  strcpy(ListUsedAminoAcids, "YWVTSRQPNMLKIHGFEDA");
  char OneLetter;
  for (p = NameSequenceProtein; p != NULL; p = p->next)
  {
    for (i = 0; i < 19; i++)
    {
      ListMutations = malloc(sizeof(struct MutationProtein));
      strcpy(ListMutations->CurrentAminoAcid,p->CurrentAminoAcid);
      strcpy(ListMutations->Chain,"X");
      ListMutations->residueNumber = p->residueNumber;
      OneLetter = ListUsedAminoAcids[i]; 
      ListMutations->NewAminoAcid= OneLetter;
      
      ListMutations->next   = Before;
      Before                = ListMutations;
    }
  }
  return(ListMutations); 
}


/* This removes all mutations with a value above the cuttoff value.  */
struct MutationProteinEnergy *CuttOffMutationsEnergy(struct MutationProteinEnergy *NameListMutationsEnergy, float CuttOffValue)
{
  /* first check it still exists */
  if (NameListMutationsEnergy == NULL)
    printf("Calledlist mutations with energy that does not (longer) exist at the stage of ==== eliminating values above cuttoff ======\n");
  
  /* now make a copy of all list members that have an energy value lower than the cut off value */
  struct MutationProteinEnergy *BelowCutOff;
  struct MutationProteinEnergy *Before  = NULL;
  struct MutationProteinEnergy *p;

  for (p = NameListMutationsEnergy; p != NULL; p = p->next)
  {
    if (p->EnergyChange <= CuttOffValue)
    {
      /* claim memory */
      BelowCutOff= malloc(sizeof(struct MutationProteinEnergy));
 
      /* assigning values */
      BelowCutOff->IndexInDirectory = p->IndexInDirectory;
      strcpy(BelowCutOff->CurrentAminoAcid,p->CurrentAminoAcid);
      BelowCutOff->residueNumber    = p->residueNumber; 
      BelowCutOff->NewAminoAcid     = p->NewAminoAcid;  
      BelowCutOff->EnergyChange     = p->EnergyChange;  
      BelowCutOff->EnergyChangeSD   = p->EnergyChangeSD;
      
      /* and prepare for the next residue */
      BelowCutOff-> next = Before;
      Before             = BelowCutOff;
    }
  }

  /* put them in the normal order again, i.e. the inverse sequence*/
  struct  MutationProteinEnergy *InversePurgedSet;
  Before  = NULL;
  for (p = BelowCutOff; p!= NULL;p= p->next)
  {
    /* claim memory */
    InversePurgedSet= malloc(sizeof(struct MutationProteinEnergy));
 
    /* assigning values */
    InversePurgedSet->IndexInDirectory = p->IndexInDirectory;
    strcpy(InversePurgedSet->CurrentAminoAcid,p->CurrentAminoAcid);
    InversePurgedSet->residueNumber    = p->residueNumber; 
    InversePurgedSet->NewAminoAcid     = p->NewAminoAcid;  
    InversePurgedSet->EnergyChange     = p->EnergyChange;  
    InversePurgedSet->EnergyChangeSD   = p->EnergyChangeSD;
    
    /* and prepare for the next residue */
    InversePurgedSet-> next = Before;
    Before                  = InversePurgedSet;
  }
  
  /* remove the intermediate table from memory */
  BelowCutOff = RemoveMutationProteinEnergy(BelowCutOff);
  
  /* return a pointer to the memory address of the new protein sequence */
  return(InversePurgedSet);
}


/* This Keeps the best solution at every position */
struct MutationProteinEnergy *KeepBestMutationsEnergy(struct MutationProteinEnergy *NameListMutationsEnergy)
{
  /* first check it still exists */
  if (NameListMutationsEnergy == NULL)
    printf("Calledlist mutations with energy that does not (longer) exist at the stage of ==== determining the best mutation at every position ======\n");
  
  /* now make a new list with only the best energy members for all positions */
  struct MutationProteinEnergy *BestPerPosition;
  struct MutationProteinEnergy *Before  = NULL;
  struct MutationProteinEnergy *p;
  int ResidueNumberPrevious = -999;
  
  /* the following loop works such that it makes a new member and copies the values if a new list member is involved, and replaces the values if the energy is better for a next residue*/
  for (p = NameListMutationsEnergy; p != NULL; p = p->next)
  {
    /* if residue number is new, add a memeber to the list */
    if (DEBUG)
      printf("CurrentResidueNumber is %d\n", ResidueNumberPrevious);
    if (ResidueNumberPrevious != p->residueNumber) 
    {
      if (DEBUG)
        printf("Found New residue number is %d, best energy is now %.3f\n",  p->residueNumber, p-> EnergyChange);

      /* claim memory */
      BestPerPosition= malloc(sizeof(struct MutationProteinEnergy));
 
      /* assigning values */
      BestPerPosition->IndexInDirectory = p->IndexInDirectory;
      strcpy(BestPerPosition->CurrentAminoAcid,p->CurrentAminoAcid);
      BestPerPosition->residueNumber    = p->residueNumber; 
      BestPerPosition->NewAminoAcid     = p->NewAminoAcid;  
      BestPerPosition->EnergyChange     = p->EnergyChange;  
      BestPerPosition->EnergyChangeSD   = p->EnergyChangeSD;
      
      /* and prepare for the next residue */
      BestPerPosition-> next = Before;
      Before                 = BestPerPosition;
      ResidueNumberPrevious  = BestPerPosition->residueNumber;  
    }
    
    /* if residue number has lower energy than current member in list, replace it */
    if (BestPerPosition->EnergyChange > p->EnergyChange) 
    {
      if (DEBUG)
        printf("Found member with better energy than previous, best energy is now %.3f instead of %.3f\n",  p-> EnergyChange, BestPerPosition->EnergyChange);

      /* assigning values */
      BestPerPosition->IndexInDirectory = p->IndexInDirectory;
      strcpy(BestPerPosition->CurrentAminoAcid,p->CurrentAminoAcid);
      BestPerPosition->residueNumber    = p->residueNumber; 
      BestPerPosition->NewAminoAcid     = p->NewAminoAcid;  
      BestPerPosition->EnergyChange     = p->EnergyChange;  
      BestPerPosition->EnergyChangeSD   = p->EnergyChangeSD;
    }
  }

  /* put them in the normal order again, i.e. the inverse sequence*/
  struct  MutationProteinEnergy *InversePurgedSet;
  Before  = NULL;
  for (p = BestPerPosition; p!= NULL;p= p->next)
  {
    /* claim memory */
    InversePurgedSet= malloc(sizeof(struct MutationProteinEnergy));
 
    /* assigning values */
    InversePurgedSet->IndexInDirectory = p->IndexInDirectory;
    strcpy(InversePurgedSet->CurrentAminoAcid,p->CurrentAminoAcid);
    InversePurgedSet->residueNumber    = p->residueNumber; 
    InversePurgedSet->NewAminoAcid     = p->NewAminoAcid;  
    InversePurgedSet->EnergyChange     = p->EnergyChange;  
    InversePurgedSet->EnergyChangeSD   = p->EnergyChangeSD;
    
    /* and prepare for the next residue */
    InversePurgedSet-> next = Before;
    Before                  = InversePurgedSet;
  }
  
  /* remove the intermediate table from memory */
  BestPerPosition = RemoveMutationProteinEnergy(BestPerPosition);
  
  /* return a pointer to the memory address of the new protein sequence */
  return(InversePurgedSet);
}


/* This counts the number of recorded mutations with energy change  */
int CountMutationsProteinEnergy(struct MutationProteinEnergy *NameListMutationsEnergy)
{
  /* check if input is OK */
  if (NameListMutationsEnergy == NULL)
    error_message("Called mutations list which does not (longer) exist at the stage of making mutations");
  
  int NumberOfMutations = 0;
  
  
  struct MutationProteinEnergy *p;
  for (p = NameListMutationsEnergy; p!=NULL; p = p->next)
  {
    NumberOfMutations++;  
  }
  
  return(NumberOfMutations);
}


/* count number of mutations */
int CountNumberMutations(struct MutationProtein *NameListMutations)
{
  /* check if input is OK */
  if (NameListMutations == NULL)
    error_message("Called mutations list which does not (longer) exist at the stage of making mutations");
  
  int NumberOfMutations = 0;
  
  
  struct MutationProtein *p;
  for (p = NameListMutations; p!=NULL; p = p->next)
  {
    NumberOfMutations++;  
  }
  
  return(NumberOfMutations);
}


/* This removes the stored list of sequences out of memory */
struct SequenceProtein * RemoveSequenceProtein(struct SequenceProtein *NameSequenceList)
{
  struct SequenceProtein * p;
  struct SequenceProtein * Previous;
  
  p = NameSequenceList;
  do 
  {
    Previous = p;
    p = p->next;
    free(Previous);  
  } while (p != NULL);
  
  return p;
}


/* This removes the stored list of mutations out of memory */
struct MutationProtein * RemoveMutationProtein( struct MutationProtein *NameMutationList)
{
  struct MutationProtein * p;
  struct MutationProtein * Previous;
  
  p = NameMutationList;
  do 
  {
    Previous = p;
    p = p->next;
    free(Previous);  
  } while (p != NULL);
  
  return p;
}


/* This removes the stored list of mutations with energy out of memory */
struct MutationProteinEnergy * RemoveMutationProteinEnergy( struct MutationProteinEnergy * NameMutationEnergyList)
{
  struct MutationProteinEnergy * p;
  struct MutationProteinEnergy * Previous;
  
  p = NameMutationEnergyList;
  do 
  {
    Previous = p;
    p = p->next;
    free(Previous);  
  } while (p != NULL);
  
  return p;
}



void error_message(char *Message)
{
   printf("-------- ERROR -------- %s -------- \n", Message);
   exit (EXIT_FAILURE);
}


/* function short info */
void short_info()
{
   printf("\nThis program is intended to distribute the calculation of a large number of ddg calculations for single mutations by Rosetta.\n");
}


/* function example use */
void example_use()
{
   printf("a command like:\n\n");
   printf("\tDistributeRosettaddg Phase1 DimerSelection.tab 2 A 4 B 149 DimerParsedForRosetta.pdb 100 FLAGrow3 /home/wijma/mini20101104/mini/bin/fix_bb_monomer_ddg.linuxgccrelease\n\n");
   printf("Should prepare mutations of DimerParsedForRosetta.pdb in the 2 subunits A and B of the residues in the File MyPreciousSelection.tab and\n");
   printf("distribute them over directories with each a 100 different mutations, in the pdb file the subunits start at residue 4 and 149. \n");
   printf("FLAGrow3 is the FLAG file to be used and /home/wijma/mini20101104/mini/bin/fix_bb_monomer_ddg.linuxgccrelease is the name and location of the ddg software\n\n");
      
   printf("\tDistributeRosettaddg Phase2 DimerSelection.tab 2 A 4 B 149  DimerParsedForRosetta.pdb 100 -5\n\n");
   printf("Should COLLECT the calculated mutations as above and\n");
   printf("AND make EIGHT lists (sorted/unsorted): a completelist, a separate list of the mutations that are improved by more than 5 kJ mol -1, a list with the best \nmutations for every position, (and a similar list with less than 10 kJ mol-1 damage).\n\n");
}


/* This makes the subdirectories to do the Rosetta calculations in */
void MakeSubdirectories( int NumberOfSubdirectories)
{
  int i;
  char SystemMessage[1000], inbetweenstring[10];
  
  for (i = 0; i<NumberOfSubdirectories ; i++)
  { 
    strcpy(SystemMessage, "mkdir Subdirectory");
    titoa((i+1), inbetweenstring);
    strcat(SystemMessage, inbetweenstring);
    system(SystemMessage);  
  }
}


/* these standard functions (itoa/reverse) are absent with some compilers */

/* reverse:  reverse string s in place */
void treverse(char s[])
 {
     int i, j;
     char c;
 
     for (i = 0, j = strlen(s)-1; i<j; i++, j--) {
         c = s[i];
         s[i] = s[j];
         s[j] = c;
     }
 }

/* itoa function */
void titoa(int n, char s[])
 {
     int i, sign;
 
     if ((sign = n) < 0)  /* record sign */
         n = -n;          /* make n positive */
     i = 0;
     do {       /* generate digits in reverse order */
         s[i++] = n % 10 + '0';   /* get next digit */
     } while ((n /= 10) > 0);     /* delete it */
     if (sign < 0)
         s[i++] = '-';
     s[i] = '\0';
     treverse(s);
 } 
 
