/**************************************************************************************************************************************************
*                                                                                                                                                 *
*                        OverviewDisulfides.h                                                                                                     *
*                                                                                                                                                 *
*   Author: Hein Wijma, University of Groningen, 2012-4-27, h.j.wijma@rug.nl                                                                      *
*                                                                                                                                                 *
*  Can be used to get an overview of the disulfide bonds created by DisulfideBridge.mcr                                                           *
*                                                                                                                                                 *
*                                                                                                                                                 *
*                                                                                                                                                 *
*                                                                                                                                                 *
*                                                                                                                                                 *
**************************************************************************************************************************************************/

/* if set to 1, much more information is provided which makes it easier to find out what is is going (wr)on(g) */
#define DEBUG 0


#define TRUE 1
#define FALSE 0


/*==============================================================The structure definitions================================================================*/
/*=======================================================================================================================================================*/

/*This is the structure definition for the disulfide bonds with geometric information, energies, and all other information with energy and standard deviation of energy */
struct DisulfideBond
{
  /*the following parameters are read in directly from the raw table file */
  float d_CA_CA;
  float d_CB_CB;
  float d_SG_SG;
  float adCA_CB_SG;
  float aaCA_CB_SG;
  float adCB_SG_SG;
  float aaCB_SG_SG;
  float hd_Chi1;
  float ha_Chi1;
  float hd_Chi2;
  float ha_Chi2;
  float hh_SS;
  float Energy_SS;
  
  /*the name of the original directory needs to be remembered for every subheading in the original table file */
  char OriginalDirectory[1000];   

  /*the name of the pdb file is again read in directly from the original table file*/
  char Name_Pdb_File[200];
  
  /*the following have to be extracted from the name of the pdb file */
  char Name_Pdb_File_noextension[196]; /*i.e. no .pdb file extension */
  char Name_Starting_Structure[196];
  char NameDonorStrand[4];
  char NameAcceptorStrand[4];
  int ResidueNumberDonorStrand;
  int ResidueNumberAcceptorStrand;
  char CurrentAminoAcidDonor[4];
  char CurrentAminoAcidAcceptor[4];
  char AbbreviatedGeometry[50];
  
  int InterChainDisulfide;
  
  /*next node*/
  struct DisulfideBond *next;
};




/*=====================================The declarations of the functions ===================================================================================*/
/*==========================================================================================================================================================*/

/* This can print an error message */                                          
void error_message(char *Message);

/* list original command line material */
void ListOriginalCommandLine(int argc, char *argv[]);

/* this prints who made the program, etcetera */                                                           
void short_info(); 

/* this prints how to use this program */                                          
void example_use();

/* This creates the unprocessed list of all disulfide bonds and its conformations */
void CreateRawDisulfideList();

/*this reads in all the disulfide bonds from CreateRawDisulfideList*/
struct DisulfideBond * LoadListFromFile(char *NameFile);


/* This removes the stored list of disulfide bonds out of memory */
struct DisulfideBond * RemoveDisulfidesFromMemory(struct DisulfideBond * NameDisulfideList);

/* count number of disulfide bridges */
int CountNumberDisulfideBridges(struct DisulfideBond *NameListDisulfides);

/* this prints the list of disulfide bridges */
void PrintDisulfideBridges(struct DisulfideBond *NameListDisulfides);

/* this prints the list of unique disulfide bridges, provided already in sequence */
void PrintDisulfideBridgesUnique(struct DisulfideBond *NameListDisulfides);

/*simple sort algorithm, bubble sort, to sort the list of disulfide bonds*/
struct DisulfideBond * SortDisulfideList(struct DisulfideBond *NameDisulfideList);


/*============================================================the actual functions=========================================================================*/
/*========================================================================================================================================================*/ 

/* eror message */
void error_message(char *Message)
{
   printf("-------- ERROR -------- %s -------- \n", Message);
   exit (EXIT_FAILURE);
}

/* list original command line material */
void ListOriginalCommandLine(int argc, char *argv[])
{
  printf("The raw command line input is\n");
  int i;
  for (i = 0; i < argc; i++)
  {
    printf("%s\n", argv[i]);
  }
}


/* function short info */
void short_info()
{
   printf("\nThis program is intended to make an overview of the the disulfide bridges designed by DisulfideBridge.mcr.\nThe *Conformations.tab files should be in subdirectories below the current directory.\n\n");
}


/* function example use */
void example_use()
{
   printf("a command like:\n\n");
   printf("\tOverviewDisulfides \n\n");
   printf("Should analyse the created disulfide bonds.\n");
   printf("Should then create\n\t-a few tab files that give an overview of all the different disulfide bonds\n\t-and a file to_copy_list that can be used to make copies of all the unique disulfide bonds\n\n");
}


/* This creates the unprocessed list of all disulfide bonds and its conformations */
void CreateRawDisulfideList()
{
   system("head -1000 */*_ConformationsDisulfideBonds.tab|sed 's/==> //g'| sed 's/ <==//g'>  RawDisulfideListtmp.tab");
   /* temporary for programming on another computer that lacks the subdirectories with information*/
}

/*struct SequenceProtein * LoadYasaraTable(char *NameTableFile)*/
struct DisulfideBond * LoadListFromFile(char *NameFile)
{
  /* some important declarations */
  char CurrentDirectory[1000], CurrentStartingStructure[1000];                                                            /* names of the current directory and of the starting structure */
  char CharforColons[1000], CharforDirectoryAndFileName[2000], CharforInfolineAbove[10000], CharforInfolineBelow[10000];    /* program assumes standard formatting of file, no exceptions allowed */
  char c, CharTemporyForDeciding[10000], CharNumbers[60];                                                                 /* when need to decide if something useful or not, or to convert something */
  int i,j;                                                                                                              /* integers for book keeping */
  /* open the file and check if it really exists */ 
  FILE *OpenedFile;
  OpenedFile = fopen(NameFile,"r");
  if (OpenedFile == '\0') /* check that the file really exists*/
    {  
      char CompositeMessage[200] = "\0";
      strcat(CompositeMessage,"File ");
      strcat(CompositeMessage,NameFile);
      strcat(CompositeMessage, " does not exist"); 
      error_message(CompositeMessage);
    }
   
  
  /* initalize the first disulfide bond*/
  struct DisulfideBond *Before = NULL;
  
  /* start reading the file*/ 
  do
  { 
    /* assume first some unimportant lines in standard format*/
    /*fgets(CharforColons,100,OpenedFile);
    if (DEBUG)
      printf("Assumed that %s was unimportant\n", CharforColons);*/
    fgets(CharforDirectoryAndFileName, 400, OpenedFile);
    if (DEBUG)
      printf("Assumed that %s was the name of the current directory and the table file\n", CharforDirectoryAndFileName);
    fgets(CharforColons,100,OpenedFile);
    if (DEBUG)
      printf("Assumed that %s was unimportant\n", CharforColons);
    fgets(CharforInfolineAbove,1000,OpenedFile);
    if (DEBUG)
      printf("Assumed that %s was unimportant being the line of info abvoe\n", CharforInfolineAbove);
 
    /* assign the name of the directory and the name of the Current starting structure*/
    for (i = 0; CharforDirectoryAndFileName[i] != '/'; i++)
    { 
      CurrentDirectory[i] = CharforDirectoryAndFileName[i];
      CurrentDirectory[i+1] = '\0';
    }
    if (DEBUG)
      printf("As directory name was identified %s\n", CurrentDirectory); 
    for (i++, j = 0;  CharforDirectoryAndFileName[i] != '\n'; j++,i++)
    {
      CurrentStartingStructure[j] = CharforDirectoryAndFileName[i];
      CurrentStartingStructure[j+1] = '\0';
    }
    /* and now remove the _ConformationsDisulfideBonds.tab (32 chars) */
    CurrentStartingStructure[j-32] = '\0';
    if (DEBUG)
      printf("As current starting structure was identified %s\n", CurrentStartingStructure); 
    
          
    /* now check whether there is something useful here, if so use it and try again, else skip the line entirely and continue with the next set of of lines */  
    do
    {
      fscanf(OpenedFile,"%1000s", CharTemporyForDeciding);
      if (DEBUG)
        printf("Read into CharTemporaryForDeciding, %s\n", CharTemporyForDeciding);

      if (strcmp(CharTemporyForDeciding, "d") != 0)
      {
         if (DEBUG)
           printf("reached that it makes a new disulfide bond in memory\n");
         /* add the new node to the list */
         struct DisulfideBond * NewDisulfideBond       = malloc(sizeof(struct DisulfideBond));
         
         if (DEBUG)
           printf("reached that it made a new disulfide bond in memory\n");
         
         /* convert the temporary char into a float. */
         NewDisulfideBond->d_CA_CA =atof(CharTemporyForDeciding);
         if (DEBUG)
           printf("reached A1\n");
         /*and get the known starting structure and directory also at their spot in memory */
         strcpy(NewDisulfideBond->OriginalDirectory, CurrentDirectory);
         if (DEBUG)
           printf("reached A2\n");
         strcpy(NewDisulfideBond->Name_Starting_Structure, CurrentStartingStructure);
          if (DEBUG)
           printf("reached B\n");
         /* add the other data straight from fscanf */
         fscanf(OpenedFile,"%f %f %f ",&NewDisulfideBond->d_CB_CB, &NewDisulfideBond->d_SG_SG,&NewDisulfideBond->adCA_CB_SG);              
         fscanf(OpenedFile,"%f %f %f ",&NewDisulfideBond->aaCA_CB_SG,&NewDisulfideBond->adCB_SG_SG,&NewDisulfideBond->aaCB_SG_SG);
         fscanf(OpenedFile,"%f %f %f ",&NewDisulfideBond->hd_Chi1,&NewDisulfideBond->ha_Chi1,&NewDisulfideBond->hd_Chi2);
         fscanf(OpenedFile,"%f %f %f ",&NewDisulfideBond->ha_Chi2,&NewDisulfideBond->hh_SS,&NewDisulfideBond->Energy_SS);
         if (DEBUG)
           printf("reached C\n");
         /* scan the name and then dissect it */
         fscanf(OpenedFile,"%s ", NewDisulfideBond->Name_Pdb_File);
         /* get the part without the extension .pdb*/
         for (i = 0; NewDisulfideBond->Name_Pdb_File[i] != '.'; i++)
         { 
           NewDisulfideBond->Name_Pdb_File_noextension[i] = NewDisulfideBond->Name_Pdb_File[i];
           NewDisulfideBond->Name_Pdb_File_noextension[i+1] = '\0';
         }
         if (DEBUG)
           printf("reached E\n");
         /*get going from front to end the rest of the file */
         for (i = 0; NewDisulfideBond->Name_Pdb_File[i] ==NewDisulfideBond->Name_Starting_Structure[i]; i++)
         {
           if (DEBUG)
             printf("%c",NewDisulfideBond->Name_Pdb_File[i]);
         }
         /* name donor strand */
         NewDisulfideBond->NameDonorStrand[0] = NewDisulfideBond->Name_Pdb_File[++i];
         NewDisulfideBond->NameDonorStrand[1] = '\0';
         if (DEBUG)
           printf("Name Donor Strand %s\n", NewDisulfideBond->NameDonorStrand); 
         
         /* name donor original amino acid */
         i++;
         NewDisulfideBond-> CurrentAminoAcidDonor[0]= NewDisulfideBond->Name_Pdb_File[++i];
         NewDisulfideBond-> CurrentAminoAcidDonor[1]= '\0';
         if (DEBUG)
           printf("Name Original amino acid donor %s\n",NewDisulfideBond-> CurrentAminoAcidDonor); 
         
         /* number donor residue */
         for (i++, j = 0; NewDisulfideBond->Name_Pdb_File[i] != 'C'; i++, j++)
         {
           CharNumbers[j] = NewDisulfideBond->Name_Pdb_File[i];  
           CharNumbers[j+1] = '\0';
         }
         NewDisulfideBond->ResidueNumberDonorStrand = atoi(CharNumbers);
         if (DEBUG)
           printf("Just identified %d as the number of the first mutation. \n", NewDisulfideBond->ResidueNumberDonorStrand );
         
         /* name strand acceptor */
         i++;
         NewDisulfideBond->NameAcceptorStrand[0] = NewDisulfideBond->Name_Pdb_File[++i];
         NewDisulfideBond->NameAcceptorStrand[1] = '\0';
         if (DEBUG)
           printf("Name acceptor Strand %s\n", NewDisulfideBond->NameAcceptorStrand); 
         
         /* name original residue acceptor */
         i++;
         NewDisulfideBond-> CurrentAminoAcidAcceptor[0]= NewDisulfideBond->Name_Pdb_File[++i];
         NewDisulfideBond-> CurrentAminoAcidAcceptor[1]= '\0';
         if (DEBUG)
           printf("Name Original amino acid acceptor residue %s\n",NewDisulfideBond-> CurrentAminoAcidAcceptor); 
         
         /* number original residue acceptor*/
         for (i++, j = 0; NewDisulfideBond->Name_Pdb_File[i] != 'C'; i++, j++)
         {
           CharNumbers[j] = NewDisulfideBond->Name_Pdb_File[i];  
           CharNumbers[j+1] = '\0';
         }
         NewDisulfideBond->ResidueNumberAcceptorStrand = atoi(CharNumbers);
         if (DEBUG)
           printf("Just identified %d as the number of the second mutation. \n", NewDisulfideBond->ResidueNumberAcceptorStrand );
         
         /*geometry description */
         i++;
         for (i++, j = 0; NewDisulfideBond->Name_Pdb_File[i] != '.'; i++, j++)
         {
           NewDisulfideBond-> AbbreviatedGeometry[j] = NewDisulfideBond->Name_Pdb_File[i];  
           NewDisulfideBond-> AbbreviatedGeometry[j+1] = '\0';
         }
         if (DEBUG)
           printf("Just identified %s as the geometry description.\n",NewDisulfideBond->AbbreviatedGeometry);
         
         /* interchain disulfide */
         if (NewDisulfideBond->NameDonorStrand[0] == NewDisulfideBond->NameAcceptorStrand[0])
           NewDisulfideBond->InterChainDisulfide = FALSE;
         else
           NewDisulfideBond->InterChainDisulfide = TRUE;
         
         if (DEBUG)
         {
           printf("Just identified %s as the pdb file without extension.\n",NewDisulfideBond->Name_Pdb_File_noextension);
         }
         /* if debug, report small part of what is scanned */
         if (DEBUG)
           printf("JustScanned %.3f as dCA_dCA distance and %.3f as energy \n\n\n", NewDisulfideBond->d_CA_CA, NewDisulfideBond->Energy_SS );

       
         /* now connect it to the previous disulfide bond in the list */
         NewDisulfideBond->next = Before;
         Before        = NewDisulfideBond;
      }
      else 
      {  
         fgets (CharforInfolineBelow, 1000, OpenedFile);
         if (DEBUG)
           printf("Assumed this was unimportant: %s %s",CharTemporyForDeciding, CharforInfolineBelow);
      }
    }
    while (strcmp(CharTemporyForDeciding, "d") != 0);
  } while((c = fgetc(OpenedFile)) != EOF);

  
  fclose(OpenedFile);

     
  return Before;  
  
}


/*the sort algorithm is slow for large numbers (speed goes as 1/2N^2) but large number of disulfide bonds are not expected and this algorithm was easiest to program. */
struct DisulfideBond * SortDisulfideList(struct DisulfideBond *NameDisulfideList)
{
  /* first count them */
  int NumberOfBridgesToSort = 0;
  int DecideToSwitch, Intradimer_p, Intradimer_q, same_disulfide, same_conformation; 
  
  struct DisulfideBond *p;
  for (p = NameDisulfideList; p!=NULL; p = p->next)
  {
    NumberOfBridgesToSort++;  
  }
  /*The algorithm goes throught the entire list and every time the last element sorted is surely at its correct position and does not need to be sorted again */
  struct DisulfideBond *q, *r,*ThePrevious;
  
  int i;
  for (; NumberOfBridgesToSort > 0; NumberOfBridgesToSort--)
  {
    for(i = 0,p = NameDisulfideList; i <( NumberOfBridgesToSort - 1); p = p->next, i++) 
    {
      q    = p->next;
      /*reset decisions and characterizations to switch to zero */
      DecideToSwitch    = FALSE;
      Intradimer_p      = FALSE;
      Intradimer_q      = FALSE;
      same_disulfide    = FALSE;
      same_conformation = FALSE; 
      
      /* we want to have inter dimers at the end */
      if (strcmp(p->NameDonorStrand, p->NameAcceptorStrand))
        Intradimer_p      = TRUE;
      if (strcmp(q->NameDonorStrand, q->NameAcceptorStrand))
        Intradimer_q      = TRUE;
      if (Intradimer_p > Intradimer_q)
      {  
        DecideToSwitch = TRUE;
	if (DEBUG)
	  printf("Decided to switch %d since %s %s and %s %s\n", i, p->NameDonorStrand, p->NameAcceptorStrand, q->NameDonorStrand, q->NameAcceptorStrand);
      }
      
      /* next sorting steps only if no difference in intra dimer status*/
      if (Intradimer_p == Intradimer_q)
      {
        /* decide if it is the same disulfide */
        if (p->ResidueNumberDonorStrand == q->ResidueNumberDonorStrand)
	  if (p->ResidueNumberAcceptorStrand == q->ResidueNumberAcceptorStrand)
	    same_disulfide    = TRUE;  
        if (same_disulfide == FALSE)
        {
          if (p->ResidueNumberDonorStrand >  q->ResidueNumberDonorStrand)
	    DecideToSwitch    = TRUE;
          if (p->ResidueNumberDonorStrand ==  q->ResidueNumberDonorStrand)
	    if (p->ResidueNumberAcceptorStrand > q->ResidueNumberAcceptorStrand)
	      DecideToSwitch    = TRUE;
        }
        
	/* next sorting steps only if it is the same disulfide */
	if (same_disulfide == TRUE)
	{
	  /* switch based on strcmp difference*/ 
	  if (strcmp(p->AbbreviatedGeometry, q->AbbreviatedGeometry) == 0)
           same_conformation = TRUE;
          else 
	  {
	    if (strcmp(p->AbbreviatedGeometry, q->AbbreviatedGeometry) < 0)
	      DecideToSwitch = TRUE;
	  }
	    
	  /* next sorting steps only if same conformation*/
	  if (same_conformation == TRUE)
	  {
	    if (DEBUG)
              printf("%d %.3f compared to the following %.3f\n", i, p-> Energy_SS , q->Energy_SS);
            if (p->Energy_SS > q->Energy_SS) /*i.e. if the energy of the first is higher than that of the latter one */
              DecideToSwitch = TRUE;
	      
          }
        }
      }
      if (DecideToSwitch == TRUE)
      { 
  	r    = q->next;
        if (i == 0) /*special case if the first point, then need to modify the starting address of the list as well */
        {
	  NameDisulfideList        = q;
	  p->next                  = r;
	  q->next                  = p;
	  p                        = q;
        }
        else
        {  
          ThePrevious->next        = q;
	  p->next                  = r;
	  q->next                  = p;
	  p                        = q;
	}
      }
      ThePrevious = p;
     }
  }  
  
  return NameDisulfideList;
  
}


/* This removes the stored list of disulfide bonds out of memory */
struct DisulfideBond * RemoveDisulfidesFromMemory(struct DisulfideBond * NameDisulfideList)
{
  struct DisulfideBond * p;
  struct DisulfideBond * Previous;
  
  p = NameDisulfideList;
  do 
  {
    Previous = p;
    p = p->next;
    free(Previous);  
  } while (p != NULL);
  
  return p;
}


/* count number of disulfide bridges */
int CountNumberDisulfideBridges(struct DisulfideBond *NameListDisulfides)
{
  /* check if input is OK */
  if (NameListDisulfides == NULL)
    error_message("Called disulfide list which does not (longer) exist ");
  
  int NumberOfBridges = 0;
  
  
  struct DisulfideBond *p;
  for (p = NameListDisulfides; p!=NULL; p = p->next)
  {
    NumberOfBridges++;  
  }
  
  return(NumberOfBridges);
}


/* this prints the list of disulfide bridges */
void PrintDisulfideBridges(struct DisulfideBond *NameListDisulfides)
{
  if (NameListDisulfides == NULL)
    printf("Called mcystein list which is empty or does not longer exist\n");
  
  struct DisulfideBond *p, *Previous;
  int PrintedBefore = FALSE;
  for (p = NameListDisulfides; p != NULL; p = p->next)
  {
    if (PrintedBefore == TRUE)
    {  if ( (p->ResidueNumberDonorStrand == Previous->ResidueNumberDonorStrand) && (p->ResidueNumberAcceptorStrand == Previous->ResidueNumberAcceptorStrand))
         if (strcmp(p->AbbreviatedGeometry, Previous->AbbreviatedGeometry) != 0)
           printf("................................................................\n");
    } 
    if (PrintedBefore == TRUE)
    {  if (!( (p->ResidueNumberDonorStrand == Previous->ResidueNumberDonorStrand) && (p->ResidueNumberAcceptorStrand == Previous->ResidueNumberAcceptorStrand)))
         printf("----------------------------------------------------------------\n");
    } 
    if (PrintedBefore == TRUE)
    {  if ( p->InterChainDisulfide == TRUE && Previous->InterChainDisulfide != TRUE)
         printf("================================================================\n");
    } 
    printf("%s %s %d C ",p->NameDonorStrand,p->CurrentAminoAcidDonor ,p->ResidueNumberDonorStrand);
    printf("%s %s %d C ",p->NameAcceptorStrand,p->CurrentAminoAcidAcceptor ,p->ResidueNumberAcceptorStrand);
    printf("%s %.3f ",p->AbbreviatedGeometry, p->Energy_SS);  
    printf("%s/%s/%s ", p->OriginalDirectory, p->Name_Pdb_File_noextension, p->Name_Pdb_File);
    if (p->InterChainDisulfide == TRUE)
      printf ("INTER");
    printf("\n");  
    Previous = p;
    PrintedBefore = TRUE;
  }
}


/* this prints the list of unique disulfide bridges */
void PrintDisulfideBridgesUnique(struct DisulfideBond *NameListDisulfides)
{
  if (NameListDisulfides == NULL)
    printf("Called cystein list which is empty or does not longer exist\n");
  int PrintCurrent, PrintedBefore = FALSE;
  struct DisulfideBond *p, *Previous;
  for (p = NameListDisulfides; p != NULL; p = p->next)
  {
    if (PrintedBefore == FALSE)
      PrintCurrent = TRUE;
    else
    {
      PrintCurrent = TRUE;
      if ((p->ResidueNumberDonorStrand) == (Previous->ResidueNumberDonorStrand))
        if ((p->ResidueNumberAcceptorStrand) == (Previous->ResidueNumberAcceptorStrand))
	  if (strcmp(p->AbbreviatedGeometry,Previous->AbbreviatedGeometry) == 0)
            PrintCurrent = FALSE;
    }    
    if (PrintCurrent == TRUE)
    {
      /* put handy separators for if printed before */
      if (PrintedBefore == TRUE)
      {  if (!( (p->ResidueNumberDonorStrand == Previous->ResidueNumberDonorStrand) && (p->ResidueNumberAcceptorStrand == Previous->ResidueNumberAcceptorStrand)))
           printf("----------------------------------------------------------------\n");
      } 
      if (PrintedBefore == TRUE)
      {  if ( p->InterChainDisulfide == TRUE && Previous->InterChainDisulfide != TRUE)
           printf("================================================================\n");
      } 
      printf("%s %s %d C ",p->NameDonorStrand,p->CurrentAminoAcidDonor ,p->ResidueNumberDonorStrand);
      printf("%s %s %d C ",p->NameAcceptorStrand,p->CurrentAminoAcidAcceptor ,p->ResidueNumberAcceptorStrand);
      printf("%s %.3f ",p->AbbreviatedGeometry, p->Energy_SS);  
      printf("%s/%s/%s ", p->OriginalDirectory, p->Name_Pdb_File_noextension, p->Name_Pdb_File);
      if (p->InterChainDisulfide == TRUE)
        printf ("INTER");
      printf("\n");  
      Previous = p;
      PrintedBefore = TRUE;
    }
  }
}


void WriteDisulfideBridgesUnique(struct DisulfideBond *NameListDisulfides)
{
  if (NameListDisulfides == NULL)
    printf("Called cystein list which is empty or does not longer exist\n");
  
  int NumberOfConformations = 0; /*different conformations of all disulfide bonds */
  int NumberOfUniqueBonds   = 1; /*differently connected disulfide bonds */
  
  FILE *BestEnergyUniqueDisulfideBonds;
  BestEnergyUniqueDisulfideBonds = fopen("BestEnergyUniqueDisulfideBonds.tab", "w");
  
  int PrintCurrent, PrintedBefore = FALSE;
  struct DisulfideBond *p, *Previous;
   
  fprintf(BestEnergyUniqueDisulfideBonds, "Residues        Geometry          INTER/Intra   Name best energy pdb file                                                               "     );   
  fprintf(BestEnergyUniqueDisulfideBonds, "d_CA_CA    d_CB_CB    d_SG_SG adCA_CB_SG aaCA_CB_SG adCB_SG_SG aaCB_SG_SG    hd_Chi1    ha_Chi1    hd_Chi2    ha_Chi2      hh_SS  Energy_SS \n");

  for (p = NameListDisulfides; p != NULL; p = p->next)
  {
    if (PrintedBefore == FALSE)
      PrintCurrent = TRUE;
    else
    {
      PrintCurrent = TRUE;
      if ((p->ResidueNumberDonorStrand) == (Previous->ResidueNumberDonorStrand))
        if ((p->ResidueNumberAcceptorStrand) == (Previous->ResidueNumberAcceptorStrand))
	  if (strcmp(p->AbbreviatedGeometry,Previous->AbbreviatedGeometry) == 0)
            PrintCurrent = FALSE;
    }    
    if (PrintCurrent == TRUE)
    {
      NumberOfConformations++;
      /* put handy separators for if printed before */
      if (PrintedBefore == TRUE)
      {  if (!( (p->ResidueNumberDonorStrand == Previous->ResidueNumberDonorStrand) && (p->ResidueNumberAcceptorStrand == Previous->ResidueNumberAcceptorStrand)))
         {
	   fprintf(BestEnergyUniqueDisulfideBonds, "----------------------------------------------------------------\n");
	   NumberOfUniqueBonds++;
         }
      } 
      if (PrintedBefore == TRUE)
      {  if ( p->InterChainDisulfide == TRUE && Previous->InterChainDisulfide != TRUE)
           fprintf(BestEnergyUniqueDisulfideBonds, "================================================================\n");
      } 
      fprintf(BestEnergyUniqueDisulfideBonds, "%s%dC/%s%dC\t",p->CurrentAminoAcidDonor ,p->ResidueNumberDonorStrand,p->CurrentAminoAcidAcceptor ,p->ResidueNumberAcceptorStrand);
      fprintf(BestEnergyUniqueDisulfideBonds, "%-18s\t",p->AbbreviatedGeometry);  
     
      if (p->InterChainDisulfide == TRUE)
        fprintf(BestEnergyUniqueDisulfideBonds,"INTER\t");
      else
        fprintf(BestEnergyUniqueDisulfideBonds,"Intra\t");
      fprintf(BestEnergyUniqueDisulfideBonds,"%-80s\t", p->Name_Pdb_File);
      
      fprintf(BestEnergyUniqueDisulfideBonds, "%8.3f  %8.3f  %8.3f   %8.3f   %8.3f", p->d_CA_CA, p->d_CB_CB,p->d_SG_SG, p->adCA_CB_SG,p->aaCA_CB_SG);  
      fprintf(BestEnergyUniqueDisulfideBonds, "   %8.3f  %8.3f  %8.3f   %8.3f   ", p->adCB_SG_SG, p->aaCB_SG_SG, p->hd_Chi1, p->ha_Chi1);  
      fprintf(BestEnergyUniqueDisulfideBonds, "   %8.3f   %8.3f   %8.3f   %8.3f   ",  p->hd_Chi2, p->ha_Chi2, p->hh_SS,  p->Energy_SS);  

      fprintf(BestEnergyUniqueDisulfideBonds,"\n");  
      Previous = p;
      PrintedBefore = TRUE;
    }
  }
  printf("The number of Disulfide bridge conformations is approximately .........%d\n",NumberOfConformations);
  printf("The number of uniquely connected disulfide bridges is approximately....%d\n",NumberOfUniqueBonds);
  fclose(BestEnergyUniqueDisulfideBonds);
}




void WriteDisulfideBridges(struct DisulfideBond *NameListDisulfides)
{
  if (NameListDisulfides == NULL)
    printf("Called cystein list which is empty or does not longer exist\n");
  
  int NumberOfConformations = 0; /*different conformations of all disulfide bonds */
  int NumberOfUniqueBonds   = 1; /*differently connected disulfide bonds */
  
  FILE *AllDisulfideBonds;
  AllDisulfideBonds = fopen("AllDisulfideBondsSorted.tab", "w");
  
  int PrintCurrent, PrintedBefore = FALSE;
  struct DisulfideBond *p, *Previous;
   
  fprintf(AllDisulfideBonds, "Residues        Geometry          INTER/Intra   Name best energy pdb file                                                               "     );   
  fprintf(AllDisulfideBonds, "d_CA_CA    d_CB_CB    d_SG_SG adCA_CB_SG aaCA_CB_SG adCB_SG_SG aaCB_SG_SG    hd_Chi1    ha_Chi1    hd_Chi2    ha_Chi2      hh_SS  Energy_SS \n");

  for (p = NameListDisulfides; p != NULL; p = p->next)
  {
    if (PrintedBefore == FALSE)
      PrintCurrent = TRUE;
    else
    {
      PrintCurrent = TRUE;
      if ((p->ResidueNumberDonorStrand) == (Previous->ResidueNumberDonorStrand))
        if ((p->ResidueNumberAcceptorStrand) == (Previous->ResidueNumberAcceptorStrand))
	  if (strcmp(p->AbbreviatedGeometry,Previous->AbbreviatedGeometry) == 0)
            PrintCurrent = FALSE;
    }    
    NumberOfConformations++;
    /* put handy separators for if printed before */
    if (PrintedBefore == TRUE)
    {  if ( (p->ResidueNumberDonorStrand == Previous->ResidueNumberDonorStrand) && (p->ResidueNumberAcceptorStrand == Previous->ResidueNumberAcceptorStrand))
         if (strcmp(p->AbbreviatedGeometry, Previous->AbbreviatedGeometry) != 0)
           fprintf(AllDisulfideBonds,"................................................................\n");
    } 

    if (PrintedBefore == TRUE)
    {  if (!( (p->ResidueNumberDonorStrand == Previous->ResidueNumberDonorStrand) && (p->ResidueNumberAcceptorStrand == Previous->ResidueNumberAcceptorStrand)))
       {
	 fprintf(AllDisulfideBonds, "----------------------------------------------------------------\n");
	 NumberOfUniqueBonds++;
       }
    } 
    if (PrintedBefore == TRUE)
    {  if ( p->InterChainDisulfide == TRUE && Previous->InterChainDisulfide != TRUE)
         fprintf(AllDisulfideBonds, "================================================================\n");
    } 
    fprintf(AllDisulfideBonds, "%s%dC/%s%dC\t",p->CurrentAminoAcidDonor ,p->ResidueNumberDonorStrand,p->CurrentAminoAcidAcceptor ,p->ResidueNumberAcceptorStrand);
    fprintf(AllDisulfideBonds, "%-18s\t",p->AbbreviatedGeometry);  
   
    if (p->InterChainDisulfide == TRUE)
      fprintf(AllDisulfideBonds,"INTER\t");
    else
      fprintf(AllDisulfideBonds,"Intra\t");
    fprintf(AllDisulfideBonds,"%-80s\t", p->Name_Pdb_File);
    
    fprintf(AllDisulfideBonds, "%8.3f  %8.3f  %8.3f   %8.3f   %8.3f", p->d_CA_CA, p->d_CB_CB,p->d_SG_SG, p->adCA_CB_SG,p->aaCA_CB_SG);  
    fprintf(AllDisulfideBonds, "   %8.3f  %8.3f  %8.3f   %8.3f   ", p->adCB_SG_SG, p->aaCB_SG_SG, p->hd_Chi1, p->ha_Chi1);  
    fprintf(AllDisulfideBonds, "   %8.3f   %8.3f   %8.3f   %8.3f   ",  p->hd_Chi2, p->ha_Chi2, p->hh_SS,  p->Energy_SS);  

    fprintf(AllDisulfideBonds,"\n");  
    Previous = p;
    PrintedBefore = TRUE;
  }
  fclose(AllDisulfideBonds);
}

 
/* this prints the cp commands */
void WritecpCommands(struct DisulfideBond *NameListDisulfides)
{
  if (NameListDisulfides == NULL)
    printf("Called cystein list which is empty or does not longer exist\n");
  int PrintCurrent, PrintedBefore = FALSE;
  struct DisulfideBond *p, *Previous;
  
  system("mkdir UniqueDisulfides");
  
  FILE *tocplist;
  tocplist = fopen("tocplist", "w");
  
  for (p = NameListDisulfides; p != NULL; p = p->next)
  {
    if (PrintedBefore == FALSE)
      PrintCurrent = TRUE;
    else
    {
      PrintCurrent = TRUE;
      if ((p->ResidueNumberDonorStrand) == (Previous->ResidueNumberDonorStrand))
        if ((p->ResidueNumberAcceptorStrand) == (Previous->ResidueNumberAcceptorStrand))
	  if (strcmp(p->AbbreviatedGeometry,Previous->AbbreviatedGeometry) == 0)
            PrintCurrent = FALSE;
    }    
    if (PrintCurrent == TRUE)
    {
      fprintf(tocplist,"mkdir UniqueDisulfides/Subdir_templates\n");
      if ((p->ResidueNumberDonorStrand) <= (p->ResidueNumberAcceptorStrand))
      {  
        fprintf(tocplist,"mkdir UniqueDisulfides/Subdir_%s%dC_%s%dC\n",p->CurrentAminoAcidDonor,p->ResidueNumberDonorStrand,p->CurrentAminoAcidAcceptor,p->ResidueNumberAcceptorStrand);
        fprintf(tocplist,"cp %s/%s UniqueDisulfides/Subdir_%s%dC_%s%dC\n", p->OriginalDirectory, p->Name_Pdb_File,p->CurrentAminoAcidDonor,p->ResidueNumberDonorStrand,p->CurrentAminoAcidAcceptor,p->ResidueNumberAcceptorStrand);
      }
      else
      {  
        fprintf(tocplist,"mkdir UniqueDisulfides/Subdir_%s%dC_%s%dC\n",p->CurrentAminoAcidAcceptor,p->ResidueNumberAcceptorStrand,p->CurrentAminoAcidDonor,p->ResidueNumberDonorStrand);
        fprintf(tocplist,"cp %s/%s UniqueDisulfides/Subdir_%s%dC_%s%dC\n", p->OriginalDirectory, p->Name_Pdb_File,p->CurrentAminoAcidAcceptor,p->ResidueNumberAcceptorStrand,p->CurrentAminoAcidDonor,p->ResidueNumberDonorStrand);
      }
        
      fprintf(tocplist,"cp %s.pdb UniqueDisulfides/Subdir_templates/\n", p->Name_Starting_Structure);
      Previous = p;
      PrintedBefore = TRUE;
    }
    
    
    
  }
  fclose(tocplist);
  
  system("sort -r tocplist| uniq >tmp");
  system("mv tmp tocplist");
  system("chmod +x tocplist");
}


