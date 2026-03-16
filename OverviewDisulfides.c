/**************************************************************************************************************************************************
*                                                                                                                                                 *
*                        OverviewDisulfides.c                                                                                                     *
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "OverviewDisulfides.h"

int main(int argc, char *argv[])
{
  /*first important declarations*/
  int NumberOfDisulfideBridges;                           /* The number of disulfide bridges */
  
  /* provide the user with info about what the program does and examples of its use*/
  short_info();
  example_use();
 
  
  /*Analyse the existing disulfide bridges*/
  /* ============================================================ */
  if (DEBUG)
    printf("Will now start making the raw list\n"); 
  /* first create a raw list from the subdirectories */
  CreateRawDisulfideList();
  if (DEBUG)
    printf("Finished creating raw list\n");
 
  /* read in this list*/ 
  if (DEBUG)
    printf("Start reading raw list\n");
  struct DisulfideBond * InitialList = LoadListFromFile("RawDisulfideListtmp.tab");
  if (DEBUG)
    printf("Still need code here\n"); 
  NumberOfDisulfideBridges = CountNumberDisulfideBridges(InitialList);
  printf("The number of read in Disulfide bridges is ............................%d\n", NumberOfDisulfideBridges);
  
  
  /* and order it based on 
  A) is it an interchain or intrachain disulfide bond 
  B) what is the number of its first amino acid
  C) what is its basic geometry
  D) what is its energy (best is lowest is first)
  */
  
  struct DisulfideBond * SortedList = SortDisulfideList(InitialList);
  
  if (DEBUG)
  {  /* this prints the list of disulfide bridges */
    PrintDisulfideBridges(SortedList);
  
    printf("And now the unique ones\n");
    PrintDisulfideBridgesUnique(SortedList);
  }
  
  /*write a file with all the unique disulfide bonds tabulated */
  WriteDisulfideBridges(SortedList);
  WriteDisulfideBridgesUnique(SortedList);
  /*make a file to copy all disulfide bond files to a separate directory and execute this file*/
  WritecpCommands(SortedList);
  system("./tocplist");  
  /* free up memory before ending program */
  InitialList =  RemoveDisulfidesFromMemory(InitialList);
  
  
  /*finish program*/
  printf("\nFinished program succesfully\n\n");
  exit(EXIT_SUCCESS);
}
