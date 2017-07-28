/*Compile the below code with "mpicc first.c -o first -lm

    and run "mpiexec -np 16 ./first abc.txt " (Here we need to give one command line argument i.e outputfile for rank).
 
    To check output of a file " od -td -v abc.txt "
*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <mpi.h>

/*FIle name to print Ranks */
char filenameR[16];  

/*Reading ARGUMENTS */
void parameters(int argc, char **argv) {

if (argc >= 2)
{
    strcpy(filenameR,argv[1]);
   // printf("filename for rank     = %s \n",filenameR);
    
}
        
     
}


int main(int argc, char **argv)
{


int rank, size, offset,N=16;

int  error, eclass,len;
char estring[MPI_MAX_ERROR_STRING];

  //  MPI_Init(&argc,&argv);
  
/*read the arguments with this function */
parameters(argc, argv);
    
/*Declaring file name */
MPI_File fhw;

MPI_Status status;
/*Declaring MPI Intilization */

MPI_Init(&argc, &argv);

MPI_Comm_rank(MPI_COMM_WORLD, &rank);

MPI_Comm_size(MPI_COMM_WORLD, &size);
/*Error Handler */
MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

/****************************************Writing Started ******************************************************/
printf("Writing Started\n");

int buf[4]; // Size of the buffer is 4

/*Offset calculation */
offset = (4*rank)*sizeof(int);

/*      buf[0]     buf[1](rank) buf[2]       buf[3]
          0           0           0           0
         16           1           0           0
         32           2           0           0
         48           3           0           0 */


buf[0] =offset; 
buf[1] =rank;
buf[2] =0;
buf[3] =0;

/*opening a file */
error =MPI_File_open(MPI_COMM_WORLD, filenameR, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fhw); //checking error in opening  a file
  if(error != 0){
    MPI_Error_class(error, &eclass);
    MPI_Error_string(error, estring, &len);
    printf("Error %d: %s\n", eclass, estring);
    fflush(stdout); }

/*printing rank and offset */
printf("\nRank: %d, Offset: %d\n", rank, offset);

/*writing offset and each rank in a file */
error =MPI_File_write_at(fhw, offset, buf,4, MPI_INT, &status); //checking error in writing  a file
if(error != 0){
    MPI_Error_class(error, &eclass);
    MPI_Error_string(error, estring, &len);
    printf("Error %d: %s\n", eclass, estring);
    fflush(stdout); }
/*closing a file*/
error =MPI_File_close(&fhw);  //checking error for closing a file
if(error != 0){
    MPI_Error_class(error, &eclass);
    MPI_Error_string(error, estring, &len);
    printf("Error %d: %s\n", eclass, estring);
    fflush(stdout); }

printf("Writing Ends\n");


/****************************************Writing Ends ******************************************************/

MPI_Finalize();



return 0;

}


/* outPut:
lab-4@lab4-MS-7693 ~/Desktop $ mpiexec -np 16 ./first abc.txt 
lab-4@lab4-MS-7693 ~/Desktop $ od -td -v abc.txt
0000000           0           0           0           0
0000020          16           1           0           0
0000040          32           2           0           0
0000060          48           3           0           0
0000100          64           4           0           0
0000120          80           5           0           0
0000140          96           6           0           0
0000160         112           7           0           0
0000200         128           8           0           0
0000220         144           9           0           0
0000240         160          10           0           0
0000260         176          11           0           0
0000300         192          12           0           0
0000320         208          13           0           0
0000340         224          14           0           0
0000360         240          15           0           0
0000400
*/


