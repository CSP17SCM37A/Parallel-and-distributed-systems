/* Compile the below code with "mpicc hw6.c -o hw6 -lm" 
 and run "mpirun -np 4 ./hw6  output1.txt 2 output.txt" (Here we need to give three command line arguments i.e outputfile1,MB for Rank, and output file name for timings).
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


/* file for timings */
char filenameR[16]; 
/* file for rank */
char filenameT[16];
/* mb for rank variable declaration */ 
int p;

/*Reading ARGUMENTS */
void parameters(int argc, char **argv) {

if (argc == 4)
{
    
    strcpy(filenameR, argv[3]);
    printf("filename for timings  = %s \n",filenameR);
 } 

if (argc >= 3) 
{
    p = atoi(argv[2]);
   
    printf("MB PER RANK          = %i \n", p);
} 

if (argc >= 2)
{
    strcpy(filenameT, argv[1]);
    printf("filename for rank     = %s \n",filenameT);
    
}

            
     
}



int main(int argc, char **argv)
{

int i;
int rank, size, nints;
long  offset;

/* 1 mb=n */ 
long N=512*512/2;


/*writing timings variable */
double t1, t2; 

/*overall timings variable */
double initial,final;

//Error handling variable */
int  error, eclass,len;
char estring[MPI_MAX_ERROR_STRING];


/*read the arguments with this function */
parameters(argc, argv);/* Process program parameters */

FILE * f; // to write timings in output file
f = fopen(filenameR,"w");  //reading a file

if(f == NULL){
printf("Error : File not found\n");
return 0;
}
     
/**Initializing MOI_TIME() */
initial= MPI_Wtime();
/*Declaring file name */
MPI_File fhw;

MPI_Status status;

/*Declaring MPI Intilization */
MPI_Init(&argc, &argv);

MPI_Comm_rank(MPI_COMM_WORLD, &rank);

MPI_Comm_size(MPI_COMM_WORLD, &size);

/*CALCULATING TOTAL SIZE FOR EACH RANK */
N=N*p;

printf("Total Size is %ld \n",N);

/****************************************Writing Started ******************************************************/
long buf[N];


/**Initializing MPI_TIME() FOR WRITING*/
t1 = MPI_Wtime();

printf("Writing Started\n");

/*STORING DATA INTO BUFFER*/
for ( i=0;i<N;i++)
{

  buf[i] = 1;

}

/*CALCULATING OFFSET*/
offset = rank*(N)*sizeof(long);

/*opening a file */
error=MPI_File_open(MPI_COMM_WORLD, "123.txt",MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fhw);   /*Error handling in opening a file */


if(error != 0){
    MPI_Error_class(error, &eclass);
    MPI_Error_string(error, estring, &len);
    printf("Error %d: %s\n", eclass, estring);
    fflush(stdout); }
/*printing rank and offset */
printf("\nRank: %d, Offset: %ld\n", rank, offset);
/*printing rank and offset */
error=MPI_File_set_view(fhw, offset, MPI_LONG, MPI_LONG, "native", MPI_INFO_NULL); /*Error handling in opening a file */
if(error != 0){
    MPI_Error_class(error, &eclass);
    MPI_Error_string(error, estring, &len);
    printf("Error %d: %s\n", eclass, estring);
    fflush(stdout); }
/*Writing a data in to a file */
error=MPI_File_write(fhw, buf, N, MPI_LONG, &status); /*Error handling in  writing data into a file */
if(error != 0){
    MPI_Error_class(error, &eclass);
    MPI_Error_string(error, estring, &len);
    printf("Error %d: %s\n", eclass, estring);
    fflush(stdout); }
/*closing a file*/
error=MPI_File_close(&fhw);  /*Error handling in closing a file */
if(error != 0){
    MPI_Error_class(error, &eclass);
    MPI_Error_string(error, estring, &len);
    printf("Error %d: %s\n", eclass, estring);
    fflush(stdout); }


/**ENDING  MTI_TIME() FOR WRITING*/
t2 = MPI_Wtime();

printf("Writing Ends\n");

//fprintf(f,"\nElapsed time = %g ms.\n",(float)(usecstop - usecstart)/(float)1000);
printf( "Elapsed time for rank %d writing is(Seconds) %f\n",rank, t2 - t1 ); 
fprintf(f, "Elapsed time for writing is(Seconds) %f\n", t2 - t1 ); 
/****************************************Writing Ends ******************************************************/

/****************************************Reading Started ******************************************************/
double t3, t4;
long  bufsize;
MPI_File fh;

MPI_Status status1;

printf("Reading Started\n");

/**starting  MTI_TIME() FOR reading*/
t3 = MPI_Wtime();
/*calculating size to read data from a file*/

bufsize= N/size;
nints= bufsize/sizeof(long);
long buf1[nints];
/*opening a file */
error=MPI_File_open(MPI_COMM_WORLD,"123.txt",MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);   /*Error handling in opening a file */

if(error != 0){
    MPI_Error_class(error, &eclass);
    MPI_Error_string(error, estring, &len);
    printf("Error %d: %s\n", eclass, estring);
    fflush(stdout); }
/*Reading a file*/
error=MPI_File_read_at(fh, rank*bufsize, buf, nints, MPI_LONG, &status1);  /*Error handling in reading a file */
if(error != 0){
    MPI_Error_class(error, &eclass);
    MPI_Error_string(error, estring, &len);
    printf("Error %d: %s\n", eclass, estring);
    fflush(stdout); }
printf("\nrank: %d, buf[%ld]: %ld \n", rank, rank*bufsize, buf[rank*bufsize]);  

/*closing a file*/
error=MPI_File_close(&fh);   /*Error handling in closing a file */

if(error != 0){
    MPI_Error_class(error, &eclass);
    MPI_Error_string(error, estring, &len);
    printf("Error %d: %s\n", eclass, estring);
    fflush(stdout); }
/**ENDING  MTI_TIME() FOR reading*/
t4 = MPI_Wtime();
printf("Reading ends\n");
//fprintf(f,"\nElapsed time = %g ms.\n",(float)(usecstop - usecstart)/(float)1000);
fprintf(f, "Elapsed time for rank %d Reading  is(Seconds) %f\n",rank, t4 - t3 ); 
printf( "Elapsed time for rank %d Reading  is(Seconds) %f\n",rank, t4 - t3 ); 
MPI_Finalize();

/**ENDING  MTI_TIME() FOR overall time*/
final= MPI_Wtime();
/*************************************Reading Ends***********************************************************/
/*Printing timings in a file */
fprintf(f, "overall time (Seconds) %f\n", final - initial );
printf( "overall time (Seconds) %f\n", final - initial );
int bandwidth;

bandwidth= p/(final - initial);
printf( "Bandwidth id (MB/sec) %d\n", bandwidth );
fprintf(f,"Bandwidth id (MB/sec) %d\n", bandwidth );
/*************************************Deleting FIle***********************************************************/
 /*Deleting a file */ 
int ret = remove("123.txt");

   if(ret == 0) 
   {
      printf("*******File deleted successfully *******\n");
   }
   
/*************************************Deleting FIle***********************************************************/


return 0;

}
