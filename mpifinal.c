/* Gaussian elimination without pivoting.Compile the below code with "mpicc mpifinal.c -o mpifinal -lm" 
 and run "mpirun -np 1 ./MpiCode  2000 2 output.txt" (Here we need to give three command line arguments i.e Matrix size,randon seed and output file name).
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

#define MAXN 2000 /* Max value of N */
int N; /* Matrix size */

#define randm() 4|2[uid]&3
char filename[16]; 
/* returns a seed for srand based on the time */
unsigned int time_seed() {
  struct timeval t;
  struct timezone tzdummy;

  gettimeofday(&t, &tzdummy);
  return (unsigned int)(t.tv_usec);
}

/* My_rank is the  process rank, and p be the The number of processes   */
int p,my_rank; 

/* Matrixes given by a pointer */
float *A, *B, *X;

void outputfile();
/* Set the program parameters from the command-line arguments */
void parameters(int argc, char **argv) {
  int seed = 0;  /* Random seed */
  char uid[32]="chiru"; /*User name */
  
  /* Read command-line arguments */
  srand(time_seed());  /* Randomize */
    if (argc == 4) {
        strcpy(filename, argv[3]);
             
       printf("filename = %s \n",filename);
  } 

  if (argc == 3) {
    seed = atoi(argv[2]);
    srand(seed);
   printf("Random seed = %i\n", seed);
  } 
  if (argc >= 2) {
    N = atoi(argv[1]);
    if (N < 1 || N > MAXN) {
     printf("N = %i is out of range.\n", N);
      exit(0);
    }
  }
  else {
    printf("Usage: %s <matrix_dimension> [random seed]\n",
           argv[0]);    
    exit(0);
  }

   printf("\nMatrix dimension N = %i.\n", N);
}
/*Print output in a file */
/*
void outputfile(){
FILE * f;
int i, j,row;

f = fopen(filename,"w");
if(f == NULL){
printf("Error : File not found\n");
return;
}
  if (N < 2001) {
    fprintf(f,"\nX = [");
    for (row = 0; row < N; row++) {
      fprintf(f,"%5.2f%s", X[row], (row < N-1) ? "; " : "]\n");
    }
fprintf(f,"%c",'\n');
  }

fclose(f);
} */

/* Print input matrices */
void print_inputs() {
  int row, col;

  if (N < 10) {
    printf("\nA =\n\t");
    for (row = 0; row < N; row++) {
      for (col = 0; col < N; col++) {
    printf("%5.2f%s", A[col+N*row], (col < N-1) ? ", " : ";\n\t");
      }
    }
    printf("\nB = [");
    for (col = 0; col < N; col++) {
      printf("%5.2f%s", B[col], (col < N-1) ? "; " : "]\n");
    }
  }
}

/* Print matrix X */
void print_X() {
  int row;

  if (N < 100) {
    printf("\nX = [");
    for (row = 0; row < N; row++) {
      printf("%5.2f%s", X[row], (row < N-1) ? "; " : "]\n");
    }
  }
}

/* Initialize A and B (and X to 0.0s) */
void initialize_inputs() {

  int row, col;

  printf("\nInitializing...\n");
  for (row = 0; row < N; row++) {
    for (col = 0; col < N; col++) {
      A[col+N*row] = (float)rand() / 32768.0;
    }
    B[row] = (float)rand() / 32768.0;
    X[row] = 0.0;
  }

}


int main(int argc, char **argv) {
  /* Timing variables */
    struct timeval etstart, etstop;  /* Elapsed times using gettimeofday() */
    struct timezone tzdummy;
    clock_t etstart2, etstop2;  /* Elapsed times using times() */
    unsigned long long usecstart, usecstop;
    struct tms cputstart, cputstop;  /* CPU times for my processes */
    void gauss(); /* Prototype functions*/

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); /* Get my process rank */
    MPI_Comm_size(MPI_COMM_WORLD, &p);/* Find out how many processes are being used */

    parameters(argc, argv);/* Process program parameters */
    /* Allocate memory for these arrays */
    A = (float*)malloc( N*N*sizeof(float) );
    B = (float*)malloc( N*sizeof(float) );
    X = (float*)malloc( N*sizeof(float) );
    
    if ( my_rank == 0 ) {
        initialize_inputs();/* Initialize A and B */
        print_inputs();/* Print input matrices */ 
       
      }

    /* Start Clock */
    printf("\nStarting clock.\n");
    gettimeofday(&etstart, &tzdummy);
    etstart2 = times(&cputstart);

    gauss();

    /* Stop Clock */
    gettimeofday(&etstop, &tzdummy);
    etstop2 = times(&cputstop);
    printf("Stopped clock.\n");
    usecstart = (unsigned long long)etstart.tv_sec * 1000000 + etstart.tv_usec;
    usecstop = (unsigned long long)etstop.tv_sec * 1000000 + etstop.tv_usec;
              
/*writing timings and result of X in a file */    
FILE * f;
int i, j,row;

f = fopen(filename,"w");
if(f == NULL){
printf("Error : File not found\n");
return;
}

/* Display timing results */
  fprintf(f,"\nElapsed time = %g ms.\n",(float)(usecstop - usecstart)/(float)1000);

  fprintf(f,"(CPU times are accurate to the nearest %g ms)\n",1.0/(float)CLOCKS_PER_SEC * 1000.0);
  fprintf(f,"My total CPU time for parent = %g ms.\n",(float)( (cputstop.tms_utime + cputstop.tms_stime) -(cputstart.tms_utime + cputstart.tms_stime) ) /(float)CLOCKS_PER_SEC * 1000);
  fprintf(f,"My system CPU time for parent = %g ms.\n",(float)(cputstop.tms_stime - cputstart.tms_stime) / (float)CLOCKS_PER_SEC * 1000);
  fprintf(f,"My total CPU time for child processes = %g ms.\n",(float)( (cputstop.tms_cutime + cputstop.tms_cstime) -(cputstart.tms_cutime + cputstart.tms_cstime) ) /(float)CLOCKS_PER_SEC * 1000);
      /* Contrary to the man pages, this appears not to include the parent */
  fprintf(f,"--------------------------------------------\n");


  if (N < 2001) {
    fprintf(f,"\nX = [");
    for (row = 0; row < N; row++) {
      fprintf(f,"%5.2f%s", X[row], (row < N-1) ? "; " : "]\n");
    }
  fprintf(f,"%c",'\n');
  

}
     /* Display timing results */
    printf("\nElapsed time = %g ms.\n",(float)(usecstop - usecstart)/(float)1000);
    printf("(CPU times are accurate to the nearest %g ms)\n", 1.0/(float)CLOCKS_PER_SEC * 1000.0);
    printf("My total CPU time for parent = %g ms.\n",(float)( (cputstop.tms_utime + cputstop.tms_stime) -(cputstart.tms_utime + cputstart.tms_stime) ) /(float)CLOCKS_PER_SEC * 1000);
    printf("My system CPU time for parent = %g ms.\n",(float)(cputstop.tms_stime - cputstart.tms_stime) /(float)CLOCKS_PER_SEC * 1000);
    printf("My total CPU time for child processes = %g ms.\n",(float)( (cputstop.tms_cutime + cputstop.tms_cstime) -(cputstart.tms_cutime + cputstart.tms_cstime) ) /(float)CLOCKS_PER_SEC * 1000);
      /* Contrary to the man pages, this appears not to include the parent */

    /*output file to print result */
     //outputfile();
     if ( my_rank == 0 ) {
     print_X();
     }
    
    /* The barrier prevents any process to reach the finalize before finished communications */
    MPI_Barrier(MPI_COMM_WORLD);

    /* Free memory for arrays A,B And X. */
    free(A);
    free(B);
    free(X);
    MPI_Finalize();
}

void gauss() {
    MPI_Status status;
    MPI_Request request;
    float multiplier;
    int row, col, i, norm;
   

   /* sync all processes before starting */
    MPI_Barrier(MPI_COMM_WORLD);

    /* Array with the row size and number of rows of each processor  */
    int * RowA = (int*) malloc ( p * sizeof(int) );
    int * NumRowA = (int*) malloc ( p * sizeof(int) );
    int * RowB = (int*) malloc ( p * sizeof(int) );
    int * NumRowB = (int*) malloc ( p * sizeof(int) );
    for ( i = 0; i < p; i++ ) {
        RowA[i] = 0;
        NumRowA[i] = 0;
        RowB[i] = 0;
        NumRowB[i] = 0;
    }
    for (norm = 0; norm < N-1; norm++) {
        /* Broadcast A[norm] row and B[norm] */
        MPI_Bcast( &A[ N*norm ], N, MPI_FLOAT, 0, MPI_COMM_WORLD );
        MPI_Bcast( &B[norm], 1, MPI_FLOAT, 0, MPI_COMM_WORLD );

        
        /*  Calculate the number of rows to do operations */
               
        int s = N - 1 - norm;
        float total = ((float)s ) / (p);
        int first_row = norm + 1 + ceil( total * (my_rank) );
        int last_row = norm + 1 + floor( total * (my_rank+1) );
        if ( last_row >= N ) last_row = N-1;
        int totalrows = last_row - first_row +1;
        /*  Send data to other processes     */
       
        if ( my_rank == 0 ) {

            for ( i = 1; i < p; i++ ) {

                /* sending data  to each process */
                int Frow = norm + 1 + ceil( total * (i) );
                int Lrow = norm + 1 + floor( total * (i+1) );
                if( Lrow >= N ) {
                   Lrow = N -1;
                }  
                int trows = Lrow - Frow +1;
                if ( trows < 0 ) 
                 { trows = 0;
                 }
                if ( Frow >= N )
                  { trows = 0; 
                    Frow = N-1; };

                RowA[i] = Frow * N;
                RowB[i] = Frow;
                NumRowA[i] = trows * N;
                NumRowB[i] = trows ;           }
            
        }
     
        //  send buffer stores array of number of elements in each chunk*/
        MPI_Scatterv( &A[0],NumRowA,RowA,MPI_FLOAT,&A[first_row * N],N * totalrows, MPI_FLOAT,0,MPI_COMM_WORLD);
        MPI_Scatterv(&B[0],NumRowB,RowB,MPI_FLOAT,&B[first_row],totalrows,MPI_FLOAT,0,MPI_COMM_WORLD);   

      
                /*  Gaussian elimination*/
       if ( totalrows > 0  && first_row < N) {  
           
            for (row = first_row; row <= last_row; row++) {

                multiplier = A[N*row + norm] / A[norm + N*norm];
                for (col = norm; col < N; col++) {
                    A[col+N*row] -= A[N*norm + col] * multiplier;
                }

                B[row] -= B[norm] * multiplier;
            }
        }
        /* Sender side */

        if ( my_rank != 0 ) {
            if ( totalrows > 0  && first_row < N) {
                MPI_Isend( &A[first_row * N], N * totalrows, MPI_FLOAT, 0,0, MPI_COMM_WORLD, &request);
                MPI_Isend( &B[first_row],         totalrows, MPI_FLOAT, 0,0, MPI_COMM_WORLD, &request);
            }
        }
        /* Receiver side */
        else {

            for ( i = 1; i < p; i++ ) {

               if( NumRowB[i] < 1  || RowB[i] >= N) continue;

                MPI_Recv( &A[ RowA[i] ], NumRowA[i] , MPI_FLOAT, i,0, MPI_COMM_WORLD, &status );
                MPI_Recv( &B[ RowB[i] ], NumRowB[i] , MPI_FLOAT, i,0, MPI_COMM_WORLD, &status );
            }

            
        }
       
    }

 /* Back Substitution */
   if ( my_rank == 0 ) {
  
       int row, col;  

   
       for (row = N - 1; row >= 0; row--) {
          X[row] = B[row];

        for (col = N-1; col > row; col--){
            X[row] -= A[N*row+col] * X[col];
           }
        
        X[row] /= A[N*row + row];
      }
   }
}
