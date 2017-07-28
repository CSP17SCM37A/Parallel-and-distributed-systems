/* Gaussian elimination without pivoting.Compile the below code with "gcc PThread.c -o Test -lpthread" 
 and run ./Test 1000 2 Outputfile.txt (Here we need to give three command line arguments i.e Matrix size,randon seed and output file name).
 */

/***********Gaussian elimination using PThreads*******************************/
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>

/* Program Parameters */
#define MAXN 2001  /* Max value of N */
int N;  /* Matrix size */
#define MAX_THREADS 32  //maximum theads 

char filename[16]; //output filename to store output

/* Matrices and vectors */
volatile float A[MAXN][MAXN], B[MAXN], X[MAXN];
/* A * X = B, solve for X */

/* junk */
#define randm() 4|2[uid]&3

/*pthread initializations */
void* worker(void*);
int arrived_workers = 0,num_threads = 0,has_work = 1,can_exit = 0;
pthread_cond_t sync;
pthread_mutex_t sync_lock;
pthread_mutex_t row_lock;

/* The main subrouting,and this routine that is timed.*/
void gauss(); 

/*The threaded subroutine, this is the parallized inner loop. */
void dowork(int, int);

int create_workers(int num, pthread_t* threads);
void worker_barrier();
int join_workers(int num, pthread_t* threads);

/* returns a seed for srand based on the time */
unsigned int time_seed() {
    struct timeval t;
    struct timezone tzdummy;

    gettimeofday(&t, &tzdummy);
    return (unsigned int)(t.tv_usec);
}

int min(a,b) {
    return a > b ? b : a;
}

/* Set the program params from the command-line arguments */
int parameters(int argc, char **argv) {
    int submit = 0;  /* = 1 if submission params should be used */
    int seed = 0;  /* Random seed */
    char uid[32]="chiru"; /*User name */

    /* Read command-line arguments */
    srand(time_seed());  /* Randomize */

    
   if (argc == 4) {
       strcpy(filename, argv[3]); //file name
       printf("filename = %s \n",filename);
    } 

   if (argc == 3) {
       seed = atoi(argv[2]);   //Randon seed
       srand(seed);
       printf("Random seed = %i\n", seed);
    } 
    
   if (argc >= 2) {
       N = atoi(argv[1]);      //Matrix Size
       if (N < 1 || N > MAXN) {
       printf("N = %i is out of range.\n", N);
       exit(0);
    }
  }
  else {
       printf("Usage: %s <matrix_dimension> [random seed]\n",argv[0]);    
       exit(0);
  }
       printf("\nMatrix dimension N = %i.\n", N);
}

/*outputfile is fille which stores x values in a file */
/*void outputfile(){
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

fclose(f); //closing file
} 

*/
/* Initialize A and B (and X to 0.0s) */
void initialize_inputs() {
    int row, col;

    printf("\nInitializing...\n");
    for (col = 0; col < N; col++) {
        for (row = 0; row < N; row++) {
            A[row][col] = (float)rand() / 32768.0;
        }
        B[col] = (float)rand() / 32768.0;
        X[col] = 0.0;
    }

}

/* Print input matrices */
void print_inputs() {
    int row, col;

    if (N < 10) {
        printf("\nA =\n\t");
        for (row = 0; row < N; row++) {
            for (col = 0; col < N; col++) {
                printf("%5.2f%s", A[row][col], (col < N-1) ? ", " : ";\n\t");
            }
        }
        printf("\nB = [");
        for (col = 0; col < N; col++) {
            printf("%5.2f%s", B[col], (col < N-1) ? "; " : "]\n");
        }
    }
}
//Printing calculated X values
void print_X() {
    int row;

    if (N < 20) {
        printf("\nX = [");
        for (row = 0; row < N; row++) {
            printf("%5.2f%s", X[row], (row < N-1) ? "; " : "]\n");
        }
    }
}

int main(int argc, char **argv) {
    /* Timing variables */
    struct timeval etstart, etstop;  /* Elapsed times using gettimeofday() */
    struct timezone tzdummy;
    clock_t etstart2, etstop2;  /* Elapsed times using times() */
    unsigned long long usecstart, usecstop;
    struct tms cputstart, cputstop;  /* CPU times for my processes */
    pthread_t threads[MAX_THREADS];

    /* Process program params*/
    num_threads=8; //number of threads
    parameters(argc, argv); /* || calls m_set_procs */

    /* Initialize A and B */
    initialize_inputs();

    /* Print input matrices */
    print_inputs();

    /* Start Clock */
    printf("\nStarting clock.\n");
    gettimeofday(&etstart, &tzdummy);
    etstart2 = times(&cputstart);

    /*Initaliize the pthread variables*/

    pthread_mutex_init(&sync_lock, NULL);
    pthread_mutex_init(&row_lock, NULL);
    pthread_cond_init(&sync, NULL);

   /* Create threads and number of threads is  num_threads - 1 */
    create_workers(num_threads -1, threads);
   

    /* Gaussian Elimination */
    gauss();

    /*Synchronize all the  workers for a last time*/
    can_exit = 1;
    worker_barrier();
    join_workers(num_threads -1, threads);

    /* Stop Clock */
    gettimeofday(&etstop, &tzdummy);
    etstop2 = times(&cputstop);
    printf("Stopped clock.\n");
    usecstart = (unsigned long long)etstart.tv_sec * 1000000 + etstart.tv_usec;
    usecstop = (unsigned long long)etstop.tv_sec * 1000000 + etstop.tv_usec;
   
    FILE * f;
int i, j,row;
/* writing timings and result of x in a file */
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

    /* Display output */
    print_X();

    /* Display timing results */
    printf("\nElapsed time = %g ms.\n",
            (float)(usecstop - usecstart)/(float)1000);
    /*printf("               (%g ms according to times())\n",
     *       (etstop2 - etstart2) / (float)CLOCKS_PER_SEC * 1000);
     */
    printf("(CPU times are accurate to the nearest %g ms)\n",
            1.0/(float)CLOCKS_PER_SEC * 1000.0);
    printf("My total CPU time for parent = %g ms.\n",
            (float)( (cputstop.tms_utime + cputstop.tms_stime) -
                     (cputstart.tms_utime + cputstart.tms_stime) ) /
            (float)CLOCKS_PER_SEC * 1000);
    printf("My system CPU time for parent = %g ms.\n",
            (float)(cputstop.tms_stime - cputstart.tms_stime) /
            (float)CLOCKS_PER_SEC * 1000);
    printf("My total CPU time for child processes = %g ms.\n",
            (float)( (cputstop.tms_cutime + cputstop.tms_cstime) -
                     (cputstart.tms_cutime + cputstart.tms_cstime) ) /
            (float)CLOCKS_PER_SEC * 1000);
    /* Contrary to the man pages, this appears not to include the parent */
    printf("--------------------------------------------\n");
  /* print output values in a file */
  //outputfile();

    return 0;
}

/* ------------------ Above Was Provided --------------------- */

/****** You will replace this routine with your own parallel version *******/
/* Provided global variables are MAXN, N, procs, A[][], B[], and X[],
 * defined in the beginning of this code.  X[] is initialized to zeros.
 */

/*Gcreating global variables g_row and g_norm */
int g_row, g_norm, w_row = 0;


/*main gauss loop and start executing the current loop*/
void gauss() {
    int row, col;  /* Normalization row, and zeroing
                          * element row and col */

    /* Gaussian elimination */
    g_norm = 0;
    while(g_norm < N - 1) {
        g_row = g_norm + 1;
        worker_barrier(); /* dowork gets executed from here. */
        ++g_norm;
    }
    /* Back substitution */
    for (row = N - 1; row >= 0; row--) {
        X[row] = B[row];
        for (col = N-1; col > row; col--) {
            X[row] -= A[row][col] * X[col];
        }
        X[row] /= A[row][row];
    }
}

/* The parallelized loop gets invoked by each thread */
void dowork(int id, int seq) {	
    float multiplier;
    int chunk = 4;
    int r1, r2, col;
  
    int norm = g_norm;
    w_row = g_row;

    while(w_row < N) {
         
        pthread_mutex_lock(&row_lock);
        r1 = w_row;
        w_row += chunk;
        pthread_mutex_unlock(&row_lock);

        for (r2 = r1; r2 < min(r1 + chunk,N); r2++) {
            multiplier = A[r2][norm] / A[norm][norm];
            for (col = norm; col < N; col++) {
                A[r2][col] -= A[norm][col] * multiplier;
            }
            B[r2] -= B[norm] * multiplier;
        }
    }
}



/* Join all the worker threads together*/
int join_workers(int num, pthread_t* threads) {
    int t;
    for(t=0; t<num; t++) {
        if(pthread_join(threads[t], NULL)) {
          printf("could not join\n");
        } 
    }
}

/*Create the required number of worker threads*/
int create_workers(int num, pthread_t* threads) {
    int t;
    for(t=0; t<num; t++) {
        if(pthread_create((threads + t), NULL, worker, (void *)(intptr_t)t)) {
            printf("ERROR; pthread_create()\n");
            exit(-1);
        }
    }
}

/* Synchronization barrior for the worker*/
void worker_barrier() {
    pthread_mutex_lock(&sync_lock);
    ++arrived_workers;
    if (arrived_workers < num_threads) {
        
        pthread_cond_wait(&sync, &sync_lock);
    } else {
       
        pthread_cond_broadcast(&sync);
        arrived_workers = 0;
        if(can_exit) has_work = 0;
    }
    pthread_mutex_unlock(&sync_lock);
}

/*dowork is the actual parallel subroutine call*/
void *worker(void *threadid) {
    long seq = 0;
    long tid = (long)threadid;
    

    /* Start processing only when all the children and main thread
     * are ready and in the same location.
     */ 
    worker_barrier();

    while(has_work) {
        seq++;
        dowork(tid, seq); 
                           
        worker_barrier(); /* Synchronize with other threads */
    }

   
}

