/* Matrix normalization.
* Compile with nvcc normcuda.cu -o abcde.out and run with ./abcde.out 2000 1 2000 2(matrix size,number of blocks,number of threads in each thread ,randon seed)
*/

/* ****** ADD YOUR CODE AT THE END OF THIS FILE. ******
* You need not submit the provided code.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>

/* Program Parameters */
#define MAXN 8000  /* Max value of N */
int N;  /* Matrix size */

/* Matrices */
float A[MAXN*MAXN], B[MAXN*MAXN];

/* junk */
#define randm() 4|2[uid]&3
/* number of blocks and threads */
int numBlocks,numThreadsPerBlock; 
/* returns a seed for srand based on the time */
unsigned int time_seed() {
    struct timeval t;
    struct timezone tzdummy;

    gettimeofday(&t, &tzdummy);
    return (unsigned int)(t.tv_usec);
}

/* Set the program parameters from the command-line arguments */
void parameters(int argc, char **argv) {
    int seed = 0;  /* Random seed */
    //char uid[32]; /*User name */

    /* Read command-line arguments */
    srand(time_seed());  /* Randomize */

    if (argc == 5) {
        seed = atoi(argv[4]);
        srand(seed);
        printf("Random Seed = %i\n", seed);
    }
    if (argc >= 4) {
        numThreadsPerBlock = atoi(argv[3]);
        srand(seed);
        printf("Number of Threads Per Block = %i\n", numThreadsPerBlock);

        numBlocks = atoi(argv[2]);
        srand(seed);
        printf("Number of Blocks = %i\n", numBlocks);

        N = atoi(argv[1]);
        if (N < 1 || N > MAXN) {
            printf("N = %i is out of range.\n", N);
            exit(0);
        }
    }
    else {
        printf("Usage: %s <matrixDimension> <numBlocks> <numThreadsPerBlock> [randomSeed]\n",
        argv[0]);
        exit(0);
    }

    /* Print parameters */
    printf("\nMatrix dimension N = %i.\n", N);
}

/* Initialize A and B*/
void initialize_inputs() {
    int row, col;

    printf("\nInitializing...\n");
    for (col = 0; col < N; col++) {
        for (row = 0; row < N; row++) {
            A[col*N+row] = (float)rand() / 32768.0;
            B[col*N+row] = 0.0;
        }
    }

}

/* Print input matrices */
void print_inputs() {
    int row, col;

    if (N < 10) {
        printf("\nA =\n\t");
        for (row = 0; row < N; row++) {
            for (col = 0; col < N; col++) {
                printf("%5.2f%s", A[row*N+col], (col < N-1) ? ", " : ";\n\t");
            }
        }
    }
}

void print_B() {
    int row, col;

    if (N < 10) {
        printf("\nB =\n\t");
        for (row = 0; row < N; row++) {
            for (col = 0; col < N; col++) {
                printf("%1.10f%s", B[row*N+col], (col < N-1) ? ", " : ";\n\t");
            }
        }
    }
}

                                               
__global__ void normCalc (float *d_A, float *d_B, int n);

int main(int argc, char **argv) {
    /* Timing variables */
    struct timeval etstart, etstop;  /* Elapsed times using gettimeofday() */
    struct timezone tzdummy;
    //clock_t etstart2, etstop2;  /* Elapsed times using times() */
    unsigned long long usecstart, usecstop;
    struct tms cputstart, cputstop;  /* CPU times for my processes */

    float elapsed=0;
    cudaEvent_t start, stop; /* Elapsed times of gpu */

    /* Process program parameters */
    parameters(argc, argv);

    /* Initialize A and B */
    initialize_inputs();

    /* Print input matrices */
    print_inputs();

    printf("Computing in Parallel\n");

    
    float *d_A, *d_B;

    /* Start Clock */
    printf("\nStarting clock.\n");
    gettimeofday(&etstart, &tzdummy);
    times(&cputstart);

    cudaEventCreate(&start);    //creating start
    cudaEventCreate(&stop);     //creating stop
    cudaEventRecord(start, 0);  // start is zero initially
    /*Allocation */
    cudaMalloc((void **) &d_A, sizeof(float)*N*N);
    cudaMalloc((void **) &d_B, sizeof(float)*N*N);
    

    cudaMemcpy(d_A, A, sizeof(float)*N*N, cudaMemcpyHostToDevice);
       /* normalization */
    normCalc<<<numBlocks,numThreadsPerBlock>>>(d_A, d_B, N);

    cudaMemcpy(B, (d_B), sizeof(float)*N*N, cudaMemcpyDeviceToHost);
        /* Stop Clock */
    gettimeofday(&etstop, &tzdummy);
    times(&cputstop);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize (stop);

    cudaEventElapsedTime(&elapsed, start, stop) ; //elapsed time for gpu

    printf("Stopped clock.\n");
    usecstart = (unsigned long long)etstart.tv_sec * 1000000 + etstart.tv_usec;
    usecstop = (unsigned long long)etstop.tv_sec * 1000000 + etstop.tv_usec;

   
    cudaEventDestroy(start); //destroying start
    cudaEventDestroy(stop);  //destroying stop

    /* Display output */
    print_B();
    /* deallocating memory*/
    cudaFree(d_A);
    cudaFree(d_B);

    /* Display timing results */
    

    printf("\nThe elapsed time in gpu was %.2f ms\n", elapsed);
    printf("\nElapsed time = %g ms.\n",
    (float)(usecstop - usecstart)/(float)1000);

    printf("\n(CPU times are accurate to the nearest %g ms)\n",
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

    exit(0);
}





__global__ void normCalc (float *d_A, float *d_B, int n) {
    
   int col = blockIdx.x * blockDim.x + threadIdx.x;
 
 __shared__ int row;
 __shared__ float m, s;
    if (col < n){
        m = 0.0;
        for (row=0; row < n; row++)
            m += d_A[col*n+row];

        m /= (float) n;
        __syncthreads();

        s = 0.0;
        for (row=0; row < n; row++)
            s += powf(d_A[col*n+row] - m, (float)2.0);
        
        s /= (float) n;
        __syncthreads();

        s = sqrt(s);
        for (row=0; row < n; row++) {
            if (s ==(float) 0.0)
                d_B[row*n+col] = (float)0.0;
            else
                d_B[row*n+col] = (d_A[col*n+row] - m) / s;
        }
    }
}

