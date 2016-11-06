#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <sys/time.h>

static unsigned long get_nth_prime_serial(unsigned long n);
static unsigned long get_nth_prime_parallel(unsigned long n);
static unsigned long nth_prime_upper(unsigned long n);
int compare (const void * a, const void * b);

int compare (const void * a, const void * b)
{
  return ( *(float*)a - *(float*)b );
}

void main() {
    int max_thread = omp_get_max_threads();
    unsigned long ns[] = {500000, 1000000, 1500000, 2000000, 4000000, 8000000, 16000000};
    int ns_length = sizeof(ns)/sizeof(unsigned long);
    float *global_exec_times = (float *)malloc(ns_length* max_thread * sizeof(float));

    //test serial implementation 
    printf("Serial\n");
    for(int i = 0; i < ns_length; ++i) {
        printf("%lu ", ns[i]);
        //Get Median of 10 times excecution
        float local_exec_times[10];
        for(int j=0; j<10; ++j) {
            struct timeval t0;
            struct timeval t1;
            gettimeofday(&t0, 0);
            get_nth_prime_serial(ns[i]);
            gettimeofday(&t1, 0);
            local_exec_times[j] = (t1.tv_sec - t0.tv_sec) * 1000.0f + (t1.tv_usec - t0.tv_usec) / 1000.0f;
        }
        qsort(local_exec_times, 10, sizeof(float), compare);
        global_exec_times[i] = local_exec_times[4];
        printf("%f\n", local_exec_times[4]);
    }
    printf("\n");

    omp_set_dynamic(0);
    for(int num_threads = 2; num_threads <= max_thread; ++num_threads) {
        omp_set_num_threads(num_threads);
        printf("Num of thread = %d\n", num_threads);
        for(int i = 0; i < ns_length; ++i) {
            printf("%lu ", ns[i]);
            //Get Median of 10 times excecution
            float local_exec_times[10];
            for(int j=0; j<10; ++j) {
                struct timeval t0;
                struct timeval t1;
                gettimeofday(&t0, 0);
                get_nth_prime_parallel(ns[i]);
                gettimeofday(&t1, 0);
                local_exec_times[j] = (t1.tv_sec - t0.tv_sec) * 1000.0f + (t1.tv_usec - t0.tv_usec) / 1000.0f;
            }
            qsort(local_exec_times, 10, sizeof(float), compare);
            global_exec_times[(num_threads-1)*ns_length+i] = local_exec_times[4];
            printf("%f\n", local_exec_times[4]);
        }
        printf("\n");
    }
    omp_set_dynamic(1);

    //write test execution times to csv file
    FILE *f = fopen("docs/output.csv", "w+");
    if (f == NULL)
    {
        printf("Cannot write output file!\n");
        exit(1);
    }

    fprintf(f, "Mode/n,Serial,");

    for(int i = 2; i <= max_thread; ++i) {
        fprintf(f, "%d Thread,", i);
    }
    fprintf(f, "\n");

    for(int i = 0; i < ns_length; ++i) {
        fprintf(f,"%lu,",ns[i]);
        for(int j = 0; j < max_thread; ++j) {
            fprintf(f, "%f,",global_exec_times[i + j*ns_length ]);
        }
        fprintf(f,"\n");

    }
}

/**
* Only check odd numbers
**/
static unsigned long get_nth_prime_serial(unsigned long n) {
    if(n == 1) {
        return 2;
    }
    char *is_prime;
    unsigned long m = nth_prime_upper(n);
    unsigned long sqrt_m = sqrt(m);
    unsigned long sieve_size = m/2+1;
    is_prime = malloc(sieve_size * sizeof(char));
    
    for(int i = 0; i <= sieve_size; ++i) {
        is_prime[i]=1;
    }

    //Only check for i < sqrt(m)
    for(unsigned long i = 3; i <= sqrt_m; i += 2) {
        if(is_prime[i/2]) {
            //start from i*i as previous multiplication should have already been marked as false
            for(unsigned long ci = i*i; ci < m; ci += 2*i) {
                is_prime[ci/2] = 0;
            }
        }
    }

    unsigned long counter = 1;
    unsigned long it = 0;
    while(counter < n) {
        ++it;
        if(is_prime[it]) {
            ++counter;
        }
    }

    free(is_prime);
    return it*2+1;
}

static unsigned long get_nth_prime_parallel(unsigned long n) {
    if(n == 1) {
        return 2;
    }

    char *is_prime;
    unsigned long m = nth_prime_upper(n);
    unsigned long sqrt_m = sqrt(m);
    unsigned long sieve_size = m/2+1;
    is_prime = malloc(sieve_size * sizeof(char));
    
    #pragma omp parallel for
    for(int i = 0; i <= sieve_size; ++i) {
        is_prime[i]=1;
    }

    #pragma omp parallel for schedule(dynamic)
    for(unsigned long i = 3; i <= sqrt_m; i += 2) {
        if(is_prime[i/2]) {
            //start from i*i as previous multiplication should have already been marked as false
            for(unsigned long ci = i*i; ci < m; ci += 2*i) {
                is_prime[ci/2] = 0;
            }
        }
    }

    unsigned long counter = 1;
    unsigned long it = 0;
    while(counter < n) {
        ++it;
        if(is_prime[it]) {
            ++counter;
        }
    }

    free(is_prime);
    return it*2+1;
}

/**
* Function to get approximation boundary of the prime search
* https://stackoverflow.com/questions/1042717/is-there-a-way-to-find-the-approximate-value-of-the-nth-prime
**/
static unsigned long nth_prime_upper(unsigned long n) {
    unsigned short primes_small[] = {0,2,3,5,7,11};
    double fn = (double) n;
    double flogn, flog2n, upper;
    if (n < 6)
        return primes_small[n];
    flogn  = log(n);
    flog2n = log(flogn);

    if (n >= 688383)    /* Dusart 2010 page 2 */
        upper = fn * (flogn + flog2n - 1.0 + ((flog2n-2.00)/flogn));
    else if (n >= 178974)    /* Dusart 2010 page 7 */
        upper = fn * (flogn + flog2n - 1.0 + ((flog2n-1.95)/flogn));
    else if (n >=  39017)    /* Dusart 1999 page 14 */
        upper = fn * (flogn + flog2n - 0.9484);
    else                    /* Modified from Robin 1983 for 6-39016 _only_ */
        upper = fn * ( flogn  +  0.6000 * flog2n );

    if (upper >= (double) ULONG_MAX) {
        /* Adjust this as needed for your type and exception method */
        if (n <= 425656284035217743UL) return 18446744073709551557UL;
        fprintf(stderr, "nth_prime_upper overflow for n = %lu\n",n);
        exit(-1);
    }

    return (unsigned long) ceil(upper);
}
