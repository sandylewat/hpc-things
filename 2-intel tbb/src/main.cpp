#include <iostream>
#include "tbb/tbb.h"
#include <iomanip>
#include <time.h>
#include <math.h>
using std::setw;
using namespace tbb;
using namespace std;

int get_nth_prime_serial(int n);
int nth_prime_upper(int n);
int get_nth_prime_parallel(int n, int threads_size);

int main(int argc, const char * argv[]) {
    tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic);
    int max_thread = 6;
    int ns[] = {500000, 1000000, 1500000, 2000000, 4000000, 8000000, 16000000, 32000000};
    int ns_length = sizeof(ns)/sizeof(int);
    float global_exec_times[max_thread+1][ns_length];

    printf("Serial\n");
    for(int i = 0; i < ns_length; ++i) {
        printf("%d ", ns[i]);
        //Get Median of 10 times excecution
        float local_exec_times[10];
        for(int j=0; j<10; ++j) {
            tick_count t0 = tick_count::now();
            get_nth_prime_serial(ns[i]);
            tick_count t1 = tick_count::now();
            local_exec_times[j] = (t1-t0).seconds();
        }
        sort(local_exec_times, local_exec_times+10);
        global_exec_times[0][i] = local_exec_times[4];
        printf("%f\n", local_exec_times[4]);
    }
    printf("\n");

    for(int num_threads = 2; num_threads <= max_thread; ++num_threads) {
        printf("Num of thread = %d\n", num_threads);
        for(int i = 0; i < ns_length; ++i) {
            printf("%d ", ns[i]);
            //Get Median of 10 times excecution
            float local_exec_times[10];
            for(int j=0; j<10; ++j) {
                tick_count t0 = tick_count::now();
                get_nth_prime_parallel(ns[i], num_threads);
                tick_count t1 = tick_count::now();
                local_exec_times[j] = (t1-t0).seconds();
            }
            sort(local_exec_times, local_exec_times + 10);
            global_exec_times[num_threads-1][i] = local_exec_times[4];
            printf("%f\n", local_exec_times[4]);
        }
        printf("\n");
    }

    printf("Num of thread = Auto\n");
    for(int i = 0; i < ns_length; ++i) {
        printf("%d ", ns[i]);
        //Get Median of 10 times excecution
        float local_exec_times[10];
        for(int j=0; j<10; ++j) {
            tick_count t0 = tick_count::now();
            get_nth_prime_parallel(ns[i], 0);
            tick_count t1 = tick_count::now();
            local_exec_times[j] = (t1-t0).seconds();
        }
        sort(local_exec_times, local_exec_times + 10);
        global_exec_times[max_thread][i] = local_exec_times[4];
        printf("%f\n", local_exec_times[4]);
    }
    printf("\n");    

    //write test execution times to csv file
    FILE *f = fopen("doc/output.csv", "w+");
    if (f == NULL)
    {
        printf("Cannot write output file!\n");
        exit(1);
    }

    fprintf(f, "N\\Mode,Serial,");

    for(int i = 2; i <= max_thread; ++i) {
        fprintf(f, "%d Thread,", i);
    }
    fprintf(f, "Auto,\n");

    for(int i = 0; i < ns_length; ++i) {
        fprintf(f,"%d,",ns[i]);
        for(int j = 0; j <= max_thread; ++j) {
            fprintf(f, "%f,",global_exec_times[j][i]);
        }
        fprintf(f,"\n");
    }
    // cout << get_nth_prime_parallel(ns[2],2) << endl;
    return 0;
}

/**
* Function to get approximation boundary of the prime search
* https://stackoverflow.com/questions/1042717/is-there-a-way-to-find-the-approximate-value-of-the-nth-prime
**/
int nth_prime_upper(int n) {
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

    return (int) ceil(upper);
}


int get_nth_prime_serial(int n) {
    if(n == 1) {
        return 2;
    }
    bool *is_prime;
    int m = nth_prime_upper(n);
    int sqrt_m = sqrt(m);
    int sieve_size = m/2+1;
    is_prime = new bool[sieve_size];
    
    for(int i = 0; i <= sieve_size; ++i) {
        is_prime[i]=1;
    }

    //Only check for i < sqrt(m)
    for(int i = 3; i <= sqrt_m; i += 2) {
        if(is_prime[i/2]) {
            //start from i*i as previous multiplication should have already been marked as false
            for(int ci = i*i; ci < m; ci += 2*i) {
                is_prime[ci/2] = false;
            }
        }
    }

    int counter = 1;
    int it = 0;
    while(counter < n) {
        ++it;
        if(is_prime[it]) {
            ++counter;
        }
    }

    return it*2+1;
}

int get_nth_prime_parallel(int n, int threads_size) {
    if (threads_size == 0) {
        tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic);
    }
    else {
        tbb::task_scheduler_init init(threads_size);
    }
    
    if(n == 1) {
        return 2;
    }
    int m = nth_prime_upper(n);
    int sqrt_m = sqrt(m);
    int sieve_size = m/2+1;
    concurrent_vector<bool> is_prime(sieve_size, true);

    //Only check for i < sqrt(m)
    //for(int i = 3; i <= sqrt_m; i += 2) {
    parallel_for(3, sqrt_m+1, 2, [&](int i) {
        if(is_prime[i/2]) {
            //start from i*i as previous multiplication should have already been marked as false
            for(int ci=i*i; i<m;i+=2*i) {
            //parallel_for(i*i, m, 2*i, [&](int ci) {
                is_prime[ci/2] = false;
            }
        }
    });

    int counter = 1;
    int it = 0;
    while(counter < n) {
        ++it;
        if(is_prime[it]) {
            ++counter;
        }
    }

    return it*2+1;
}