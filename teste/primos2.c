#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stddef.h>   
#include <sched.h>


#define BLOCK_SIZE      (1ULL << 20)
#define TIMEOUT_SEC     59.5
#define MAX_THREADS     256

static const uint8_t WHEEL30_OFFSETS[8] = { 1,  7, 11, 13, 17, 19, 23, 29};
static const uint8_t WHEEL30_GAPS[8]    = { 6,  4,  2,  4,  2,  4,  6,  2};

static const char DIGIT_TABLE[200] = {
    '0','0','0','1','0','2','0','3','0','4','0','5','0','6','0','7','0','8','0','9',
    '1','0','1','1','1','2','1','3','1','4','1','5','1','6','1','7','1','8','1','9',
    '2','0','2','1','2','2','2','3','2','4','2','5','2','6','2','7','2','8','2','9',
    '3','0','3','1','3','2','3','3','3','4','3','5','3','6','3','7','3','8','3','9',
    '4','0','4','1','4','2','4','3','4','4','4','5','4','6','4','7','4','8','4','9',
    '5','0','5','1','5','2','5','3','5','4','5','5','5','6','5','7','5','8','5','9',
    '6','0','6','1','6','2','6','3','6','4','6','5','6','6','6','7','6','8','6','9',
    '7','0','7','1','7','2','7','3','7','4','7','5','7','6','7','7','7','8','7','9',
    '8','0','8','1','8','2','8','3','8','4','8','5','8','6','8','7','8','8','8','9',
    '9','0','9','1','9','2','9','3','9','4','9','5','9','6','9','7','9','8','9','9',
};

static inline char* u64_to_str(char *buf, uint64_t n) {
    char tmp[20];
    char *p = tmp + 20;

    while (n >= 100) {
        uint64_t q  = n / 100;
        uint32_t r  = (uint32_t)(n - q * 100) * 2;
        n = q;
        *--p = DIGIT_TABLE[r + 1];
        *--p = DIGIT_TABLE[r];
    }
    if (n >= 10) {
        uint32_t r = (uint32_t)n * 2;
        *--p = DIGIT_TABLE[r + 1];
        *--p = DIGIT_TABLE[r];
    } else {
        *--p = (char)('0' + n);
    }

    ptrdiff_t len = (tmp + 20) - p;
    memcpy(buf, p, (size_t)len);
    return buf + len;
}


static uint64_t *g_base_primes = NULL;
static uint64_t  g_base_count  = 0;
static uint64_t  g_base_limit  = 0;

static void ensure_base_primes(uint64_t new_limit) {
    uint64_t need = (uint64_t)sqrt((double)new_limit) + 2;
    if (need <= g_base_limit) return;
    g_base_limit = need;

    free(g_base_primes);


    uint8_t *sieve = malloc(need + 1);
    if (!sieve) { g_base_count = 0; g_base_primes = NULL; return; }

    memset(sieve, 1, need + 1);
    sieve[0] = sieve[1] = 0;
    for (uint64_t i = 3; i * i <= need; i += 2)
        if (sieve[i])
            for (uint64_t j = i * i; j <= need; j += 2 * i)
                sieve[j] = 0;

    uint64_t cnt = (need >= 2) ? 1 : 0;
    for (uint64_t i = 3; i <= need; i += 2)
        if (sieve[i]) cnt++;

    g_base_primes = malloc(cnt * sizeof(uint64_t));
    uint64_t idx  = 0;
    if (need >= 2) g_base_primes[idx++] = 2;
    for (uint64_t i = 3; i <= need; i += 2)
        if (sieve[i]) g_base_primes[idx++] = i;

    g_base_count = cnt;
    free(sieve);
}


typedef struct {
    uint64_t  start;
    uint64_t  end;
    uint64_t *base_primes;
    uint64_t  base_count;

  
    uint8_t  *sieve_buf;
    uint64_t  sieve_buf_cap;

    uint64_t  count;
    char     *out_buf;
    size_t    out_bytes;
} ThreadData;


static void run_sieve(ThreadData *data) {
    uint64_t low   = data->start;
    uint64_t high  = data->end;
    uint64_t range = high - low + 1;


    if (range > data->sieve_buf_cap) {
        free(data->sieve_buf);
        data->sieve_buf     = malloc(range);
        data->sieve_buf_cap = range;
    }
    uint8_t *is_prime = data->sieve_buf;


    memset(is_prime, 1, range);
    if (low == 0) is_prime[0] = 0;
    if (low <= 1 && high >= 1) is_prime[1 - low] = 0;


    for (uint64_t i = 0; i < data->base_count; i++) {
        uint64_t p = data->base_primes[i];
        if (p * p > high) break;

        uint64_t start = ((low + p - 1) / p) * p;
        if (start < p * p) start = p * p;
        if (start > high)  continue;

        uint64_t step = (p == 2) ? p : 2 * p;
        if (p > 2 && (start & 1) == 0) start += p;

        for (uint64_t j = start; j <= high; j += step) {

            __builtin_prefetch(&is_prime[j - low + step], 1, 0);
            is_prime[j - low] = 0;
        }
    }

    char    *ptr   = data->out_buf;
    uint64_t count = 0;

    for (uint64_t special = 2; special <= 5; special++) {
        if (special < low || special > high) continue;
        if (!is_prime[special - low])        continue;
        ptr = u64_to_str(ptr, special);
        *ptr++ = ','; *ptr++ = ' ';
        count++;
    }

    uint64_t base30 = (low / 30) * 30;     
    if (base30 < 6) base30 = 6;            

    for (uint64_t w = base30; w <= high; w += 30) {
        for (int k = 0; k < 8; k++) {
            uint64_t n = w + WHEEL30_OFFSETS[k];
            if (n < low || n < 7) continue;   
            if (n > high)         goto done;   
            if (!is_prime[n - low]) continue;
            ptr = u64_to_str(ptr, n);
            *ptr++ = ','; *ptr++ = ' ';
            count++;
        }
        (void)WHEEL30_GAPS;   
                               
                                
                               
                                
                                
    }
done:
    data->count     = count;
    data->out_bytes = (size_t)(ptr - data->out_buf);
}

typedef struct {
    pthread_t       tid;
    pthread_mutex_t mtx;
    pthread_cond_t  cv_work;
    pthread_cond_t  cv_done;
    ThreadData     *data;
    int             ready;
    int             done;
    int             shutdown;
} Worker;

static Worker g_workers[MAX_THREADS];
static int    g_num_threads = 0;

static void* worker_loop(void *arg) {
    Worker *w  = (Worker*)arg;
    int     id = (int)(w - g_workers);


    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(id % (int)sysconf(_SC_NPROCESSORS_ONLN), &cpuset);
    pthread_setaffinity_np(pthread_self(), sizeof(cpuset), &cpuset);

    pthread_mutex_lock(&w->mtx);
    for (;;) {
        while (!w->ready && !w->shutdown)
            pthread_cond_wait(&w->cv_work, &w->mtx);
        if (w->shutdown) { pthread_mutex_unlock(&w->mtx); break; }

        w->ready = 0;
        pthread_mutex_unlock(&w->mtx);

        run_sieve(w->data);

        pthread_mutex_lock(&w->mtx);
        w->done = 1;
        pthread_cond_signal(&w->cv_done);
    }
    return NULL;
}

static void pool_init(int n) {
    g_num_threads = n;
    for (int i = 0; i < n; i++) {
        Worker *w = &g_workers[i];
        pthread_mutex_init(&w->mtx,     NULL);
        pthread_cond_init (&w->cv_work, NULL);
        pthread_cond_init (&w->cv_done, NULL);
        w->ready = w->done = w->shutdown = 0;
        w->data  = NULL;
        pthread_create(&w->tid, NULL, worker_loop, w);
    }
}

static void pool_dispatch(int id, ThreadData *data) {
    Worker *w = &g_workers[id];
    pthread_mutex_lock(&w->mtx);
    w->data  = data;
    w->ready = 1;
    w->done  = 0;
    pthread_cond_signal(&w->cv_work);
    pthread_mutex_unlock(&w->mtx);
}

static void pool_wait(int id) {
    Worker *w = &g_workers[id];
    pthread_mutex_lock(&w->mtx);
    while (!w->done)
        pthread_cond_wait(&w->cv_done, &w->mtx);
    pthread_mutex_unlock(&w->mtx);
}

static void pool_shutdown(void) {
    for (int i = 0; i < g_num_threads; i++) {
        Worker *w = &g_workers[i];
        pthread_mutex_lock(&w->mtx);
        w->shutdown = 1;
        pthread_cond_signal(&w->cv_work);
        pthread_mutex_unlock(&w->mtx);
        pthread_join(w->tid, NULL);
    }
}


int main(void) {
    
    setvbuf(stdout, NULL, _IOFBF, 1 << 20);

    int num_threads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    if (num_threads < 1)        num_threads = 1;
    if (num_threads > MAX_THREADS) num_threads = MAX_THREADS;

    uint64_t chunk_size = BLOCK_SIZE / (uint64_t)num_threads;
    if (chunk_size < 30) chunk_size = 30;
    chunk_size = (chunk_size / 30) * 30;   


    size_t per_thread_out = (size_t)(chunk_size + 32) * 22;
    char  *global_out_buf  = malloc(per_thread_out * (size_t)num_threads);

    ThreadData *tdata = calloc((size_t)num_threads, sizeof(ThreadData));
    for (int i = 0; i < num_threads; i++) {
        tdata[i].sieve_buf     = malloc(chunk_size + 64);
        tdata[i].sieve_buf_cap = chunk_size + 64;
        tdata[i].out_buf       = global_out_buf + (size_t)i * per_thread_out;
    }

    pool_init(num_threads);

    struct timespec ts_start, ts_now;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    uint64_t current_start = 2;
    uint64_t total_primes  = 0;

    for (;;) {
        clock_gettime(CLOCK_MONOTONIC, &ts_now);
        double elapsed = (double)(ts_now.tv_sec  - ts_start.tv_sec)
                       + (double)(ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;
        if (elapsed >= TIMEOUT_SEC) break;

        uint64_t current_end = current_start + BLOCK_SIZE - 1;
        ensure_base_primes(current_end);  

        
        uint64_t pos = current_start;
        for (int i = 0; i < num_threads; i++) {
            uint64_t seg_end = (i == num_threads - 1)
                               ? current_end
                               : pos + chunk_size - 1;
            if (seg_end > current_end) seg_end = current_end;

            tdata[i].start       = pos;
            tdata[i].end         = seg_end;
            tdata[i].base_primes = g_base_primes;
            tdata[i].base_count  = g_base_count;
            tdata[i].count       = 0;
            tdata[i].out_bytes   = 0;

            pool_dispatch(i, &tdata[i]);
            pos = seg_end + 1;
        }

        for (int i = 0; i < num_threads; i++)
            pool_wait(i);

     
        for (int i = 0; i < num_threads; i++) {
            if (tdata[i].out_bytes > 2) {
              
                size_t len = tdata[i].out_bytes - 2;
                fwrite(tdata[i].out_buf, 1, len, stdout);
             
                if (i < num_threads - 1 && tdata[i + 1].out_bytes > 2)
                    fwrite(", ", 1, 2, stdout);
            }
            total_primes += tdata[i].count;
        }
        fputc('\n', stdout);

        current_start = current_end + 1;
    }
    char  final_buf[32];
    char *end_ptr = u64_to_str(final_buf, total_primes);
    *end_ptr++ = '\n';
    fwrite(final_buf, 1, (size_t)(end_ptr - final_buf), stdout);
    fflush(stdout);

    pool_shutdown();
    for (int i = 0; i < num_threads; i++)
        free(tdata[i].sieve_buf);
    free(tdata);
    free(global_out_buf);
    free(g_base_primes);
    return 0;
}