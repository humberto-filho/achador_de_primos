#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/sysinfo.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

#define BLOCK_SIZE 16384  // 2^14 mantido constante
#define TIMEOUT_SEC 59.5

typedef struct {
    uint64_t start;
    uint64_t end;
    uint64_t *primes;
    uint64_t count;
    uint64_t *base_primes;
    uint64_t base_count;
} ThreadData;

void* segmented_sieve_thread(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    uint64_t low = data->start;
    uint64_t high = data->end;
    uint64_t range = high - low + 1;
    
    if (range == 0 || low > high) {
        data->count = 0;
        data->primes = NULL;
        pthread_exit(NULL);
    }
    
    uint8_t *is_prime = calloc(range, sizeof(uint8_t));
    if (!is_prime) {
        data->count = 0;
        pthread_exit(NULL);
    }
    
    memset(is_prime, 1, range);
    if (low == 1) is_prime[0] = 0;
    
    for (uint64_t i = 0; i < data->base_count; i++) {
        uint64_t p = data->base_primes[i];
        if (p * p > high) break;
        
        uint64_t start = ((low + p - 1) / p) * p;
        if (start < p * p) start = p * p;
        
        for (uint64_t j = start; j <= high; j += p) {
            is_prime[j - low] = 0;
        }
    }
    
    data->count = 0;
    for (uint64_t i = 0; i < range; i++) {
        if (is_prime[i]) data->count++;
    }
    
    if (data->count > 0) {
        data->primes = malloc(data->count * sizeof(uint64_t));
        uint64_t idx = 0;
        for (uint64_t i = 0; i < range; i++) {
            if (is_prime[i]) {
                data->primes[idx++] = low + i;
            }
        }
    } else {
        data->primes = NULL;
    }
    
    free(is_prime);
    pthread_exit(NULL);
}

uint64_t* generate_base_primes(uint64_t limit, uint64_t *count) {
    if (limit < 2) {
        *count = 0;
        return NULL;
    }
    
    uint64_t n = (uint64_t)sqrt((double)limit) + 1;
    uint8_t *is_prime = calloc(n + 1, sizeof(uint8_t));
    if (!is_prime) {
        *count = 0;
        return NULL;
    }
    
    memset(is_prime, 1, n + 1);
    is_prime[0] = is_prime[1] = 0;
    
    for (uint64_t i = 2; i * i <= n; i++) {
        if (is_prime[i]) {
            for (uint64_t j = i * i; j <= n; j += i) {
                is_prime[j] = 0;
            }
        }
    }
    
    *count = 0;
    for (uint64_t i = 2; i <= n; i++) {
        if (is_prime[i]) (*count)++;
    }
    
    uint64_t *primes = malloc(*count * sizeof(uint64_t));
    if (primes) {
        uint64_t idx = 0;
        for (uint64_t i = 2; i <= n; i++) {
            if (is_prime[i]) primes[idx++] = i;
        }
    }
    
    free(is_prime);
    return primes;
}

int main() {
    long num_threads = sysconf(_SC_NPROCESSORS_ONLN);
    if (num_threads < 1) num_threads = 1;
    
    uint64_t chunk_size = BLOCK_SIZE / num_threads;
    chunk_size = (chunk_size / 2) * 2;
    if (chunk_size < 2) chunk_size = 2;
    
    clock_t start_time = clock();
    uint64_t current_start = 2;
    uint64_t total_primes = 0;
    
    while (1) {
        double elapsed = (double)(clock() - start_time) / CLOCKS_PER_SEC;
        if (elapsed >= TIMEOUT_SEC) break;
        
        uint64_t current_end = current_start + BLOCK_SIZE - 1;
        
        uint64_t base_count;
        uint64_t *base_primes = generate_base_primes(current_end, &base_count);
        
        pthread_t threads[num_threads];
        ThreadData thread_data[num_threads];
        
        uint64_t pos = current_start;
        uint64_t remaining = BLOCK_SIZE;
        
        for (int i = 0; i < num_threads; i++) {
            thread_data[i].start = pos;
            
            uint64_t this_size = (i == num_threads - 1) ? remaining : chunk_size;
            this_size = (this_size / 2) * 2;
            if (this_size == 0) this_size = 2;
            
            thread_data[i].end = pos + this_size - 1;
            thread_data[i].base_primes = base_primes;
            thread_data[i].base_count = base_count;
            thread_data[i].count = 0;
            thread_data[i].primes = NULL;
            
            pthread_create(&threads[i], NULL, segmented_sieve_thread, &thread_data[i]);
            
            pos += this_size;
            remaining -= this_size;
        }
        
        for (int i = 0; i < num_threads; i++) {
            pthread_join(threads[i], NULL);
        }
        
        uint64_t batch_count = 0;
        for (int i = 0; i < num_threads; i++) {
            batch_count += thread_data[i].count;
        }
        
        size_t buffer_size = batch_count * 16 + 10;
        char *buffer = malloc(buffer_size);
        char *ptr = buffer;
        
        for (int i = 0; i < num_threads; i++) {
            for (uint64_t j = 0; j < thread_data[i].count; j++) {
                int written = sprintf(ptr, "%lu", thread_data[i].primes[j]);
                ptr += written;
                
                int is_last_global = (i == num_threads - 1 && j == thread_data[i].count - 1);
                if (!is_last_global) {
                    *ptr++ = ',';
                    *ptr++ = ' ';
                }
            }
            free(thread_data[i].primes);
        }
        *ptr = '\0';
        
        printf("%s\n", buffer);
        free(buffer);
        free(base_primes);
        
        total_primes += batch_count;
        current_start = current_end + 1;
        
        elapsed = (double)(clock() - start_time) / CLOCKS_PER_SEC;
        if (elapsed >= TIMEOUT_SEC) break;
    }
    
    // ===== ALTERAÇÕES NO FINAL =====
    double total_elapsed = (double)(clock() - start_time) / CLOCKS_PER_SEC;
    
    printf("Quantidade: %lu\n", total_primes);
    printf("Tempo elapsado: %.2f segundos\n", total_elapsed);
    // ================================
    
    return 0;
}