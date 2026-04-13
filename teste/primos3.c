#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/sysinfo.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>

#define BLOCK_SIZE (1 << 20)  // 1MB blocks para reduzir overhead
#define TIMEOUT_SEC 59.5
#define OUTPUT_BUFFER_SIZE (1 << 24) // 16MB buffer de saída

typedef struct {
    uint64_t start;
    uint64_t end;
    uint64_t *primes;
    uint64_t count;
    uint64_t *base_primes;
    uint64_t base_count;
    int thread_id;
    int active;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
} ThreadData;

// Fast int to string (aproximadamente 3-5x mais rápido que sprintf)
static inline char* fast_utoa(uint64_t value, char *buffer) {
    static const char digits[] = "00010203040506070809"
                                 "10111213141516171819"
                                 "20212223242526272829"
                                 "30313233343536373839"
                                 "40414243444546474849"
                                 "50515253545556575859"
                                 "60616263646566676869"
                                 "70717273747576777879"
                                 "80818283848586878889"
                                 "90919293949596979899";
    
    char *ptr = buffer;
    char temp[32];
    int len = 0;
    
    // Processa de 2 em 2 dígitos
    while (value >= 100) {
        uint32_t idx = (value % 100) * 2;
        value /= 100;
        temp[len++] = digits[idx + 1];
        temp[len++] = digits[idx];
    }
    
    // Últimos 1-2 dígitos
    if (value < 10) {
        temp[len++] = '0' + value;
    } else {
        uint32_t idx = value * 2;
        temp[len++] = digits[idx + 1];
        temp[len++] = digits[idx];
    }
    
    // Reverte para buffer de saída
    while (len--) {
        *ptr++ = temp[len];
    }
    *ptr = '\0';
    
    return ptr;
}

// Sieve segmentado apenas para ímpares (wheel factorization básico)
void segmented_sieve_odd(uint64_t low, uint64_t high, uint64_t *base_primes, 
                         uint64_t base_count, uint64_t *out_primes, uint64_t *out_count) {
    if (low < 3) low = 3;
    if (low > high) {
        *out_count = 0;
        return;
    }
    
    // Ajustar para ímpar
    if ((low & 1) == 0) low++;
    
    uint64_t range = (high - low) / 2 + 1; // Apenas ímpares
    uint8_t *is_prime = malloc(range);
    if (!is_prime) {
        *out_count = 0;
        return;
    }
    
    memset(is_prime, 1, range);
    
    for (uint64_t i = 0; i < base_count; i++) {
        uint64_t p = base_primes[i];
        if (p == 2) continue; // Pular 2, estamos só em ímpares
        if (p * p > high) break;
        
        // Primeiro ímpar múltiplo de p >= low
        uint64_t start = ((low + p - 1) / p) * p;
        if ((start & 1) == 0) start += p; // Garantir ímpar
        if (start == p) start += 2 * p; // Pular o próprio primo
        
        // Marcar múltiplos ímpares: start, start+2p, start+4p...
        for (uint64_t j = start; j <= high; j += 2 * p) {
            is_prime[(j - low) / 2] = 0;
        }
    }
    
    // Contar e coletar em um único passe
    uint64_t count = 0;
    for (uint64_t i = 0; i < range; i++) {
        if (is_prime[i]) {
            out_primes[count++] = low + 2 * i;
        }
    }
    
    free(is_prime);
    *out_count = count;
}

// Worker thread pool
void* worker_thread(void* arg) {
    ThreadData *data = (ThreadData*)arg;
    
    while (1) {
        pthread_mutex_lock(&data->mutex);
        while (!data->active && data->thread_id >= 0) {
            pthread_cond_wait(&data->cond, &data->mutex);
        }
        
        if (data->thread_id < 0) { // Sinal de shutdown
            pthread_mutex_unlock(&data->mutex);
            break;
        }
        
        // Processar
        segmented_sieve_odd(data->start, data->end, data->base_primes, 
                           data->base_count, data->primes, &data->count);
        
        data->active = 0;
        pthread_mutex_unlock(&data->mutex);
    }
    
    return NULL;
}

// Gera ou estende base_primes (cache persistente)
uint64_t* get_base_primes(uint64_t needed_limit, uint64_t **cached_primes, 
                          uint64_t *cached_count, uint64_t *cached_capacity) {
    if (!*cached_primes) {
        *cached_capacity = 100000;
        *cached_primes = malloc(*cached_capacity * sizeof(uint64_t));
        *cached_count = 0;
        
        // Seed com pequenos primos
        (*cached_primes)[0] = 2;
        (*cached_primes)[1] = 3;
        (*cached_primes)[2] = 5;
        *cached_count = 3;
    }
    
    uint64_t current_max = (*cached_primes)[*cached_count - 1];
    if (current_max * current_max >= needed_limit) return *cached_primes;
    
    // Estender crivo até sqrt(needed_limit)
    uint64_t new_limit = (uint64_t)sqrt((double)needed_limit) + 100;
    if (new_limit <= current_max) return *cached_primes;
    
    uint8_t *is_prime = calloc(new_limit - current_max, sizeof(uint8_t));
    memset(is_prime, 1, new_limit - current_max);
    
    // Marcar com primos existentes
    for (uint64_t i = 0; i < *cached_count; i++) {
        uint64_t p = (*cached_primes)[i];
        uint64_t start = ((current_max + 1 + p - 1) / p) * p;
        for (uint64_t j = start; j <= new_limit; j += p) {
            if (j > current_max) is_prime[j - current_max - 1] = 0;
        }
    }
    
    // Adicionar novos primos
    for (uint64_t i = 0; i < new_limit - current_max; i++) {
        if (is_prime[i]) {
            uint64_t p = current_max + 1 + i;
            if (*cached_count >= *cached_capacity) {
                *cached_capacity *= 2;
                *cached_primes = realloc(*cached_primes, *cached_capacity * sizeof(uint64_t));
            }
            (*cached_primes)[(*cached_count)++] = p;
            
            // Marcar múltiplos do novo primo
            for (uint64_t j = p * p; j <= new_limit; j += p) {
                if (j > current_max) is_prime[j - current_max - 1] = 0;
            }
        }
    }
    
    free(is_prime);
    return *cached_primes;
}

double get_time() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// Função auxiliar para write com verificação de erro
static inline void safe_write(int fd, const void *buf, size_t count) {
    size_t written = 0;
    while (written < count) {
        ssize_t n = write(fd, (const char*)buf + written, count - written);
        if (n < 0) {
            if (errno == EINTR) continue;
            perror("write failed");
            break;
        }
        written += (size_t)n;
    }
}

int main() {
    long num_threads = sysconf(_SC_NPROCESSORS_ONLN);
    if (num_threads < 1) num_threads = 2;
    
    printf("Iniciando com %ld threads, timeout %.1fs\n", num_threads, TIMEOUT_SEC);
    
    // Pré-alocar estruturas
    uint64_t chunk_size = BLOCK_SIZE / num_threads;
    chunk_size = (chunk_size / 2) * 2;
    if (chunk_size < 2) chunk_size = 2;
    
    ThreadData *thread_data = calloc(num_threads, sizeof(ThreadData));
    pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
    
    // Pré-alocar buffers de primos para cada thread (tamanho máximo estimado)
    uint64_t max_primes_per_thread = BLOCK_SIZE / 10;
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].primes = malloc(max_primes_per_thread * sizeof(uint64_t));
        thread_data[i].thread_id = i;
        thread_data[i].active = 0;
        pthread_mutex_init(&thread_data[i].mutex, NULL);
        pthread_cond_init(&thread_data[i].cond, NULL);
        pthread_create(&threads[i], NULL, worker_thread, &thread_data[i]);
    }
    
    // Cache de primos base
    uint64_t *base_primes = NULL;
    uint64_t base_count = 0, base_capacity = 0;
    
    // Buffer de saída grande
    char *output = malloc(OUTPUT_BUFFER_SIZE);
    char *out_ptr = output;
    size_t out_remaining = OUTPUT_BUFFER_SIZE;
    
    double start = get_time();
    uint64_t current_start = 2;
    uint64_t total_primes = 0;
    int first_output = 1;
    
    while (1) {
        if (get_time() - start >= TIMEOUT_SEC) break;
        
        uint64_t current_end = current_start + BLOCK_SIZE - 1;
        
        // Obter/estender primos base
        base_primes = get_base_primes(current_end, &base_primes, &base_count, &base_capacity);
        
        // Distribuir trabalho
        uint64_t pos = current_start;
        uint64_t remaining = BLOCK_SIZE;
        
        for (int i = 0; i < num_threads && remaining > 0; i++) {
            uint64_t this_size = (i == num_threads - 1) ? remaining : chunk_size;
            this_size = (this_size / 2) * 2;
            if (this_size == 0) this_size = 2;
            
            pthread_mutex_lock(&thread_data[i].mutex);
            thread_data[i].start = pos;
            thread_data[i].end = pos + this_size - 1;
            thread_data[i].base_primes = base_primes;
            thread_data[i].base_count = base_count;
            thread_data[i].active = 1;
            pthread_cond_signal(&thread_data[i].cond);
            pthread_mutex_unlock(&thread_data[i].mutex);
            
            pos += this_size;
            remaining -= this_size;
        }
        
        // Aguardar conclusão
        int completed;
        do {
            completed = 1;
            for (int i = 0; i < num_threads; i++) {
                pthread_mutex_lock(&thread_data[i].mutex);
                if (thread_data[i].active) completed = 0;
                pthread_mutex_unlock(&thread_data[i].mutex);
            }
            if (!completed) sched_yield();
        } while (!completed && (get_time() - start) < TIMEOUT_SEC);
        
        // Coletar resultados e formatar saída
        uint64_t batch_count = 0;
        for (int i = 0; i < num_threads; i++) {
            batch_count += thread_data[i].count;
        }
        
        // Verificar se cabe no buffer, se não, flush
        if (out_remaining < batch_count * 22) {
            safe_write(STDOUT_FILENO, output, out_ptr - output);
            out_ptr = output;
            out_remaining = OUTPUT_BUFFER_SIZE;
        }
        
        // Formatar com fast_utoa
        for (int i = 0; i < num_threads; i++) {
            for (uint64_t j = 0; j < thread_data[i].count; j++) {
                if (!first_output) {
                    *out_ptr++ = ',';
                    *out_ptr++ = ' ';
                    out_remaining -= 2;
                }
                first_output = 0;
                
                char temp[32];
                char *end = fast_utoa(thread_data[i].primes[j], temp);
                size_t len = end - temp;
                
                if (out_remaining > len) {
                    memcpy(out_ptr, temp, len);
                    out_ptr += len;
                    out_remaining -= len;
                }
            }
        }
        *out_ptr++ = '\n';
        out_remaining--;
        
        total_primes += batch_count;
        current_start = current_end + 1;
    }
    
    // Flush final
    if (out_ptr > output) {
        safe_write(STDOUT_FILENO, output, out_ptr - output);
    }
    
    // Shutdown threads
    for (int i = 0; i < num_threads; i++) {
        pthread_mutex_lock(&thread_data[i].mutex);
        thread_data[i].thread_id = -1; // Sinaliza shutdown
        pthread_cond_signal(&thread_data[i].cond);
        pthread_mutex_unlock(&thread_data[i].mutex);
        pthread_join(threads[i], NULL);
        pthread_mutex_destroy(&thread_data[i].mutex);
        pthread_cond_destroy(&thread_data[i].cond);
        free(thread_data[i].primes);
    }
    
    double elapsed = get_time() - start;
    
    // Saída final
    printf("Quantidade de primos: %lu\n", total_primes);
    printf("Tempo gasto: %.2f\n", elapsed);
    
    free(thread_data);
    free(threads);
    free(base_primes);
    free(output);
    
    return 0;
}