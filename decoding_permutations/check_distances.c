/* File to check whether the mapping is still a DPM as conjectured at the
 * end of the paper. Confirmed true till q <= 8, n <= 6.
 */
/* Compile as
 * gcc -O3 -march=native -ffast-math -funroll-loops check_distances.c -lm -o check_distances
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

enum { FALSE, TRUE };

int    check_distances          (int n, int Q, int q,
                                 int *index, int *store_moduli);
void   qary_to_permutations     (int *vector, int *Pvector,
                                 int n, int N, int q, int Q,
                                 int *index, int *store_moduli);
static inline void generate_next_qary_vector(int *vector, int n, int q);
static inline int  hamming_distance         (int *x, int *y, int n);
static inline void print_vector             (int *vector, int n);
static inline void set_qary_vector          (int *vector, int *values, int n);

int main(int argc, char **argv)
{
    int i, j, n, N, q, Q, ret_value;
    int *index, *store_moduli;
    // n = length of q-ary vector
    // N = length of permutation vector
    // q = alphabet size in F_q
    // Q = number of permutation symbols required in each step
    // ret_value = TRUE/FALSE, holds whether the test succeeded
    // index = value of moduli's. (q, q+Q, q+2Q, ..., q+(n-1)Q)
    // store_moduli = for an N length vector it stores the value i%index[j]
    //                in (n*i+j)-th position, i=0,...,N-1, j=0,...,n-1.

    switch (argc)
    {
        case 2:
            {
                q = atoi(argv[1]);
                n = q+1;
                break;
            }
        case 3:
            {
                q = atoi(argv[1]);
                n = atoi(argv[2]);
                break;
            }
        default:
            {
                printf ("Usage: %s <q> [<n>]\nDefault n = q+1\n\n", argv[0]);
                return 1;
            }
    }

    Q = (int)ceil(log(q)/log(2));           // number of extra symbols used
    N = q + (n-1)*Q;                        // length of permutation vector

    /* Generate the moduli indices */
    index = (int *)malloc(sizeof(int) * n);
    index[0] = q;
    for (i = 0, j = 0; i < N; i++)
        if (i == index[j])
        {
            j++;
            index[j] = index[j-1] + Q;
        }
    //printf("Index = ");
    //print_vector(index, n);
    //printf(", N = %d\n", N);

    /* Store the values that each modulus can take */
    store_moduli = (int *)malloc(sizeof(int) * (N+q)*n);
    for (i = 0; i < N+q; i++)
        for (j = 0; j < n; j++)
            store_moduli[i*n+j] = i%index[j];

    ret_value = check_distances(n, Q, q, index, store_moduli);
    if (ret_value == FALSE)
        printf("Contradiction found\n");
    else
        printf("Everything a-OK\n");

    free(store_moduli);
    free(index);
    return ret_value;
}

/* Checks whether the q-ary to permutation mapping is valid
 * Input:
 *  n - the length of the q-ary vector
 *  Q - the number of permutation symbols used by each q-ary symbol
 *  q - the alphabet size
 *  index - the moduli
 *  store_moduli - the output of the moduli operations
 *
 * Output:
 *  TRUE or FALSE
 */
int check_distances(int n, int Q, int q, int *index, int *store_moduli)
{
    int *pvector, *qvector, *Pvector, *Qvector, *zero_vector;
    int N = Q*(n-1) + q;
    register unsigned long long i, j;
    unsigned long long num_vectors;

    // pvector = q-ary vector of length n
    // qvector = q-ary vector of length n
    // Pvector = permutation vector of length n corresponding to pvector
    // Qvector = permutation vector of length n corresponding to qvector
    // zero_vector = zero vector of length n
    num_vectors = (unsigned long long)pow(q, n);
    pvector = (int *)malloc(sizeof(int) * n);
    qvector = (int *)malloc(sizeof(int) * n);
    Pvector = (int *)malloc(sizeof(int) * N);
    Qvector = (int *)malloc(sizeof(int) * N);
    zero_vector = (int *)malloc(sizeof(int) * n);

    // Initial allocation
    for (i = 0; i < n; i++)
        zero_vector[i] = 0;
    set_qary_vector(pvector, zero_vector, n);
    set_qary_vector(qvector, zero_vector, n);
    free(zero_vector);

    /* The two for loops below loops through all pairs of vectors in the
     * q-ary space {0,1,...,q-1}^n. The first vector in the pair is pvector
     * and the second vector in the pair is qvector. qvector starts from
     * the vector following pvector, and continues on till it reaches the
     * last possible vector.
     */
    for (i = 0; i < num_vectors-1; i++)
    {
        for (j = i+1; j < num_vectors; j++)
        {
            generate_next_qary_vector(qvector, n, q);
            qary_to_permutations(pvector, Pvector, n, N, q, Q, index,
                                 store_moduli);
            qary_to_permutations(qvector, Qvector, n, N, q, Q, index,
                                 store_moduli);
            if (hamming_distance(Pvector, Qvector, N)
                < hamming_distance(pvector, qvector, n))
            {
                print_vector(pvector, n);
                printf(" -> ");
                print_vector(Pvector, N);
                printf("\n");
                print_vector(qvector, n);
                printf(" -> ");
                print_vector(Qvector, N);
                printf("\n");
                free(pvector);
                free(qvector);
                free(Pvector);
                free(Qvector);
                return FALSE;
            }
        }
        generate_next_qary_vector(pvector, n, q);
        set_qary_vector(qvector, pvector, n);
        print_vector(pvector, n);
        printf("\r");
    }
    printf("\n");
    free(pvector);
    free(qvector);
    free(Pvector);
    free(Qvector);
    return TRUE;
}

/* Given a vector of length n over an alphabet of size q, change the vector
 * to the next vector. If the vector is (v1, v2, ..., vn), then this
 * function increases it to (v1+1, v2, ..., vn) if possible. If v1 is q-1,
 * then it resets v1 to 0, and continues with v2.
 * Input:
 *  vector - the n length vector
 *  n - the length of the vector
 *  q - the alphabet size
 *
 * Output:
 *  None
 */
static inline void generate_next_qary_vector(int *vector, int n, int q)
{
    int i;
    for (i = 0; i < n; i++)
    {
        if (vector[i] == q-1)
            vector[i] = 0;
        else
        {
            vector[i]++;
            break;
        }
    }
}

/* Given two vectors of length n, determine the hamming distance between
 * them.
 * Input:
 *  x - the n length vector
 *  y - the n length vector
 *  n - the length of the vectors
 *
 * Output:
 *  int denoting the hamming distance.
 */
static inline int hamming_distance(int *x, int *y, int n)
{
    int d = 0, i;
    for (i=0; i<n; i++,x++,y++)
        d += (*x != *y);
    return d;
}

/* Given a vector of length n, print it nicely
 * Input:
 *  vector - the n length vector
 *  n - the length of the vector
 *
 * Output:
 *  None
 */
static inline void print_vector(int *vector, int n)
{
    int i;
    printf("( ");
    for (i = 0; i < n; i++)
        printf("%d ", vector[i]);
    printf(")");
}

/* Given q-ary vector, determine the permutation vector
 * Input:
 *  vector - the q-ary vector of length n
 *  Pvector - the permutation vector of length N
 *  n - the length of the q-ary vector
 *  N - the length of the permutation vector
 *  Q - the number of permutation symbols used by each q-ary symbol
 *  index - the moduli
 *  store_moduli - the output of the moduli operations
 *
 * Output:
 *  None
 */
void qary_to_permutations(int *vector, int *Pvector, int n, int N, int q,
                          int Q, int *index, int *store_moduli)
{
    int i, j, indx = 0;
    int *vec;
    for (i = 0; i < N; i++, Pvector++)
    {
        *Pvector = i;             // Initialization of Pvector
        if (i == index[indx])
            indx++;
        for (j = indx, vec = vector + indx; j < n; j++, vec++)
            *Pvector = *(store_moduli + (n*(*Pvector + *vec) + j));
    }
}

/* Copy the entries of values to the entries of vector
 * Input:
 *  vector - the n length vector whose entries are to be modified
 *  values - the n length vector from where the values are to be copied
 *  n - the length of the vector
 *
 * Output:
 *  None
 */
static inline void set_qary_vector(int *vector, int *values, int n)
{
    int i;
    for (i = 0; i < n; i++)
        vector[i] = values[i];
}

