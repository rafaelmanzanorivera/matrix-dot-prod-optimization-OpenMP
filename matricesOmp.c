//
//  matricesOmp.c
//  openmp
//
//  Created by Rafael Manzano on 24/9/18.
//  Copyright © 2018 Rafael Manzano. All rights reserved.
//


/*#define PARALLEL_CONF_OMP2M (_schedule,_parts,_i,_j,_k,_temporal) \
 {\
 _Pragma("omp parallel for schedule(static,parts)  private(i,j,k,temporal) shared(m1,m2,m3)")\
 }\*/


#include "matricesOmp.h"
#define MAT_SIZE 1500

#define PARALLEL_CONF_OMP3 __pragma("omp parallel for collapse(2) schedule(static,parts)  private(i,j,k,temporal) shared(m1,m2,m3) num_threads(nThreads)")

#define PARALLEL_CONF_OMP2_OUTER_FOR(type,parts,i,j,k,temporal,m1,m2,m3,nThreads) __pragma("omp parallel for schedule(#type,#parts)  private(#i,#j,#k,#temporal) shared(#m1,#m2,#m3) num_threads(#nThreads)")
#define PARALLEL_CONF_OMP2_INNER_FOR //__pragma("omp parallel for schedule(static,parts)  private(j,k,temporal) shared(m1,m2,m3) num_threads(nThreads)")
#define STATIC static
#define GUIDED guided
#define DYNAMIC dynamic



Matriz_t* crearMatriz(int numFilas, int numColumnas)
{
    int i = 0;
    
    Matriz_t* nMatriz = NULL;
    nMatriz = (Matriz_t*)malloc(sizeof(Matriz_t));
    
    nMatriz->numFilas = numFilas;
    nMatriz->numColumnas = numColumnas;
    
    nMatriz->datos = (float**)malloc(sizeof(float*)*numFilas);
    
    for (i = 0; i<numFilas; i++) {
        nMatriz->datos[i] = (float*)malloc(sizeof(float)*numColumnas);
    }
    //memset((void*)(nMatriz->datos[0][0]),0, sizeof(numFilas*numColumnas*sizeof(int)));
    return nMatriz;
}


void iniciaMatriz(Matriz_t *m)
{
    int i = 0;
    int j = 0;
    
    float t = 0.0;
    for (i = 0; i<m->numFilas; i++)
    {
        for (j = 0; j<m->numColumnas; j++)
        {
            m->datos[i][j] = (j + 1) / 6.0;
            t += 1.0;
        }
    }
}

int sonIguales(Matriz_t *m, Matriz_t*m2)
{
    int i = 0;
    int j = 0;
    
    for (i = 0; i<m->numFilas; i++)
    {
        for (j = 0; j<m->numColumnas; j++)
        {
            if (m->datos[i][j] != m2->datos[i][j]) return 0;
        }
    }
    return 1;
}

void imprimirMatriz(Matriz_t *m)
{
    int i, j;
    
    for (i = 0; i < m->numFilas; i++)
    {
        for (j = 0; j < m->numColumnas; j++)
        {
            printf("%f ", m->datos[i][j]);
        }
        printf("\n");
    }
    
    printf("\n");
}


double multiplicarMatricesSerie(Matriz_t *m1, Matriz_t *m2, Matriz_t *m3)
{
    int i, j, k;
    float temporal = 0.0;
    
    double start_time, run_time;
    
    start_time = omp_get_wtime();
    
    for (i = 0; i < m1->numFilas; i++) //i para las filas de la matriz resultante
    {
        for (j = 0; j < m2->numColumnas; j++) // k para las columnas de la matriz resultante
        {
            temporal = 0;
            for (k = 0; k < m1->numColumnas; k++) //j para realizar la multiplicacion de
            {                                   //los elementos   de la matriz
                temporal += m1->datos[i][j] * m2->datos[k][j];
            }
            m3->datos[i][j] = temporal;
        }
    }
    
    run_time = omp_get_wtime() - start_time;
    
    return run_time;
    
    
}

double multiplicarMatricesParaleloStatic(Matriz_t *m1, Matriz_t *m2, Matriz_t *m3, int nThreads, int parts)
{
    int i, j, k;
    float temporal = 0.0;
    
    double start_time, run_time;
    
    start_time = omp_get_wtime();
    
    
    
    
#ifdef _WIN32
    //PARALLEL_CONF_OMP2_OUTER_FOR(STATIC, parts, i, j, k, temporal, m1, m2, m3, nThreads);
#pragma omp parallel for schedule(static, parts)  private(j, k, temporal) shared(m1, m2, m3) num_threads(nThreads)
#else
    PARALLEL_CONF_OMP3;
#endif
    
    for (i = 0; i < m1->numFilas; i++)
    {
#ifdef _WIN32
        PARALLEL_CONF_OMP2_INNER_FOR;
#endif
        for (j = 0; j < m2->numColumnas; j++)
        {
            temporal = 0.0;
            
            
            for (k = 0; k < m1->numColumnas; k++)
            {
                temporal += m1->datos[i][j] * m2->datos[k][j];
            }
            
            m3->datos[i][j] = temporal;
        }
    }
    
    run_time = omp_get_wtime() - start_time;
    
    return run_time;
}

double multiplicarMatricesParaleloDynamic(Matriz_t *m1, Matriz_t *m2, Matriz_t *m3, int nThreads, int parts)
{
    int i, j, k;
    float temporal = 0.0;
    
    double start_time, run_time;
    
    start_time = omp_get_wtime();
    
    
    #ifdef _WIN32
        //PARALLEL_CONF_OMP2_OUTER_FOR(DYNAMIC, parts, i, j, k, temporal, m1, m2, m3, nThreads);
    #pragma omp parallel for schedule(dynamic, parts)  private(j, k, temporal) shared(m1, m2, m3) num_threads(nThreads)
    #else
        PARALLEL_CONF_OMP3;
    #endif
    
    for (i = 0; i < m1->numFilas; i++) //i para las filas de la matriz resultante
    {
    #ifdef _WIN32
            PARALLEL_CONF_OMP2_INNER_FOR;
    #endif
        for (j = 0; j < m2->numColumnas; j++) // k para las columnas de la matriz resultante
        {
            temporal = 0.0;
            
            
            for (k = 0; k < m1->numColumnas; k++)
            {
                temporal += m1->datos[i][j] * m2->datos[k][j];
            }
            
            m3->datos[i][j] = temporal;
        }
    }
    
    run_time = omp_get_wtime() - start_time;
    
    return run_time;
}

double multiplicarMatricesParaleloGuided(Matriz_t *m1, Matriz_t *m2, Matriz_t *m3, int nThreads, int parts)
{
    int i, j, k;
    float temporal = 0.0;
    double start_time, run_time;
    
    start_time = omp_get_wtime();
    
    #ifdef _WIN32
        //PARALLEL_CONF_OMP2_OUTER_FOR(GUIDED, parts, i, j, k, temporal, m1, m2, m3, nThreads);
    #pragma omp parallel for schedule(guided, parts)  private(j, k, temporal) shared(m1, m2, m3) num_threads(nThreads)
    #else
        PARALLEL_CONF_OMP3;
    #endif
    
    
    for (i = 0; i < m1->numFilas; i++)
    {
    #ifdef _WIN32
            PARALLEL_CONF_OMP2_INNER_FOR;
    #endif
        for (j = 0; j < m2->numColumnas; j++)
        {
            temporal = 0.0;
            
            for (k = 0; k < m1->numColumnas; k++)
            {
                temporal += m1->datos[i][j] * m2->datos[k][j];
            }
            m3->datos[i][j] = temporal;
        }
    }
    
    run_time = omp_get_wtime() - start_time;
    
    return run_time;
}

void multiplicarYComparar(Matriz_t *m1, Matriz_t *m2)
{
    Matriz_t *mresSerial, *mresParallel;
    mresSerial = crearMatriz(MAT_SIZE, MAT_SIZE);
    mresParallel = crearMatriz(MAT_SIZE, MAT_SIZE);
    
    int i, j;
    int n_threads_max = 9;
    int n_threads_step = 2;
    int n_parts_max = 200;
    int n_parts_step = 60;
    
    double tSerial, tParallel;
    
    printf("Multiplicando en Serie...\n");
    
    tSerial = multiplicarMatricesSerie(m1, m2, mresSerial);
    printf("Tiempo en serie %lf\n", tSerial);
    
    //TODO mediana cinco ejecuciones
    
    printf("\n\nMultiplicando en Paralelo Schedule Static\n");
    for (int i = 2; i<n_threads_max; i += n_threads_step)
    {
        printf("\n%d threads -- \n\n", i);
        for (j = 1; j < n_parts_max; j += n_parts_step)
        {
            tParallel = multiplicarMatricesParaleloStatic(m1, m2, mresParallel, i, j);
            printf("LS= %d -> tm= %lf, ", j, tParallel);
            printf("sp= %lf ", tSerial / tParallel);
            printf("(%d)\n", sonIguales(mresSerial, mresParallel));
        }
        
    }
    
    printf("\n\nMultiplicando en Paralelo Schedule Dynamic\n");
    for (int i = 2; i<n_threads_max; i += n_threads_step)
    {
        printf("\n%d threads -- \n\n", i);
        for (j = 1; j < n_parts_max; j += n_parts_step)
        {
            tParallel = multiplicarMatricesParaleloDynamic(m1, m2, mresParallel, i, j);
            printf("LS= %d -> tm= %lf, ", j, tParallel);
            printf("sp= %lf ", tSerial / tParallel);
            printf("(%d)\n", sonIguales(mresSerial, mresParallel));
        }
        
    }
    
    printf("\n\nMultiplicando en Paralelo Schedule Guided\n");
    for (int i = 2; i<n_threads_max; i += n_threads_step)
    {
        printf("\n%d threads -- \n\n", i);
        for (j = 1; j < n_parts_max; j += n_parts_step)
        {
            tParallel = multiplicarMatricesParaleloGuided(m1, m2, mresParallel, i, j);
            printf("LS= %d -> tm= %lf, ", j, tParallel);
            printf("sp= %lf ", tSerial / tParallel);
            printf("(%d)\n", sonIguales(mresSerial, mresParallel));
        }
        
    }
    
    
    
    
    
    
    
}

int main(int argc, char** argv)
{
    Matriz_t *m1, *m2;
    //
    //    m1 = crearMatriz(10, 10);
    //    m2 = crearMatriz(10, 10);
    //    mres = crearMatriz(10, 10);
    //
    printf("Creando Matrices...\n");
    
    m1 = crearMatriz(MAT_SIZE, MAT_SIZE);
    m2 = crearMatriz(MAT_SIZE, MAT_SIZE);
    
    
    printf("Inicializando Matrices...\n");
    
    iniciaMatriz(m1);
    iniciaMatriz(m2);
    
    multiplicarYComparar(m1, m2);
    
    
    getchar();
    
    
    
    
    
    
    
}
