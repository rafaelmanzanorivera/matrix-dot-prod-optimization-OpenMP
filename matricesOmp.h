//
//  matricesOmp.h
//  openmp
//
//  Created by Rafael Manzano on 24/9/18.
//  Copyright © 2018 Rafael Manzano. All rights reserved.
//

#ifndef matricesOmp_h
#define matricesOmp_h

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

typedef struct Matriz_t
{
    int numFilas;
    int numColumnas;
    float **datos;
    
}Matriz_t;

Matriz_t* crearMatriz(int numFilas, int numColumnas);

void iniciaMatriz(Matriz_t *m);

void multiplicarMatrices(Matriz_t *m1, Matriz_t *m2, Matriz_t *m3);

void multiplicarYComparar(Matriz_t *m1, Matriz_t *m2);

void imprimirMatriz(Matriz_t *m);


#endif /* matricesOmp_h */
