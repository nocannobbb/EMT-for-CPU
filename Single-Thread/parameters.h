#pragma once
//
// Created by c50008023 on 2019/8/5.
//
#include<iostream>
#include <time.h>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <cassert>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <assert.h>
#include <armadillo>
#include <direct.h>



using namespace arma;

#define	ETA 5
#define MU	2
#define RAND_LOWER_01 0
#define RAND_UPPER_01 1.0
#define LOWER	0
#define UPPER	1.0

//added in later version
#define CLEAN 1
//#define PI 3.1415926535
const float PI = acos(-1);



//the decision of explicit transfer and the implicit transfer
const int IMPLICIT_TRANSFER = 0;
const float GEN_LOWER = 0.0;
const float GEN_UPPER = 1.0;
const int RUN_SIZE = 4;


//get from gpu version
// Parameters for MFEA-2
const int TASK_SIZE = 1000;
const int POPULATION_SIZE = 512;
const int MAX_DIMENSION = 50;
const int GENERATION_SIZE = 1000;
const int TRANSFER_INTERVAL = 50;


/* Parameters for GPU */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
const int ISLAND_SIZE = 128;
const int TRANSFER_SIZE = 20;
const int ISLAND_NUM = POPULATION_SIZE / ISLAND_SIZE; // number of islands
const int MIGRATION_SIZE = 20;
const int MIGRATION_INTERVAL = 20;


/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


/* Parameters for optimization functions */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
const float BIAS = 1e-3;
const int FUNC_SIZE = 7;
const float FUNC_DOMAIN[FUNC_SIZE] = { 100, 50, 50, 50, 100, 0.5, 500 };
const int BIAS_SIZE = 5;
const float BIAS_DOMAIN[BIAS_SIZE] = { 0, 0.1, -0.1, 0.2, -0.2 };
/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */



struct Chromosome {
	float phenotype;
	float genotype[MAX_DIMENSION];

	//bool operator< (const Chromosome& rhs) const {
	//	return phenotype < rhs.phenotype;
	//}
};

struct Auxiliary {
	float lower;
	float upper;
	float bias;
	int dimension;
	float M[MAX_DIMENSION * MAX_DIMENSION];
};

struct Task {
	mat Q;
	mat P;
	Auxiliary auxiliary;
	Chromosome* chromosome[POPULATION_SIZE];
};


//evaluate
extern float func0(float* elem, int& dimension);
extern float func1(float* elem, int& dimension);
extern float func2(float* elem, int& dimension);
extern float func3(float* elem, int& dimension);
extern float func4(float* elem, int& dimension);
extern float func5(float* elem, int& dimension);
extern float func6(float* elem, int& dimension);
extern float functionChoice(int& task_id, float* elem, int& dimension);
//extern void transform(float* M, float* elem, float bias, int dimension);
//void evaluateTasks(Task* tasks, int offset);
extern void cpy(float* target, float* source, int dimension);




//fileOperation
extern std::string nts(int number);
extern void transMatrixRead(std::string input_file, float* matrix, int dimension);
