#pragma once
//
// Created by c50008023 on 2019/8/5.
//

#include"parameters.h"

using namespace std;

class Solver
{
public:

	Solver();							//init the array that needs to be calculated from 0-20
	~Solver();
	void Initial(int task, int island);						//init the beginning population
	void InitialForOnce();
	void InitForPQ(int task);
	void UpdateByTask(int task, int island);		 //pick the highest probability as father and mother from population
	void SBX(int task, float* p1, float* p2, float lower, float upper, int dimension);	 //crossover
	void Mutation(float* offspring, float lower, float upper);
	int getRunId();

	void SingleUpdateByTask(int task);
	void evaluateByTask(int task_id, int island);
	void generateSolution(int task, int island);
	void exGenerateSolution(int task, int island);
	void wholeProcess(string file_path, int run_id);
	void migration(int task, int island);
	void explicitTransfer(int target,int island);
	void implicitTransfer(int target,int island);
	static bool sortFunc(Chromosome a, Chromosome b);

	float ChangeToTaskScale(int task, const float gen_value);
	double getGenerateTime();
	double getTransferTime();
	double getMigrationTime();
	double getSortingTime();

	Task tasks[TASK_SIZE];
	//string input_file = "./inputs/";
//	string input_file = "D:\\Visual Stdio\\Island Model\\inputs\\";
	string input_file;


    //thread thread_array_transfer[TASK_SIZE * ISLAND_NUM * TRANSFER_SIZE];
	//thread thread_array_generation[TASK_SIZE * ISLAND_NUM * ISLAND_SIZE / 2];
	//thread thread_array_migration[TASK_SIZE * ISLAND_NUM * MIGRATION_SIZE];

	//task for thread
	//thread thread_array_task[TASK_SIZE];
	//thread thread_array_generation[TASK_SIZE];
	//thread thread_array_transfer[TASK_SIZE];
	//thread thread_array_migration[TASK_SIZE];

	thread init_array[TASK_SIZE * ISLAND_NUM];
	thread thread_array_evaluate[TASK_SIZE * ISLAND_NUM];
	//island for thread
	thread thread_array_island[TASK_SIZE * ISLAND_NUM];
	thread thread_array_generation[TASK_SIZE * ISLAND_NUM];
	thread thread_array_transfer[TASK_SIZE * ISLAND_NUM];
	thread thread_array_migration[TASK_SIZE * ISLAND_NUM];
    thread thread_array_task[TASK_SIZE];

private:
	int run_id;

	//time
	double start_time;
	double migration_time;
	double transfer_time;
	double evolve_time;
	double sorting_time;
	double evaluate_time;
    double generate_time;

	float getRandomFloat(float minRange, float maxRange);
	int getRandomInt(int minRange, int maxRange);
	void CheckDomain(float& res1, float lower, float upper);
	static bool compareFunc(Chromosome* chrom1, Chromosome* chrom2);
	void transferByM(float* elem, float* M, float* res, int row, int col);
	int rouletteWheel();
	void calculateM(float* M, int target, int source);
	void transform(float* res, float* M, float* elem, float bias, int dimension, int task);

};
