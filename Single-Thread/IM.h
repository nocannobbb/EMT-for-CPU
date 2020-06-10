#pragma once
//
// Created by c50008023 on 2019/8/5.
//

#include"parameters.h"

using namespace std;


class Solver
{
public:

	Solver(int run_id);							//init the array that needs to be calculated from 0-20
	~Solver();
	void Initial(float& sorting_time, float& evaluate_time);						//init the beginning population
	void InitialForPQ();
	float Fitness(Chromosome* chrom);			//calculate the fitness, and return int
	void UpdateByTask(int task);		 //pick the highest probability as father and mother from population
	void SBX(int task, float* p1, float* p2, float lower, float upper, int dimension);	 //crossover
	void Mutation(float* offspring, float lower, float upper);
	int getRunId();
	float ChangeToTaskScale(int task, const float gen_value);

    void evaluateByTask(int task_id);
	void generateSolution(string file_path);
	void migration();
	void ExplicitTransfer();
	void ImplicitTransfer();
	static bool sortFunc(Chromosome a, Chromosome b);


	Task tasks[TASK_SIZE];
	string input_file = "D:\\Visual Stdio\\Island Model\\inputs\\";   //own pc
//	string input_file = "D:\\EMT\\DemoCode\\MFEA-ISLAND\\inputs\\";     //huawei
    //string input_file = "C:\\Users\\user5\\Desktop\\cm\\Island_Model-master\\MFEA-ISLAND\\inputs\\";    //



private:
    int run_id;

	float ToDec(float* bit);
	float getRandomFloat(float minRange, float maxRange);
	int getRandomInt(int minRange, int maxRange);
	void CheckDomain(float &res1, float lower, float upper);
	static bool compareFunc(Chromosome* chrom1, Chromosome* chrom2);
	void transferByM(float* elem, float* M, float* res, int row, int col);
	int rouletteWheel();
	void calculateM(float* M, int target, int source);
	void transform(float *res, float* M, float* elem, float bias, int dimension, int task);
};
