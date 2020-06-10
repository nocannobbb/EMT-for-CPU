//
// Created by c50008023 on 2019/8/5.
//

#include "IM.h"
#include "parameters.h"
#define OFFSET_RANDOM 0.999999;


using namespace arma;

Solver::Solver()
{
	//no use
	srand(time(NULL));

	for (int task = 0; task < TASK_SIZE; ++task)
	{
		for (int popu = 0; popu < POPULATION_SIZE; ++popu)
		{
			tasks[task].chromosome[popu] = new Chromosome();
		}
	}



}

Solver::~Solver()
{
	for (int task = 0; task < TASK_SIZE; ++task)
	{
		for (int popu = 0; popu < POPULATION_SIZE; ++popu)
		{
			delete tasks[task].chromosome[popu];
			tasks[task].chromosome[popu] = NULL;
		}
	}
}

void Solver::InitForPQ(int task)
{
    tasks[task].Q = zeros(MAX_DIMENSION + 1, POPULATION_SIZE);
    tasks[task].P = zeros(MAX_DIMENSION, POPULATION_SIZE);

    sort(tasks[task].chromosome,tasks[task].chromosome + POPULATION_SIZE, compareFunc);
    for(int gen = 0; gen < MAX_DIMENSION; ++gen)
    {
        for(int popu = 0; popu < POPULATION_SIZE; ++popu)
        {
            tasks[task].P(gen,popu) = tasks[task].chromosome[popu]->genotype[gen];
            tasks[task].Q(gen,popu) = tasks[task].chromosome[popu]->genotype[gen];
        }
    }

}

void Solver::evaluateByTask(int task_id, int island) {

	for (int popu = 0; popu < ISLAND_SIZE; ++popu)
	{
		int i = island * ISLAND_SIZE + popu;
		float gene[MAX_DIMENSION];
		//        cpy(gene, tasks[task_id].chromosome[i]->genotype, tasks[task_id].auxiliary.dimension);
		transform(gene, tasks[task_id].auxiliary.M, tasks[task_id].chromosome[i]->genotype,
			tasks[task_id].auxiliary.bias, tasks[task_id].auxiliary.dimension, task_id);

		tasks[task_id].chromosome[i]->phenotype = functionChoice(task_id,
			gene, tasks[task_id].auxiliary.dimension);

	}
}


void Solver::Initial(int task, int island)
{
	//srand(random_num * task * task * island * island);

	LARGE_INTEGER nFrequency;

	if (::QueryPerformanceFrequency(&nFrequency))
	{
		LARGE_INTEGER nStartCounter;

		::QueryPerformanceCounter(&nStartCounter);

		::srand((unsigned)nStartCounter.LowPart);

	}
	else // 否则使用当前系统时间初始化随机种子(毫秒级)
	{
		struct timeb stb;

		::ftime(&stb);

		::srand((unsigned)stb.millitm);

	}

	for (int popu = 0; popu < ISLAND_SIZE; ++popu)
	{
		for (int j = 0; j < MAX_DIMENSION; ++j)
		{
		    if(IMPLICIT_TRANSFER)
            {
                tasks[task].chromosome[island*ISLAND_SIZE + popu]->genotype[j] = this->getRandomFloat(GEN_LOWER, GEN_UPPER);
            }
            else{
                tasks[task].chromosome[island*ISLAND_SIZE + popu]->genotype[j] = this->getRandomFloat(tasks[task].auxiliary.lower, tasks[task].auxiliary.upper);
            }

		}
	}

	//sort(tasks[task].chromosome, tasks[task].chromosome + POPULATION_SIZE, compareFunc);
}

void Solver::InitialForOnce()
{
	for (int task = 0; task < TASK_SIZE; ++task)
	{
		//initialize auxiliary
		tasks[task].auxiliary.bias = BIAS_DOMAIN[task % BIAS_SIZE] * FUNC_DOMAIN[task % FUNC_SIZE];
		tasks[task].auxiliary.lower = -FUNC_DOMAIN[task % FUNC_SIZE];
		tasks[task].auxiliary.upper = FUNC_DOMAIN[task % FUNC_SIZE];
		tasks[task].auxiliary.dimension = MAX_DIMENSION;
		//M
		transMatrixRead(input_file + "matrix/" + nts(task + 1) + ".txt", tasks[task].auxiliary.M, tasks[task].auxiliary.dimension);
	}



	//initialize task phenotype output file
//	for (int task = 0; task < TASK_SIZE; ++task)
//	{
//		string output_file = "D:\\Visual Stdio\\Island Model\\IslandModel\\IslandModel\\output\\task " + to_string(task) + ".txt";
//		//string output_file = "C:\\Users\\user5\\Desktop\\cm\\Island_Model\\Island_Model\\output\\problem_" + to_string(task) + "_" +to_string(run_id) + ".txt";
//		ofstream outfile;
//		outfile.open(output_file, ios::out | ios::trunc);
//		outfile << "worst individual : " << tasks[task].chromosome[POPULATION_SIZE - 1]->phenotype << endl;
//		outfile.close();
//	}

	//intialize time output file
//	string output_file = "D:\\Visual Stdio\\Island Model\\IslandModel\\IslandModel\\output\\time\\" + to_string(run_id) + ".txt";
//	//string output_file = "C:\\Users\\user5\\Desktop\\cm\\Island_Model\\Island_Model\\output\\time\\" + to_string(run_id) + ".txt";
//	ofstream outfile;
//	outfile.open(output_file, ios::out | ios::trunc);
//	outfile.close();
}

float Solver::ChangeToTaskScale(int task, const float gen_value)
{
	float task_lower = tasks[task].auxiliary.lower;
	float task_upper = tasks[task].auxiliary.upper;
	float gen_lower = GEN_LOWER;
	float gen_upper = GEN_UPPER;

	float regular_value = (gen_value - gen_lower) / (gen_upper - gen_lower);
	float eva_value = regular_value * (task_upper - task_lower) + task_lower;
	return eva_value;
}



bool Solver::compareFunc(Chromosome* chrom1, Chromosome* chrom2)
{
	return chrom1->phenotype < chrom2->phenotype;
}

bool Solver::sortFunc(Chromosome a, Chromosome b)
{
	return a.phenotype < b.phenotype;
}

void Solver::UpdateByTask(int task,int island)
{
	//sort by every island
//	int i = 0;
//	while(i > 1);
	sort(tasks[task].chromosome + island * ISLAND_SIZE, tasks[task].chromosome + island * ISLAND_SIZE + ISLAND_SIZE, compareFunc);
	//cout << "island " << island + 1 << " : " << tasks[task].chromosome[island*ISLAND_SIZE]->phenotype << endl;
}

int Solver::rouletteWheel() {
	int rand1 = getRandomInt(0, ISLAND_SIZE);
	int rand2 = getRandomInt(0, ISLAND_SIZE);
	return min(rand1, rand2);
}


void Solver::Mutation(float* offspring, float lower, float upper)
{
	for (int i = 0; i < MAX_DIMENSION; ++i)
	{
		float randf = getRandomFloat(RAND_LOWER_01, RAND_UPPER_01);
		if (randf < 1.0 / MAX_DIMENSION)
		{
			float u = getRandomFloat(RAND_LOWER_01, RAND_UPPER_01);
			float temp = offspring[i];
			if (u <= 0.5)
			{
				float del = powf(2.0 * u, 1.0 / (1.0 + ETA)) - 1.0;
				float res = temp + del * temp;
				CheckDomain(res, lower, upper);
				offspring[i] = res;
			}
			else
			{
				float del = 1 - powf(2.0 * (1 - u), 1.0 / (1.0 + ETA));
				float res = temp + del * (1.0 - temp);
				CheckDomain(res, lower, upper);
				offspring[i] = res;
			}

		}
	}

	//cout << "debug mutation :" << Fitness(&offspring) << endl;
}

//simulated binary crossover,using p1 and p2 to generate off1,off2
void Solver::SBX(int task, float *p1, float *p2, float lower, float upper, int dimension)
{
	float cf;
	for (int i = 0; i < dimension; ++i)
	{
		float u = getRandomFloat(RAND_LOWER_01, RAND_UPPER_01);
		if (u <= 0.5) {
			cf = powf(2.0 * u, 1.0 / (MU + 1));
		}
		else {
			cf = powf(1.0 / (2.0 * (1 - u)), 1.0 / (MU + 1));
		}
//		float res1 = 0.5 * ((1 + cf) * (tasks[task].chromosome[p1]->genotype[i]) + (1 - cf) * (tasks[task].chromosome[p2]->genotype[i]));
//		float res2 = 0.5 * ((1 + cf) * (tasks[task].chromosome[p2]->genotype[i]) + (1 - cf) * (tasks[task].chromosome[p1]->genotype[i]));

        float res1 = 0.5 * ((1 + cf) * (p1[i]) + (1 - cf) * (p2[i]));
		float res2 = 0.5 * ((1 + cf) * (p2[i]) + (1 - cf) * (p1[i]));

		CheckDomain(res1, lower, upper);	//check border, it's necessary
		CheckDomain(res2, lower, upper);
		p1[i] = res1;
		p2[i] = res2;

	}
}

int Solver::getRunId()
{
	return this->run_id;
}


void Solver::CheckDomain(float& res, float lower, float upper)
{
	if (res < lower)	res = lower;
	if (res > upper)	res = upper;
}

//get the random float value from minRange to maxRange
//float Solver::getRandomFloat(float minRange, float maxRange)
//{
//
//	int right = floor(maxRange);
//	int left = ceil(minRange);
//	int intPart;
//	if(right == left)
//    {
//        intPart = 0;
//    }
//    else{
//        intPart = rand() % (right - left) + left;
//    }
//
//	float floatPart = rand() % 10;
//	floatPart /= 10.0;
//	float p = (float)intPart + floatPart;
//
//	return p;
//}


float Solver::getRandomFloat(float minRange, float maxRange)
{
	float fResult;

	if (minRange > maxRange)
	{
		int temp = minRange;
		minRange = maxRange;
		maxRange = temp;
	}

	fResult = minRange + (maxRange - minRange) * rand() / (RAND_MAX + 1);

	return fResult;
}

int Solver::getRandomInt(int minRange, int maxRange)
{
	return rand() % (maxRange - minRange) + minRange;
}


void Solver::transform(float* res, float* M, float* elem, float bias, int dimension, int task) {
	//	float res[MAX_DIMENSION];
	float temp[MAX_DIMENSION];

	cpy(temp, elem, MAX_DIMENSION);

	if(IMPLICIT_TRANSFER)
    {
        for (int i = 0; i < dimension; ++i)
        {
            temp[i] = ChangeToTaskScale(task, elem[i]);
        }
    }

	for (int i = 0; i < dimension; ++i) {
		float sum = 0;
		for (int j = 0; j < dimension; ++j) {
			sum += M[i * dimension + j] * (temp[j] + bias);
		}
		res[i] = sum;
	}
	//	cpy(elem, res, dimension);


}

void Solver::exGenerateSolution(int task, int island)
{
    Chromosome offspring1;
	Chromosome offspring2;
	int p1, p2, dimension;
	float lower, upper, bias;
	float M[MAX_DIMENSION * MAX_DIMENSION] = { 0.0 };


	cpy(M, tasks[task].auxiliary.M, MAX_DIMENSION * MAX_DIMENSION);
	bias = tasks[task].auxiliary.bias;
	lower = tasks[task].auxiliary.lower;
	upper = tasks[task].auxiliary.upper;
	dimension = tasks[task].auxiliary.dimension;


	for (int popu = island * ISLAND_SIZE; popu < ISLAND_SIZE / 2 + island * ISLAND_SIZE; ++popu)
	{

		//float evolve_start = clock();
		p1 = this->rouletteWheel();
		p2 = this->rouletteWheel();
		if (p1 == p2) p2 = (p2 + 1) % ISLAND_SIZE;

		p1 += island * ISLAND_SIZE;
		p2 += island * ISLAND_SIZE;

		cpy(offspring1.genotype, tasks[task].chromosome[p1]->genotype, dimension);
		cpy(offspring2.genotype, tasks[task].chromosome[p2]->genotype, dimension);

		if(IMPLICIT_TRANSFER)
        {
            this->SBX(task, offspring1.genotype, offspring2.genotype, GEN_LOWER, GEN_UPPER, dimension);
            this->Mutation(offspring1.genotype, GEN_LOWER, GEN_UPPER);
            this->Mutation(offspring2.genotype, GEN_LOWER, GEN_UPPER);
        }
        else{
            this->SBX(task, offspring1.genotype, offspring2.genotype, lower, upper, dimension);
            this->Mutation(offspring1.genotype, lower, upper);
            this->Mutation(offspring2.genotype, lower, upper);
        }


//		float gene1[MAX_DIMENSION], gene2[MAX_DIMENSION];
//		//float evaluate_start = clock();
//		transform(gene1, M, offspring1.genotype, bias, dimension, task);
//		transform(gene2, M, offspring2.genotype, bias, dimension, task);
//
//		float res1 = functionChoice(task, gene1, dimension);
//		float res2 = functionChoice(task, gene2, dimension);
//		//evaluate_time += clock() - evaluate_start;
//
//		//evolve_start = clock();
//		//all assigned to p1, because we do for loop for p1, size: 1 / 2 * ISLAND_SIZE
//		if (tasks[task].chromosome[popu]->phenotype - BIAS > res1)
//		{
//			cpy(tasks[task].chromosome[popu]->genotype, offspring1.genotype, dimension);
//			tasks[task].chromosome[popu]->phenotype = res1;
//		}
//		if (tasks[task].chromosome[popu + ISLAND_SIZE / 2]->phenotype - BIAS > res2)
//		{
//			cpy(tasks[task].chromosome[popu + ISLAND_SIZE / 2]->genotype, offspring2.genotype, dimension);
//			tasks[task].chromosome[popu + ISLAND_SIZE / 2]->phenotype = res2;
//		}
		//evolve_time += clock() - evolve_start;
	}
}

void Solver::generateSolution(int task, int island)
{
	Chromosome offspring1;
	Chromosome offspring2;
	int p1, p2, dimension;
	float lower, upper, bias;
	float M[MAX_DIMENSION * MAX_DIMENSION] = { 0.0 };


	cpy(M, tasks[task].auxiliary.M, MAX_DIMENSION * MAX_DIMENSION);
	bias = tasks[task].auxiliary.bias;
	lower = tasks[task].auxiliary.lower;
	upper = tasks[task].auxiliary.upper;
	dimension = tasks[task].auxiliary.dimension;


	for (int popu = island * ISLAND_SIZE; popu < ISLAND_SIZE / 2 + island * ISLAND_SIZE; ++popu)
	{

		//float evolve_start = clock();
		p1 = this->rouletteWheel();
		p2 = this->rouletteWheel();
		if (p1 == p2) p2 = (p2 + 1) % ISLAND_SIZE;

		p1 += island * ISLAND_SIZE;
		p2 += island * ISLAND_SIZE;

		cpy(offspring1.genotype, tasks[task].chromosome[p1]->genotype, dimension);
		cpy(offspring2.genotype, tasks[task].chromosome[p2]->genotype, dimension);

		if(IMPLICIT_TRANSFER)
        {
            this->SBX(task, offspring1.genotype, offspring2.genotype, GEN_LOWER, GEN_UPPER, dimension);
            this->Mutation(offspring1.genotype, GEN_LOWER, GEN_UPPER);
            this->Mutation(offspring2.genotype, GEN_LOWER, GEN_UPPER);
        }
        else{
            this->SBX(task, offspring1.genotype, offspring2.genotype, lower, upper, dimension);
            this->Mutation(offspring1.genotype, lower, upper);
            this->Mutation(offspring2.genotype, lower, upper);
        }


		//evolve_time += (clock() - evolve_start);
		/*cpy(p1_temp, offspring1.genotype, dimension);
		cpy(p2_temp, offspring2.genotype, dimension);*/

		float gene1[MAX_DIMENSION], gene2[MAX_DIMENSION];
		//float evaluate_start = clock();
		transform(gene1, M, offspring1.genotype, bias, dimension, task);
		transform(gene2, M, offspring2.genotype, bias, dimension, task);

		float res1 = functionChoice(task, gene1, dimension);
		float res2 = functionChoice(task, gene2, dimension);
		//evaluate_time += clock() - evaluate_start;

		//evolve_start = clock();
		//all assigned to p1, because we do for loop for p1, size: 1 / 2 * ISLAND_SIZE
		if (tasks[task].chromosome[popu]->phenotype - BIAS > res1)
		{
			cpy(tasks[task].chromosome[popu]->genotype, offspring1.genotype, dimension);
			tasks[task].chromosome[popu]->phenotype = res1;
		}
		if (tasks[task].chromosome[popu + ISLAND_SIZE / 2]->phenotype - BIAS > res2)
		{
			cpy(tasks[task].chromosome[popu + ISLAND_SIZE / 2]->genotype, offspring2.genotype, dimension);
			tasks[task].chromosome[popu + ISLAND_SIZE / 2]->phenotype = res2;
		}
		//evolve_time += clock() - evolve_start;
	}

}

void Solver::wholeProcess(string file_path, int run_id)
{
    start_time = clock();
    migration_time = 0;
	transfer_time = 0;
	evolve_time = 0;
	sorting_time = 0;
    generate_time = 0;

    double migartion_start = 0;
    double transfer_start = 0;
    double generate_start = 0;
    double sorting_start = 0;
    double cur = 0;
	//initial
	this->InitialForOnce();

	//multi thread initial
	for (int task = 0; task < TASK_SIZE; ++task)
	{
		for (int island = 0; island < ISLAND_NUM; ++island)
		{
			init_array[task*ISLAND_NUM + island] = thread(&Solver::Initial, this, task, island);
		}
	}
	for (auto& t : init_array)
	{
		t.join();
	}


	for (int task = 0; task < TASK_SIZE; ++task)
	{
		for (int island = 0; island < ISLAND_NUM; ++island)
		{
			thread_array_evaluate[task * ISLAND_NUM + island] = thread(&Solver::evaluateByTask , this , task , island);
		}
	}
	for(auto&t : thread_array_evaluate)
    {
        t.join();
    }

    if(!IMPLICIT_TRANSFER)
    {
        cur = clock();
        //init for PQ
        thread_array_task[TASK_SIZE];
        for(int task = 0; task < TASK_SIZE; ++task)
        {
            thread_array_task[task] = thread(&Solver::InitForPQ, this , task);
        }
        for(auto& t: thread_array_task)
        {
            t.join();
        }
        transfer_time += clock() - cur;
    }


	//beginning
	for (int task = 0; task < min(7, TASK_SIZE); ++task)
	{
		cout << "island beginning  " << task + 1 << " : " << tasks[task].chromosome[0]->phenotype << endl;
		//for (int island = 0; island < ISLAND_NUM; ++island)
		//{
		//	cout << "task:    " << task + 1 << " :     " <<
		//		tasks[task].chromosome[island * ISLAND_SIZE]->phenotype << endl;
		//}
	}



	//must put it outside, can't put it into for loop, or the random values will not change
	srand(time(NULL));
	//for loop each value in the front half,and randomly choose a item in the whole array,then do SBX between the two

	for (int generation = 1; generation <= GENERATION_SIZE; ++generation)
	{
		 //record the one-minute results
//	    if(clock() - start_time >= 60000)
//        {
//            //output best individual into files
//            for(int task = 0; task < TASK_SIZE; ++task)
//            {
//                string output_file = file_path + "TASKS_RUN_" + to_string(run_id) + ".txt";
//                ofstream outfile;
//                outfile.open(output_file, ios::app);
//                //best individual of all populations
//                sort(tasks[task].chromosome,tasks[task].chromosome + POPULATION_SIZE,compareFunc);
//                outfile << "best individual : " << tasks[task].chromosome[0]->phenotype << endl;
//
//                outfile.close();
//
//            }
//            return;
//        }


		//migration
		if (generation % MIGRATION_INTERVAL == 0)
		{
		    migartion_start = clock();
			for (int task = 0; task < TASK_SIZE; ++task)
			{
				for (int island = 0; island < ISLAND_NUM; ++island)
				{
					//int popid = task * ISLAND_NUM * MIGRATION_SIZE + island * MIGRATION_SIZE + subpopid;
					thread_array_migration[task*ISLAND_NUM + island] = thread(&Solver::migration, this, task, island);
				}
			}
			for (auto& t : thread_array_migration)
			{
				t.join();
			}

			for (int task = 0; task < TASK_SIZE; ++task)
			{
				for (int island = 0; island < ISLAND_NUM; ++island)
				{
					//this->UpdateByTask(task);
					thread_array_island[task*ISLAND_NUM + island] = thread(&Solver::UpdateByTask, this, task, island);
				}
			}
			for (auto& t : thread_array_island)
			{
				t.join();
			}
            migration_time += (clock() - migartion_start);
		}


		//transfer
		if (generation % TRANSFER_INTERVAL == 0)
		{
			transfer_start = clock();
			for (int target = 0; target < TASK_SIZE; ++target)
			{
				for (int island = 0; island < ISLAND_NUM; ++island)
				{
					//int popid = target * ISLAND_NUM * TRANSFER_SIZE + island * TRANSFER_SIZE + popu;
					if(IMPLICIT_TRANSFER)
                    {
                        thread_array_transfer[target*ISLAND_NUM + island] = thread(&Solver::implicitTransfer, this, target, island);
                    }
                    else{
                        thread_array_transfer[target*ISLAND_NUM + island] = thread(&Solver::explicitTransfer, this, target, island);
                    }
				}
			}
			for (auto& t : thread_array_transfer)
			{
				t.join();
			}

			//evaluate in transfer
			for (int task = 0; task < TASK_SIZE; ++task)
			{
				for (int island = 0; island < ISLAND_NUM; ++island)
				{
					//this->evaluateByTask(task);
					thread_array_island[task*ISLAND_NUM + island] = thread(&Solver::evaluateByTask, this, task, island);
				}
			}
			for (auto& t : thread_array_island)
			{
				t.join();
			}

			//sort in transfer
			for (int task = 0; task < TASK_SIZE; ++task)
			{
				for (int island = 0; island < ISLAND_NUM; ++island)
				{
					thread_array_island[task*ISLAND_NUM + island] = thread(&Solver::UpdateByTask, this, task, island);
				}
			}
			for (auto& t : thread_array_island)
			{
				t.join();
			}

			transfer_time += (clock() - transfer_start);
		}


        //exgenerateSolution
//		exgenerate_start = clock();
//		for (int task = 0; task < TASK_SIZE; ++task)
//		{
//			for (int island = 0; island < ISLAND_NUM; ++island)
//			{
//				//int popid = task * ISLAND_NUM * ISLAND_SIZE / 2 + island * ISLAND_SIZE / 2 + popu;
//				thread_array_generation[task * ISLAND_NUM + island] = thread(&Solver::exGenerateSolution, this, task, island);
//			}
//		}
//		for (auto& t : thread_array_generation)
//		{
//			t.join();
//		}
//        exgenerate_time += (clock() - exgenerate_start);


		// SBX Mutation Part
		generate_start = clock();
		for (int task = 0; task < TASK_SIZE; ++task)
		{
			for (int island = 0; island < ISLAND_NUM; ++island)
			{
				//int popid = task * ISLAND_NUM * ISLAND_SIZE / 2 + island * ISLAND_SIZE / 2 + popu;
				thread_array_generation[task * ISLAND_NUM + island] = thread(&Solver::generateSolution, this, task, island);
			}
		}
		for (auto& t : thread_array_generation)
		{
			t.join();
		}
        generate_time += (clock() - generate_start);

        sorting_start = clock();
//        for(int task = 0; task < TASK_SIZE; ++task){
//            SingleUpdateByTask(task);
//        }
		for (int task = 0; task < TASK_SIZE; ++task)
		{
			for (int island = 0; island < ISLAND_NUM; ++island)
			{
				thread_array_island[task*ISLAND_NUM + island] = thread(&Solver::UpdateByTask, this, task, island);
			}
		}
		for (auto& t : thread_array_island)
		{
			t.join();
		}
		sorting_time += (clock() - sorting_start);


//		if (generation % 20 == 0)
//		{
//			float sum = 0.0;
//			for (int task = 0; task < TASK_SIZE; ++task)
//			{
//				for (int island = 0; island < ISLAND_NUM; ++island)
//				{
//					/*cout << "task:  " << task + 1 << "    island:    " << island + 1 << "   :    "
//						<< tasks[task].chromosome[island * ISLAND_SIZE]->phenotype << endl;*/
//					sum += tasks[task].chromosome[island * ISLAND_SIZE]->phenotype;
//				}
//			}
//			cout << "sum of best results:  " << sum << endl;
//		}

		cout << "generation " << generation << " finished ! " << endl;

		//output generation phenotype into files
//        for (int task = 0; task < TASK_SIZE; ++task)
//		{
////			string output_file = "D:\\Visual Stdio\\Island Model\\IslandModel\\IslandModel\\output\\task " + to_string(task) + ".txt";
//			string output_file = "C:\\Users\\user5\\Desktop\\cm\\Island_Model\\Island_Model\\output\\problem_" + to_string(task) + "_"
//                    + to_string(run_id) + ".txt";
//			ofstream outfile;
//			outfile.open(output_file, ios::app);
//			//best individual in each island
//			for (int island = 0; island < ISLAND_NUM; ++island)
//			{
//				outfile << tasks[task].chromosome[island*ISLAND_SIZE]->phenotype << endl;
//			}
//
//			outfile.close();
//		}

	}


	//output time of each part into files
	//string output_file = "C:\\Users\\user5\\Desktop\\cm\\Island_Model\\Island_Model\\time\\" + to_string(run_id) + ".txt";
//	string output_file = file_path + "TASKS_RUN_" + to_string(run_id) + ".txt";
//	ofstream outfile;
//	outfile.open(output_file, ios::app);
//	outfile << "evolve_time : " << evolve_time << endl;
//	outfile << "migration_time : " << migration_time << endl;
//	outfile << "transfer_time : " << transfer_time << endl;
//	outfile << "sorting_time : " << sorting_time << endl;
//	outfile << "evaluate_time : " << evaluate_time << endl;
//	outfile << "generate_time : " << generate_time << endl;
//	outfile << "exgenerate_time : " << exgenerate_time << endl;
//	outfile.close();




    FILE* f;
    string output_file = file_path + "TASKS_RUN_" + to_string(run_id) + ".txt";
    f = fopen((char*)output_file.c_str(),"w");
    if(f == NULL)
    {
        cout << "----------------------------- end --------------------------------" << endl;
        system("pause");
    }
    fprintf(f,"reproduction_time : %f\n",generate_time/second);
    fprintf(f,"sorting_time : %f\n",sorting_time/second);
    fprintf(f,"migration_time : %f\n",migration_time/second);
    fprintf(f,"transfer_time : %f\n",transfer_time/second);

    fclose(f);



	//ending output
	for (int task = 0; task < min(7, TASK_SIZE); task++)
	{
		cout << "task final  " << task + 1 << " : " << tasks[task].chromosome[0]->phenotype << endl;
		printf("lower = %f  upper = %f\n", tasks[task].auxiliary.lower, tasks[task].auxiliary.upper);
		for (int i = 0; i < tasks[task].auxiliary.dimension; ++i) {
			printf("%f ", tasks[task].chromosome[0]->genotype[i]);
		}
		printf("\n");
	}

//	for (int task = 0; task < TASK_SIZE; ++task)
//	{
//		cout << "task best: " << task + 1 << tasks[task].chromosome[0]->phenotype << endl;
//	}

}

void Solver::SingleUpdateByTask(int task)
{

	//sort by every island
	for (int island = 0; island < ISLAND_NUM; island++)
	{
		sort(tasks[task].chromosome + island * ISLAND_SIZE, tasks[task].chromosome + island * ISLAND_SIZE + ISLAND_SIZE, compareFunc);
		//cout << "island " << island + 1 << " : " << tasks[task].chromosome[island*ISLAND_SIZE]->phenotype << endl;
	}

}


double Solver::getGenerateTime()
{
    return this->generate_time;
}

double Solver::getMigrationTime()
{
    return this->migration_time;
}

double Solver::getTransferTime()
{
    return this->transfer_time;
}

double Solver::getSortingTime()
{
    return this->sorting_time;
}

void Solver::migration(int task, int island)
{

	for (int subpopid = 0; subpopid < MIGRATION_SIZE; ++subpopid)
	{
		int offset = ISLAND_SIZE - MIGRATION_SIZE;
		int source = island * ISLAND_SIZE + subpopid;
		int target = (source + ISLAND_SIZE + offset) % POPULATION_SIZE;

		if (tasks[task].chromosome[target]->phenotype - BIAS > tasks[task].chromosome[source]->phenotype)
		{
			tasks[task].chromosome[target]->phenotype = tasks[task].chromosome[source]->phenotype;
//			cpy(tasks[task].chromosome[target]->genotype, tasks[task].chromosome[source]->genotype, tasks[task].auxiliary.dimension);
		}
	}
}



void Solver::transferByM(float* elem, float* M, float* res, int row, int col)
{
	for (int i = 0; i < row; ++i) {
		float temp = 0;
		for (int j = 0; j < col; ++j) {
			int idx = i * (MAX_DIMENSION + 1) + j;
			temp = temp + M[idx] * elem[j];
		}
		res[i] = temp;
	}
}


void Solver::calculateM(float* M, int target, int source)
{
//	mat Q(MAX_DIMENSION + 1, POPULATION_SIZE);
//	mat P(MAX_DIMENSION, POPULATION_SIZE);
//
//
//	for (int gen = 0; gen < MAX_DIMENSION; ++gen)
//	{
//		for (int popu = 0; popu < POPULATION_SIZE; ++popu)
//		{
//			Q(gen, popu) = tasks[target].chromosome[popu]->genotype[gen];
//			P(gen, popu) = tasks[source].chromosome[popu]->genotype[gen];
//		}
//	}

	//add bias to Q
	for (int bias_item = 0; bias_item < POPULATION_SIZE; ++bias_item)
	{
		tasks[target].Q(MAX_DIMENSION, bias_item) = 1;
	}

	mat M_Mat(MAX_DIMENSION, MAX_DIMENSION + 1);

	mat term1 = (tasks[source].P) * (tasks[target].Q.t());
	mat term2 = (tasks[target].Q) * (tasks[target].Q.t());

	//add bias to (Q * Qt)
//	for (int row = 0; row < term2.n_rows; ++row)
//	{
//	    term2(row, row) += BIAS;
//	}
	term2 = term2 + eye(term2.n_rows, term2.n_cols) * BIAS;

	term2 = inv(term2);
	M_Mat = term1 * term2;

	for (int i = 0; i < M_Mat.n_rows; ++i)
	{
		for (int j = 0; j < M_Mat.n_cols; ++j)
		{
			M[i * M_Mat.n_rows + j] = M_Mat.at(i, j);
		}
	}
}


void Solver::explicitTransfer(int target, int island)
{
	float M[MAX_DIMENSION * (MAX_DIMENSION + 1)];
	//float res[MAX_DIMENSION];


	float MAux[MAX_DIMENSION * MAX_DIMENSION] = { 0.0 };

	cpy(MAux, tasks[target].auxiliary.M, MAX_DIMENSION * MAX_DIMENSION);
	float bias = tasks[target].auxiliary.bias;
	int dimension = tasks[target].auxiliary.dimension;

	int source;
	do {
		source = rand() % TASK_SIZE;
	} while (source == target);


	int target_dimension = tasks[target].auxiliary.dimension;
	int source_dimension = tasks[source].auxiliary.dimension;

    this->calculateM(M, target, source);

	for (int popu = island * ISLAND_SIZE; popu < island * ISLAND_SIZE + TRANSFER_SIZE; ++popu)
	{

		int target_popu = popu + (ISLAND_SIZE - TRANSFER_SIZE);
		//				int target_popu = ISLAND_SIZE - (TRANSFER_SIZE - (popu - island * ISLAND_SIZE)) + island * ISLAND_SIZE;
			//cout << "target : " << target_popu << endl;
		float temp[MAX_DIMENSION], temp2[MAX_DIMENSION], genotype_copy[MAX_DIMENSION + 1];
		cpy(genotype_copy, tasks[source].chromosome[popu]->genotype, target_dimension);
		genotype_copy[MAX_DIMENSION] = 1;

		this->transferByM(genotype_copy, M, temp, target_dimension, source_dimension + 1);
		for (int i = 0; i < target_dimension; ++i)
		{
			this->CheckDomain(temp[i], tasks[target].auxiliary.lower, tasks[target].auxiliary.upper);
		}
		transform(temp2, MAux, temp, bias, dimension, target);
		float res = functionChoice(target, temp2, target_dimension);
		if (tasks[target].chromosome[target_popu]->phenotype > res + BIAS) {
			cpy(tasks[target].chromosome[target_popu]->genotype, temp, target_dimension);
			tasks[target].chromosome[target_popu]->phenotype = res;
		}
	}

}


void Solver::implicitTransfer(int target, int island)
{

	float MAux[MAX_DIMENSION * MAX_DIMENSION] = { 0.0 };

    int source;
	do {
		source = rand() % TASK_SIZE;
	} while (source == target);

	cpy(MAux, tasks[target].auxiliary.M, MAX_DIMENSION * MAX_DIMENSION);
	float bias = tasks[target].auxiliary.bias;
	float lower = GEN_LOWER;
	float upper = GEN_UPPER;
	int target_dimension = tasks[target].auxiliary.dimension;
    int source_dimension = tasks[source].auxiliary.dimension;



	float gen_res[TRANSFER_SIZE][MAX_DIMENSION],transfered[TRANSFER_SIZE][MAX_DIMENSION];

	for(int transfer_size = 0; transfer_size < TRANSFER_SIZE; ++transfer_size)
    {
        int pop_id = island * ISLAND_SIZE + transfer_size;
        cpy(transfered[transfer_size],tasks[target].chromosome[pop_id]->genotype,target_dimension);
        cpy(gen_res[transfer_size],tasks[source].chromosome[pop_id]->genotype,source_dimension);
    }



	for (int popu = island * ISLAND_SIZE; popu < island * ISLAND_SIZE + TRANSFER_SIZE; ++popu)
	{
		int target_popu = popu + (ISLAND_SIZE - TRANSFER_SIZE);
		//				int target_popu = ISLAND_SIZE - (TRANSFER_SIZE - (popu - island * ISLAND_SIZE)) + island * ISLAND_SIZE;
			//cout << "target : " << target_popu << endl;
		float temp[MAX_DIMENSION];

		int id = popu - island * ISLAND_SIZE;
		for (int i = 0; i < target_dimension; ++i)
		{
			this->CheckDomain(gen_res[id][i], lower, upper);
		}

		SBX(target,gen_res[id],transfered[id], lower, upper, target_dimension);
		Mutation(gen_res[id],lower,upper);

		transform(temp, MAux, gen_res[id], bias, target_dimension,target);
		float res = functionChoice(target, temp, target_dimension);
		if (tasks[target].chromosome[target_popu]->phenotype > res + BIAS) {
			cpy(tasks[target].chromosome[target_popu]->genotype, gen_res[id], target_dimension);
			tasks[target].chromosome[target_popu]->phenotype = res;
		}
	}

}
