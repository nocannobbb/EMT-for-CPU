#include "parameters.h"



// Sphere [-100, 100]
float func0(float* elem, int& dimension) {

	float sum = 0;
	for (int i = 0; i < dimension; ++i) {
//        cout << "elem:  " << elem[i] << endl;
		sum += elem[i] * elem[i];
	}
//	cout << "func0: " << sum << endl;
	return sum;
}
// Rosenbrock [-50, 50]
float func1(float* elem, int& dimension) {
	float sum = 0;
	for (int i = 0; i < dimension - 1; ++i) {
		sum += 100 * (elem[i] * elem[i] - elem[i + 1]) * (elem[i] * elem[i] - elem[i + 1]) + (elem[i] - 1) * (elem[i] - 1);
	}
	//cout << "func1: " << sum << endl;
	return sum;
}
// Ackley [-50, 50]
float func2(float* elem, int& dimension) {
	float sum = 0;
	float term1 = 0, term2 = 0;
	for (int i = 0; i < dimension; ++i) {
		term1 += elem[i] * elem[i];
		term2 += cos(elem[i] * 2.0 * PI);
	}
	term1 = -0.2 * sqrtf(term1 / (1.0 * dimension));
	term2 = term2 / (1.0 * dimension);
	sum = -20 * expf(term1) - expf(term2) + 22.71828;
	//cout << "func2: " << sum << endl;
	return sum;
}
// Rastrgin [-50, 50]
float func3(float* elem, int& dimension) {
	float sum = 0;
	for (int i = 0; i < dimension; ++i) {
		sum += elem[i] * elem[i] + 10.0 - 10 * cos(elem[i] * 2.0 * PI);
	}
	//cout << "func3: " << sum << endl;
	return sum;
}
// Griewank [-100, 100]
float func4(float* elem, int& dimension) {
	float sum = 1;
	float term1 = 0, term2 = 1.0;
	for (int i = 0; i < dimension; ++i) {
		term1 += elem[i] * elem[i];
		term2 *= cosf(elem[i] / sqrtf(i + 1.0));
	}
	sum += term1 / 4000.0 - term2;
	//cout << "func4: " << sum << endl;
	return sum;
}
// Weierstrass [-0.5, 0.5]
float func5(float* elem, int& dimension) {
	float sum = 0, term = 0;
	float a = 0.5, b = 3, k_max = 20;

	for (int i = 0; i < dimension; ++i) {
		//cout << "elem:  " << elem[i] << endl;
		float temp = 0;
		for (int j = 0; j < k_max; ++j) {
			temp += powf(a, j) * cos(2.0 * PI * powf(b, j) * (elem[i] + 0.5));
		}
		sum += temp;
	}
	for (int i = 0; i < k_max; ++i) {
		term += powf(a, i) * cos(powf(b, i) * PI);
	}
	sum = sum - dimension * term;
	//cout << "func5: " << sum << endl;
	return sum;
}

// Schwefel [-500, 500]
float func6(float* elem, int& dimension) {
	float sum = 418.9829 * dimension;
	for (int i = 0; i < dimension; ++i) {
//        cout << "elem:  " << elem[i] << endl;
		sum -= elem[i] * sinf(sqrtf(fabsf(elem[i])));
	}
//	cout << "func6: " << sum << endl;
	return sum;
}
float functionChoice(int& task_id, float* elem, int& dimension) {

	switch (task_id % 7)
	{
	case 0:
		return func0(elem, dimension);
	case 1:
		return func1(elem, dimension);
	case 2:
		return func2(elem, dimension);
	case 3:
		return func3(elem, dimension);
	case 4:
		return func4(elem, dimension);
	case 5:
		return func5(elem, dimension);
	case 6:
		return func6(elem, dimension);
	default:
		return func0(elem, dimension);
		//            printf("bugs in functionChoice!!!\n");
	}
}

void cpy(float* target, float* source, int dimension)
{
	for (int i = 0; i < dimension; i++)
	{
		target[i] = source[i];
	}
}

 //To transform the elem by vector translation and rotation. The save of M is same with C & C++ mode.
//void transform(float* M, float* elem, float bias, int dimension) {
//	float res[MAX_DIMENSION];
//	for (int i = 0; i < dimension; ++i) {
//		float sum = 0;
//		for (int j = 0; j < dimension; ++j) {
//			sum += M[i * dimension + j] * (elem[j] + bias);
//		}
//		res[i] = sum;
//	}
//	cpy(elem, res, dimension);
//}

//void evaluateTasks(Task* tasks, int offset) {
//	int task_id = blockIdx.x;
//	int island_id = blockIdx.y;
//	int subpop_id = threadIdx.x;
//	int pop_id = subpop_id + island_id * ISLAND_SIZE;
//	int idx = (island_id + task_id * ISLAND_NUM) * ISLAND_SIZE + subpop_id;
//
//
//	int dimension = tasks[task_id].auxiliary.dimension;
//	float bias = tasks[task_id].auxiliary.bias;
//
//	float elem[ISLAND_SIZE][MAX_DIMENSION];
//	float M[MAX_DIMENSION * MAX_DIMENSION];
//	cpy(elem[subpop_id], tasks[task_id].chromosome[pop_id]->genotype, dimension);
//
//	int thread_size = blockDim.x;
//	for (int i = subpop_id; i < dimension * dimension; i += thread_size) {
//		M[i] = tasks[task_id].auxiliary.M[i];
//	}
//	__syncthreads();
//
//	//    // debug
//	//    __syncthreads();
//	//    if(idx == 0)
//	//    {
//	//        print(elem[subpop_id], dimension);
//	//        print(M, dimension*dimension);
//	//    }
//	//    __syncthreads();
//
//
//
//	transform(M, elem[subpop_id], bias, dimension);
//
//	//    // debug
//	//    __syncthreads();
//	//    if(idx == 0)
//	//    {
//	//        print(M, dimension*dimension);
//	//        print(elem[subpop_id], dimension);
//	//    }
//	//    __syncthreads();
//
//
//
//
//	float res = functionChoice(task_id, elem[subpop_id], dimension);
//	tasks[task_id].chromosome[pop_id]->phenotype = res;
//
//	__syncthreads();
//	if (idx == 0 && !CLEAN) {
//		printf("InitTasks finished\n");
//	}
//}


