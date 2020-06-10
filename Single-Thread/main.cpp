//
// Created by c50008023 on 2019/8/5.
//
//#pragma warning(disable:4996)
#include "IM.h"

Solver solver(20);

void CreateTimeDictionaies()
{
	string task_size[3] = { "100", "500", "1000" };
	string transfer_type[2] = { "e", "i" };
	//create dictionaries of time
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			string file_path = "D:\\Visual Stdio\\Island Model\\outputs\\time\\" + task_size[i] + "_tasks_" + transfer_type[j];
			_mkdir((char*)file_path.data());
		}
	}
}

int main()
{
	//CreateTimeDictionaies();
	double startTime = clock();

	string file_path = "D:\\Visual Stdio\\Island Model\\outputs\\time\\100_tasks_e\\";

	solver.generateSolution(file_path);
	double endTime = clock();

	double total_time = endTime - startTime;
	//string output_file = "C:\\Users\\user5\\Desktop\\cm\\Island_Model\\Island_Model\\time\\" + to_string(solver.getRunId()) + ".txt";
	string output_file = file_path + "TASKS_RUN_" + to_string(solver.getRunId()) + ".txt";
	ofstream outfile;
	outfile.open(output_file,ios::app);
	outfile << "total_time : " << total_time << endl;

	outfile.close();


	cout << "Total Time: " << total_time << endl;
	system("pause");
}
