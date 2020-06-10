//
// Created by c50008023 on 2019/8/5.
//
//#pragma warning(disable:4996)
#include "IM.h"

Solver solver;

void CreateTimeDictionaies()
{
	string task_size[3] = { "100", "500", "1000" };
	string transfer_type[2] = { "e", "i" };
	//create dictionaries of time
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			string file_path = "C:\\Users\\Administrator\\Desktop\\cm\\multi-thread\\outputs\\problem\\" + task_size[i] + "_tasks_" + transfer_type[j];
//			string file_path = "D:\\Visual Stdio\\Island Model\\outputs\\time\\" + task_size[i] + "_tasks_" + transfer_type[j];
			_mkdir((char*)file_path.data());
		}
	}
}



int main()
{
//      prepare for file dict and file path
//    CreateTimeDictionaies();
    char *directory;
    if((directory = getcwd(NULL, 0)) == NULL) {
        perror("getcwd error");
    }
//    cout << directory << endl;
    string dir = directory;
    string replace = "ImplicitTransfer";
    string replace_1 = "multi-thread\\ImplicitTransfer";
    solver.input_file = dir.replace(dir.find(replace_1), replace_1.size(), "inputs\\");
    dir = directory;
    dir.replace(dir.find(replace), replace.size(), "outputs\\time\\");

	string file_path = dir + nts(TASK_SIZE) + "_tasks_" + (IMPLICIT_TRANSFER?"i":"e")+"\\";
	cout << file_path << ' ' << solver.input_file << endl;
	for (int run_id = 0; run_id < 5; ++run_id)
	{
		double startTime = clock();

		solver.wholeProcess(file_path,run_id);
		double endTime = clock();

		double total_time = endTime - startTime;

//          write into files
        FILE* f;
        string output_file = file_path + "TASKS_RUN_" + to_string(run_id) + ".txt";
        f = fopen((char*)output_file.c_str(),"a");

        double other_time = total_time - solver.getGenerateTime() - solver.getMigrationTime()
                - solver.getSortingTime() - solver.getTransferTime();

        fprintf(f,"total_time : %f\n",total_time/second);
        fprintf(f,"other_time : %f\n",other_time/second);
        fclose(f);

		cout << "Total Time: " << total_time << endl;
	}

	system("pause");
}
