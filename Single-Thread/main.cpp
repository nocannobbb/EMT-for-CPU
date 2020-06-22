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
//			string file_path = "D:\\Visual Stdio\\Island Model\\outputs\\time\\" + task_size[i] + "_tasks_" + transfer_type[j];
			string file_path = "C:\\Users\\Administrator\\Desktop\\cm\\single-thread\\outputs\\problem\\" + task_size[i] + "_tasks_" + transfer_type[j];
			_mkdir((char*)file_path.data());
		}
	}
}

int main()
{
//	CreateTimeDictionaies();

    char *directory;
    if((directory = getcwd(NULL, 0)) == NULL) {
        perror("getcwd error");
    }
//    cout << directory << endl;
//    string dir = directory;
//    string replace = "SingleThread";
//    string replace_1 = "single-thread\\SingleThread";
//    solver.input_file = dir.replace(dir.find(replace_1), replace_1.size(), "inputs\\");
//    dir = directory;
//    dir.replace(dir.find(replace), replace.size(), "outputs\\time\\");

    string dir = "C:\\Users\\Administrator\\Desktop\\cm\\single-thread\\outputs\\problem\\";
    solver.input_file = "C:\\Users\\Administrator\\Desktop\\cm\\inputs\\";

	string file_path = dir + nts(TASK_SIZE) + "_tasks_" + (IMPLICIT_TRANSFER?"i":"e")+"\\";

	cout << solver.input_file << endl;
	cout << file_path << endl;
	for(int run_id = 0; run_id < 1; ++run_id)
    {
        double startTime = clock();
        solver.generateSolution(file_path,run_id);
        double endTime = clock();

        //  write into files
        double total_time = endTime - startTime;
        //string output_file = "C:\\Users\\user5\\Desktop\\cm\\Island_Model\\Island_Model\\time\\" + to_string(solver.getRunId()) + ".txt";
//        string output_file = file_path + "TASKS_RUN_" + to_string(run_id) + ".txt";
//        ofstream outfile;
//        outfile.open(output_file,ios::app);
//        outfile << "total_time : " << total_time / CLOCKS_PER_SEC << endl;
//
//        outfile.close();

        cout << "Total Time: " << total_time << endl;
    }

	system("pause");
}
