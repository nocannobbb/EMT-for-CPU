#include "parameters.h"

using namespace std;


string nts(int number) {
	char buf[20];
	memset(buf, 0, sizeof(buf));
	sprintf(buf, "%d", number);
	string res = buf;
	return res;
}

void transMatrixRead(string input_file, float* matrix, int dimension) {
	FILE* fp;
	fp = fopen(input_file.c_str(), "r");
	if (fp == NULL) {
		puts("open file failed!");
		cout << input_file << endl;
		//system("pause");
		exit(-1);
	}
	for (int i = 0; i < dimension * dimension; ++i) {
		fscanf(fp, "%f", &matrix[i]);
	}
	fclose(fp);
}
