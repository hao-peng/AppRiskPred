/**********
Change the format of file from cvx to dat for training and testing

Format of training data:
First line: M D. M is the number of categories and D is the number of permissions
Then M groups of data follows.
For group i,
First line: N_i. N_i is the number of apps in category i
Then N_i lines follows, which contains D binary numbers indicating wheather the app uses the corresponding permission

Format of test data:
First line: Ns D. Ns is the number of test data and D is the number of permissions.
Then Ns lines follows, which contains D binaray numbers indicating weather the app uses the corresponding permission

Author: Hao Peng (penghATpurdue.edu)
Last updated time: Feb 25, 2014
**********/
#include<stdio.h>
#include<string.h>
#include<string>
#include<stdlib.h>
#include<map>


using namespace std;

map<string, int> m;

int numCat = 0;
char buf[5000];
int mode = 0;
FILE *outs[100];
int counts[100];
int l;

int main(int argc, char **argv) {
	if (argc != 4) {
		printf("usage: changeFormat input output mode\n");
		printf("mode: 0 for training, 1 for testing\n");
		return 0;
	}
	mode = atoi(argv[3]);
	FILE *in, *out;
	in = fopen(argv[1], "r");
	if (in == 0) {
		printf("error\n");
		return 0;
	}

	if (mode != 0) {
		char filename[1000];
		sprintf(filename, "%s.tmp", argv[2]);
		out = fopen(filename, "w");
	}
	memset(counts, 0, sizeof(counts));
	//skip the first line
	fgets(buf, 5000, in);
	while(fgets(buf, 5000,in) != 0) {
		char *p;
		p = strtok(buf, ",");

		if (mode == 0) {
			if (m.find(p) == m.end()) {
				m[p] = numCat;
				char filename[100];
				sprintf(filename, "%s.%d", argv[2], numCat);
				outs[numCat] = fopen(filename, "w");
				if (outs[numCat] == 0) {
					printf("error\n");
					return 0;
				}
				numCat++;
			}
			int index = m[p];
			//outs[index];
			l = 0;
			if ((p = strtok(NULL, ",")) != 0) {
				while((p = strtok(NULL, ",")) != NULL) {
					l++;
					fprintf(outs[index], "%d ", atoi(p));
				}
				fprintf(outs[index], "\n");
				counts[index]++;
			}
		} else if (mode == 2) {
			if (m.find(p) == m.end()) {
				m[p] = numCat;
				numCat++;
			}
			int index = m[p];
			fprintf(out, "%d ", index);
			l = 0;
			if ((p = strtok(NULL, ",")) != 0) {
				while((p = strtok(NULL, ",")) != NULL) {
					l++;
					fprintf(out, "%d ", atoi(p));
				}
				fprintf(out, "\n");
				counts[0]++;
			}
		} else {
			l = 0;
			if ((p = strtok(NULL, ",")) != 0) {
				while((p = strtok(NULL, ",")) != NULL) {
					l++;
					fprintf(out, "%d ", atoi(p));
				}
				fprintf(out, "\n");
				counts[0]++;
			}
		}
	}
	fclose(in);

	if (mode == 0) {
		for (int i = 0; i < numCat; i++) {
			fclose(outs[i]);
		}
		char command[1000];
		sprintf(command, "echo \"%d %d\" > %s", numCat, l, argv[2]);
		system(command);
		for (int i = 0; i < numCat; i++) {
			sprintf(command, "echo \"%d\" >> %s", counts[i], argv[2]);
			system(command);
			sprintf(command, "cat %s.%d >> %s", argv[2], i, argv[2]);
			system(command);
			sprintf(command, "rm %s.%d", argv[2], i);
			system(command);
		}
	} else {
		fclose(out);
		char command[1000];
		sprintf(command, "echo \"%d %d\" > %s", counts[0], l, argv[2]);
		system(command);
		sprintf(command, "cat %s.tmp >> %s", argv[2], argv[2]);
		system(command);
		sprintf(command, "rm %s.tmp", argv[2]);
		system(command);
	}

	

	return 0;
}
