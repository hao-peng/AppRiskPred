#include<stdlib.h>
#include<stdio.h>
#include<time.h>

char head[5000];
char cont[5000];

FILE *trunks[1000];
int main(int argc, char **argv) {
	FILE *fp;
	int fold = 0;
	if (argc != 4) {
		printf("usage: %s input output fold\n", argv[0]);
		return 0;
	}
	fp = fopen(argv[1],"r");
	if (fp == 0) {
		printf("error in open file %s\n", argv[1]);
		return 0;
	}
	
	fold = atoi(argv[3]);
	
	if (fold < 1) {
		printf("wrong number of folds\n");
		return 0;
	}

	for (int i = 0; i < fold; i++) {
		sprintf(head, "%s.%d", argv[2], i);
		trunks[i] = fopen(head, "w");
		if (trunks[i] == 0) {
			printf("error in open file %s\n", head);
			return 0;
		}	
	}

	fgets(head, 5000, fp);
	for(int i = 0; i < fold; i++) {
		fputs(head, trunks[i]);
	}
	srand(time(0));

	while (fgets(cont, 5000, fp) != 0) {
		int i = rand() % fold;
		fputs(cont, trunks[i]);
	}

	fclose(fp);
	for (int i = 0; i < fold; i++) {
		fclose(trunks[i]);
	}
	for (int i = 0; i < fold; i++) {
		sprintf(head, "rm -f %strain", argv[2]);
		system(head);
		for (int j = 0; j < fold; j++) {
			if (j != i) {
				sprintf(head, "cat %s.%d >> %d.a", argv[2], j, i);
				system(head);
			}
		}
		sprintf(head, "cat %s.%d > %d.b", argv[2], i, i);
		system(head);
	}
	for (int i = 0; i < fold; i++) {
		sprintf(head, "rm %s.%i", argv[2], i);
		system(head);
	}
	return 0;
}
