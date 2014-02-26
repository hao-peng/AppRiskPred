/**********
Combine the good apps with the bad apps for testing data
Author: Hao Peng (penghATpurdue.edu)
Last updated time: Feb 25, 2014
**********/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int main(int argc, char **argv) {
	if (argc != 4) {
		printf("Usage: %s good.output bad.output combined\n", argv[0]);
		return 0;
	}
	double p;
	FILE *fp, *cb;
	cb = fopen(argv[3], "w");
	if (cb == 0) {
		printf("cannot open file\n");
		return 0;
	}
	fp = fopen(argv[1], "r");
	if (fp == 0) {
		printf("cannot open file\n");
		return 0;
	}
	while (fscanf(fp, "%lf", &p) != EOF) {
		fprintf(cb, "1\n"); 
	}
	fclose(fp);

	fp = fopen(argv[2], "r");
	if (fp == 0) {
		printf("cannot open file\n");
		return 0;
	}
	while (fscanf(fp, "%lf", &p) != EOF) {
		fprintf(cb, "0\n"); 
	}
	fclose(fp);
	fclose(cb);
	return 0;
}
