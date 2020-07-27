// ESM 511.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#define _CRT_SECURE_NO_WARNINGS

int main()
{
	int i = 0, j = 0, k = 0, l = 0, ni1 = 0, nj1 = 0, ni2 = 0, nj2 = 0, dir = 0, vacsite = 0, test = 0;
	const int numt2 = 6554;
	const int numvac = 1638;
	const int s1 = 128;
	const int s2 = 128;
	int nn[s1][s2][4][2];
	int sites[numt2][2];
	int sitesvac[numvac][2];
	float f[s1][s2];
	float fi[s1][s2];
	float fKMC[s1][s2];
	float fKMCi[s1][s2];
	const int numtri = 10000;
	float RN = 0, J = 0.8, prob = 0, dE = 0, Esum = 0, Esqrsum = 0, Eavg = 0, Esqravg = 0, num = 0, ftemp = 0, E1 = 0, E2 = 0;
	float Cv = 0;
	float CvKMC = 0;
	float E[numtri];
	float EKMC[numtri];
	float Eav[numtri];
	float EavKMC[numtri];
	FILE *resultsE;
	FILE *resultsEav;
	FILE *resultsEKMC;
	FILE *resultsEavKMC;
	FILE *initialconfig;
	FILE *finalconfig;
	FILE *initialconfigKMC;
	FILE *finalconfigKMC;
	FILE *CvwithCvKMC;

	srand(time(NULL));
	/*Randomly generate values to fill fi*/

	for (i = 0; i < s1; i++)										/*Fill entire array with 1's*/
	{
		for (j = 0; j < s2; j++)
		{
			fi[i][j] = 1;
		}
	}

	for (i = 0; i < numt2; i++)										/*Generate random values for i and j to put in -1's*/
	{
		sites[i][0] = rand() % s1;
		sites[i][1] = rand() % s2;
		for (j = 0; j < i; j++)										/*For each iteration, make sure that combo has not yet been used*/
		{
			if ((sites[i][0] == sites[j][0]) && (sites[i][1] == sites[j][1]))
				i--;
		}
	}


	for (i = 0; i < numt2; i++)										/*Assign the -1's*/
	{
		fi[sites[i][0]][sites[i][1]] = -1;
	}

	for (i = 0; i < s1; i++)										/*Fill the matrix that will be updated in the algorithm*/
	{
		for (j = 0; j < s2; j++)
		{
			f[i][j] = fi[i][j];
		}
	}


	for (i = 0; i < s1; i++)
	{
		for (j = 0; j < s2; j++)
		{
			if (i == 0)												/*row and column of nn above*/
			{
				nn[i][j][0][0] = s1 - 1;
				nn[i][j][0][1] = j;
			}
			else
			{
				nn[i][j][0][0] = i - 1;
				nn[i][j][0][1] = j;
			}
			if (j == s2 - 1)										/*row and column of nn right*/
			{
				nn[i][j][1][0] = i;
				nn[i][j][1][1] = 0;
			}
			else
			{
				nn[i][j][1][0] = i;
				nn[i][j][1][1] = j + 1;
			}
			if (i == s1 - 1)										/*row and column of nn below*/
			{
				nn[i][j][2][0] = 0;
				nn[i][j][2][1] = j;
			}
			else
			{
				nn[i][j][2][0] = i + 1;
				nn[i][j][2][1] = j;
			}
			if (j == 0)												/*row and column of nn left*/
			{
				nn[i][j][3][0] = i;
				nn[i][j][3][1] = s2 - 1;
			}
			else
			{
				nn[i][j][3][0] = i;
				nn[i][j][3][1] = j - 1;
			}

		}
	}
	
	for (k = 0; k < numtri; k++)
	{
		for (l = 0; l < s1*s2; l++)
		{
			ni1 = rand() % s1;
			nj1 = rand() % s2;
			ni2 = rand() % s1;
			nj2 = rand() % s2;
			RN = (float)rand() / (float)RAND_MAX;
			/*Randomly generate values for ni1, nj1, ni2, nj2, and RN*/
			E2 = -1 * J * (f[ni2][nj2] * (f[nn[ni1][nj1][0][0]][nn[ni1][nj1][0][1]] + f[nn[ni1][nj1][1][0]][nn[ni1][nj1][1][1]] + f[nn[ni1][nj1][2][0]][nn[ni1][nj1][2][1]] + f[nn[ni1][nj1][3][0]][nn[ni1][nj1][3][1]]) + f[ni1][nj1] * (f[nn[ni2][nj2][0][0]][nn[ni2][nj2][0][1]] + f[nn[ni2][nj2][1][0]][nn[ni2][nj2][1][1]] + f[nn[ni2][nj2][2][0]][nn[ni2][nj2][2][1]] + f[nn[ni2][nj2][3][0]][nn[ni2][nj2][3][1]]));
			E1 = -1 * J * (f[ni1][nj1] * (f[nn[ni1][nj1][0][0]][nn[ni1][nj1][0][1]] + f[nn[ni1][nj1][1][0]][nn[ni1][nj1][1][1]] + f[nn[ni1][nj1][2][0]][nn[ni1][nj1][2][1]] + f[nn[ni1][nj1][3][0]][nn[ni1][nj1][3][1]]) + f[ni2][nj2] * (f[nn[ni2][nj2][0][0]][nn[ni2][nj2][0][1]] + f[nn[ni2][nj2][1][0]][nn[ni2][nj2][1][1]] + f[nn[ni2][nj2][2][0]][nn[ni2][nj2][2][1]] + f[nn[ni2][nj2][3][0]][nn[ni2][nj2][3][1]]));
			dE = E2 - E1;
			if (dE <= 0)
			{
				ftemp = f[ni1][nj1];
				f[ni1][nj1] = f[ni2][nj2];
				f[ni2][nj2] = ftemp;

			}
			else
			{
				prob = exp(-dE);
				if (RN < prob)
				{
					ftemp = f[ni1][nj1];
					f[ni1][nj1] = f[ni2][nj2];
					f[ni2][nj2] = ftemp;

				}
				else
					continue;
			}
		}

		E[k] = 0;
		for (i = 0; i < s1; i++)
		{
			for (j = 0; j < s2; j++)
			{
				E[k] = E[k] - J * f[i][j] * (f[nn[i][j][0][0]][nn[i][j][0][1]] + f[nn[i][j][1][0]][nn[i][j][1][1]] + f[nn[i][j][2][0]][nn[i][j][2][1]] + f[nn[i][j][3][0]][nn[i][j][3][1]]);
			}
		}

	}
	/*Calculate average of E and average of E^2*/

	Esum = 0;
	Esqrsum = 0;
	l = 1;
	for (k = numtri - 1000; k < numtri; k++)
	{
		Esum = Esum + E[k];
		Esqrsum = Esqrsum + (E[k] * E[k]);
		Eav[k] = Esum / l;
		l++;
	}
	Eavg = Esum / (1000);
	Esqravg = Esqrsum / (1000);
	Cv = Esqravg - (Eavg * Eavg);

	Esum = 0;
	Esqrsum = 0;
	for (k = 0; k < numtri; k++)
	{
		Esum = Esum + E[k];
		Esqrsum = Esqrsum + (E[k] * E[k]);
		Eav[k] = Esum / (k + 1);
	}

	/*KMC----------------------------------------------------------------------------------------------------*/
	for (i = 0; i < s1; i++)
	{
		for (j = 0; j < s2; j++)
		{
			fKMC[i][j] = fi[i][j];
		}
	}

	for (i = 0; i < numvac; i++)										/*Generate random values for i and j to put in vacancies*/
	{
		sitesvac[i][0] = rand() % s1;
		sitesvac[i][1] = rand() % s2;
		for (j = 0; j < i; j++)											/*For each iteration, make sure that combo has not yet been used*/
		{
			if ((sitesvac[i][0] == sitesvac[j][0]) && (sitesvac[i][1] == sitesvac[j][1]))
				i--;
		}
	}


	for (i = 0; i < numvac; i++)										/*Assign the 0's*/
	{
		fKMC[sitesvac[i][0]][sitesvac[i][1]] = 0;
	}

	for (i = 0; i < s1; i++)
	{
		for (j = 0; j < s2; j++)
		{
			fKMCi[i][j] = fKMC[i][j];
		}
	}

	for (k = 0; k < numtri; k++)
	{
		for (l = 0; l < s1*s2; l++)
		{
			/*Randomly generate NEW values for ni1, nj1, and RN, and the direction of attempted movement*/
			ni1 = rand() % s1;
			nj1 = rand() % s2;
			RN = (float)rand() / (float)RAND_MAX;
			dir = rand() % 4;

			i = nn[ni1][nj1][dir][0];
			j = nn[ni1][nj1][dir][1];
			E1 = -J*(fKMC[ni1][nj1] * (fKMC[nn[ni1][nj1][0][0]][nn[ni1][nj1][0][1]] + f[nn[ni1][nj1][1][0]][nn[ni1][nj1][1][1]] + f[nn[ni1][nj1][2][0]][nn[ni1][nj1][2][1]] + f[nn[ni1][nj1][3][0]][nn[ni1][nj1][3][1]]) + fKMC[i][j] * (fKMC[nn[i][j][0][0]][nn[i][j][0][1]] + fKMC[nn[i][j][1][0]][nn[i][j][1][1]] + fKMC[nn[i][j][2][0]][nn[i][j][2][1]] + fKMC[nn[i][j][3][0]][nn[i][j][3][1]]));
			E2 = -J*(fKMC[ni1][nj1] * (fKMC[nn[i][j][0][0]][nn[i][j][0][1]] + fKMC[nn[i][j][1][0]][nn[i][j][1][1]] + fKMC[nn[i][j][2][0]][nn[i][j][2][1]] + fKMC[nn[i][j][3][0]][nn[i][j][3][1]]) + fKMC[ni1][nj1] * (fKMC[nn[i][j][0][0]][nn[ni1][nj1][0][1]] + f[nn[ni1][nj1][1][0]][nn[ni1][nj1][1][1]] + f[nn[ni1][nj1][2][0]][nn[ni1][nj1][2][1]] + f[nn[ni1][nj1][3][0]][nn[ni1][nj1][3][1]]));
			dE = E2 - E1;
			vacsite = fKMC[nn[ni1][nj1][dir][0]][nn[ni1][nj1][dir][1]];
			if (vacsite == 0 || fKMC[ni1][nj1] == 0)
			{
				if (dE <= 0)
				{
					ftemp = fKMC[nn[ni1][nj1][dir][0]][nn[ni1][nj1][dir][1]];
					fKMC[nn[ni1][nj1][dir][0]][nn[ni1][nj1][dir][1]] = fKMC[ni1][nj1];
					fKMC[ni1][nj1] = ftemp;
				}
				else
				{
					prob = exp(-dE);
					if (RN < prob)
					{
						ftemp = fKMC[nn[ni1][nj1][dir][0]][nn[ni1][nj1][dir][1]];
						fKMC[nn[ni1][nj1][dir][0]][nn[ni1][nj1][dir][1]] = fKMC[ni1][nj1];
						fKMC[ni1][nj1] = ftemp;
					}
				}
			}
		}
		EKMC[k] = 0;
		for (i = 0; i < s1; i++)
		{
			for (j = 0; j < s2; j++)
			{
				EKMC[k] = EKMC[k] - J * fKMC[i][j] * (fKMC[nn[i][j][0][0]][nn[i][j][0][1]] + fKMC[nn[i][j][1][0]][nn[i][j][1][1]] + fKMC[nn[i][j][2][0]][nn[i][j][2][1]] + fKMC[nn[i][j][3][0]][nn[i][j][3][1]]);
			}
		}
	}
	/*Calculate average of E and average of E^2*/

	Esum = 0;
	Esqrsum = 0;
	l = 1;
	for (k = numtri - 1000; k < numtri; k++)
	{
		Esum = Esum + EKMC[k];
		Esqrsum = Esqrsum + (EKMC[k] * EKMC[k]);
		EavKMC[k] = Esum / l;
		l++;
	}
	Eavg = Esum / 1000;
	Esqravg = Esqrsum / 1000;
	CvKMC = Esqravg - (Eavg * Eavg);

	Esum = 0;
	Esqrsum = 0;
	for (k = 0; k < numtri; k++)
	{
		Esum = Esum + EKMC[k];
		Esqrsum = Esqrsum + (EKMC[k] * EKMC[k]);
		EavKMC[k] = Esum / (k + 1);
	}
	
	resultsE = fopen("resultsE.txt", "w");
	for (k = 0; k <= numtri; k++)
		fprintf(resultsE, "%f, ", E[k]);
	fclose(resultsE);

	resultsEav = fopen("resultsEav.txt", "w");
	for (k = 0; k <= numtri; k++)
		fprintf(resultsEav, "%f, ", Eav[k]);
	fclose(resultsEav);

	resultsEKMC = fopen("resultsEKMC.txt", "w");
	for (k = 0; k <= numtri; k++)
		fprintf(resultsEKMC, "%f, ", EKMC[k]);
	fclose(resultsEKMC);

	resultsEavKMC = fopen("resultsEavKMC.txt", "w");
	for (k = 0; k <= numtri; k++)
		fprintf(resultsEavKMC, "%f, ", EavKMC[k]);
	fclose(resultsEavKMC);

	initialconfig = fopen("initialconfig.txt", "w");
	for (i = 0; i < s1; i++)
	{
		for (j = 0; j < s2; j++)
		{
			fprintf(initialconfig, "%.1f,", fi[i][j]);
		}
		fprintf(initialconfig, "\n");
	}
	fclose(initialconfig);

	finalconfig = fopen("finalconfig.txt", "w");
	for (i = 0; i < s1; i++)
	{
		for (j = 0; j < s2; j++)
		{
			fprintf(finalconfig, "%.1f,", f[i][j]);
		}
		fprintf(finalconfig, "\n");
	}
	fclose(finalconfig);

	initialconfigKMC = fopen("initialconfigKMC.txt", "w");
	for (i = 0; i < s1; i++)
	{
		for (j = 0; j < s2; j++)
		{
			fprintf(initialconfigKMC, "%.1f,", fKMCi[i][j]);
		}
		fprintf(initialconfigKMC, "\n");
	}
	fclose(initialconfigKMC);

	finalconfigKMC = fopen("finalconfigKMC.txt", "w");
	for (i = 0; i < s1; i++)
	{
		for (j = 0; j < s2; j++)
		{
			fprintf(finalconfigKMC, "%.1f,", fKMC[i][j]);
		}
		fprintf(finalconfigKMC, "\n");
	}
	fclose(finalconfigKMC);
	
	CvwithCvKMC = fopen("Cv.txt", "w");
	fprintf(CvwithCvKMC, "Cv: %f\nCvKMC: %f,", Cv, CvKMC);
	fclose(CvwithCvKMC);
}


