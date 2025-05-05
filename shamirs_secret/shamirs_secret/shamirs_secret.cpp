// shamirs_secret.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


double cal_y_poly(double x, int size_poly, double* poly)
{
	double y = 0;
	for (int i = 0; i < size_poly; i++)
	{
		y += poly[i] * (pow(x, (double)i));
	}
	return y;
}

double* get_shared(int secret, int n, int k) // where n = number of shares and k = minimum shares require to reverse
{
	if (k > n) return NULL;
	double* poly = (double*)malloc(k * sizeof(double));
	poly[0] = secret;
	srand(time(0)); // seeding the rand function with current time
	for (int i = 1; i < k; i++)
	{
		poly[i] = rand() % (secret + 1); // there are many better way for getting random numbers
		//printf("%f\n", poly[i]);
	}

	double* shared = (double*)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++)
	{
		shared[i] = cal_y_poly(i + 1, k, poly);
	}
	free(poly);
	return shared;
}

double* mul_poly(double *x, size_t x_size, double *y, size_t y_size)
{
	double* poly = (double*)malloc(sizeof(double) * (x_size + y_size - 1));
	memset(poly, 0, (x_size + y_size - 1) * sizeof(double));

	for (size_t i = 0; i < x_size; i++)
	{
		for (size_t j = 0; j < y_size; j++)
		{
			poly[i + j] += x[i] * y[j];
		}
	}
	return poly;
}


double* get_LBP(int ref_point_idx, int size_array, double* x_array)
{
	// Lagrange basis polynomials
	// ref_point_idx = ref_point_index , A point for which LBP need to calculate

	double den = 1; //denominator
	double* poly_LBP = &den; // for initial multiplication i.e 1 * (x-a)
	int poly_LBP_size = 1;
	for (int i = 0; i < size_array; i++)
	{
		if (i == ref_point_idx) continue;
		double poly_num[2] = { 0 }; // numerator, this always be like (x-a)
		poly_num[1] = 1;
		poly_num[0] = (-1) * x_array[i];

		double* ptr_temp = mul_poly(poly_num, 2, poly_LBP, poly_LBP_size++);
		if (poly_LBP && poly_LBP != &den)  free(poly_LBP);
		poly_LBP = ptr_temp;

		den *= x_array[ref_point_idx] - x_array[i];
	}

	for (int i = 0; i < size_array; i++)
	{
		poly_LBP[i] = poly_LBP[i] / den; // ignoring decimal precision
	}

	return poly_LBP;
}

double* get_secret_poly(int size_shared, double* x_shared, double* y_shared)
{
	// x_shared, y_shared = point (x,y)

	double* poly = (double*)malloc(sizeof(double) * size_shared); // verify size
	memset(poly, 0, sizeof(double) * size_shared); // zero poly

	for (int i = 0; i < size_shared; i++)
	{
		double x = x_shared[i];
		double y = y_shared[i]; // take one point (x,y)

		double* poly_LBP = get_LBP(i, size_shared, x_shared); // ignoring NULL check

		/*for (int j = 0; j < size_shared; j++)
		{
			printf("%f\n", poly_LBP[j]);
		}*/

		for (int j = 0; j < size_shared; j++)
		{
			poly_LBP[j] = poly_LBP[j] * y;
			poly[j] += poly_LBP[j];
		}

		if(poly_LBP) free(poly_LBP);
	}

	return poly;
}


int main()
{
	int secret = 1234; // the secret
	int n = 6; // secret will be split into 6 shares
	int k = 3; // where any subset of 3 shares is sufficient to reconstruct the secret

	double* shared = get_shared(secret, n, k); // where y = shared[index] and x = index+1

	printf("the shared are :\n");
	for (int i = 0; i < n; i++)
	{
		printf("%d, %f\n", i + 1, shared[i]);
	}

	printf("Now let choose three points 2, 4, 5 to get secret back\n");

	double x_points[3] = { 2,4,5 };
	double y_points[3] = { 0 };
	y_points[0] =  shared[1]; 
	y_points[1] =  shared[3]; 
	y_points[2] =  shared[4]; 

	printf("the points are :\n");
	for (int i = 0; i < 3; i++)
	{
		printf("%f, %f\n", x_points[i], y_points[i]);
	}

	double* secret_poly = get_secret_poly(k, x_points, y_points);

	printf("the secret is %f \n", secret_poly[0]);

	free(shared);
	free(secret_poly); // we can ignore it here
	getc(stdin);
	return 0;
}


