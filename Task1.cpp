#include <iostream>
#include <fstream>
#include <omp.h>
#include <random>
#include <string>

using namespace std;


int *full_array(int *arr, int n) 
{
	random_device rd;
	mt19937 mersenne(rd());
	for (int i = 0; i < n; i++)
	{
		arr[i] = 0 + mersenne() % 100;
	}

	return arr;
}

int parallel_array_sum(int* arr, int n, int m)
{
	int sum = 0;
	#pragma omp parallel reduction(+:sum) num_threads(m)
	{
		int threads_num = omp_get_num_threads();
		int id = omp_get_thread_num();
		for (int i = id * n / threads_num; i < (id + 1) * n / threads_num; i++)
		{
			sum += arr[i];
		}
		cout << id << " - " << sum << endl;
	}

	return sum;
}

int straight_array_sum(int* arr, int n)
{
	int sum1 = 0;
	for (int i = 0; i < n; i++)
	{
		sum1 += arr[i];
	}

	return sum1;
}


int parallel_and_straight_time(int *arr, int n, int m)
{
	float t = clock();
	int sum = 0;
	#pragma omp parallel reduction(+:sum) num_threads(m)
	{
		int threads_num = omp_get_num_threads();
		int id = omp_get_thread_num();
		for (int i = id * n / threads_num; i < (id + 1) * n / threads_num; i++)
		{
			sum += arr[i];
		}
	}
	t = clock() - t;
	printf("Parallel %f seconds\n", t);

	float t1 = clock();
	int sum1 = 0;
	for (int i = 0; i < n; i++)
	{
		sum1 += arr[i];
	}
	t1 = clock() - t1;
	printf("Straight %f seconds\n", t1);

	printf("sum parallel = %i, sum = %i - ", sum, sum1);
	cout << (sum == sum1 ? "true\n" : "false\n");
	return (t-t1);

}

int find_M(int m)
{
	float res = 1.;
	int n = 10;
	while (res >= 0.)
	{
		n *= 10;
		int* arr = new int[n];
		arr = full_array(arr, n);
		res = parallel_and_straight_time(arr, n, m);
		cout << res << endl;
	}
	printf("M = %i\n", n);
	return n;
}

int write_to_file(string s, string filepath)
{
	ofstream out;          
	out.open(filepath);
	if (out.is_open())
	{
		out << s << endl;
	}
	return 1;
}

string read_from_file(string filepath)
{
	string line;

	ifstream in(filepath);

	if (in.is_open())
	{
		getline(in, line);
	}
	in.close();

	return line;
}

int main()
{
	cout << "Task 1" << endl;
	int m = 3;
	int n = 10;
	int* arr = new int[n];
	arr = full_array(arr, n);
	int sum = parallel_array_sum(arr, n, m);
	cout << "sum parallel = " << sum << endl;
	int sum1 = straight_array_sum(arr, n);
	cout << "sum straight = " << sum1 << endl;
	cout << endl;

	cout << "Task 2" << endl;
	int M = find_M(3);
	write_to_file(to_string(M), "M.txt");
	string text = read_from_file("M.txt");
	M = stoi(text);
	cout << "M from file - " << M << endl;
	cout << endl;

}

