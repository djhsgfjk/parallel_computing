#include <iostream>
#include <fstream>
#include <omp.h>
#include <random>
#include <string>

using namespace std;

void write_to_file(string s, string filepath)
{
	ofstream out;
	out.open(filepath);
	if (out.is_open())
	{
		out << s << endl;
	}
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


int *full_array(int n) 
{
	int* arr = new int[n];

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
		int id = omp_get_thread_num();
		for (int i = (id * n) / m ; i < ((id + 1) * n) / m; i++)
		{
			sum += arr[i];
		}
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


bool parallel_and_straight_time(int *arr, int n, int m)
{
	clock_t t = clock();
	int sum = parallel_array_sum(arr, n, m);
	t = clock() - t;
	printf("Parallel %f seconds\n", double(t) / CLOCKS_PER_SEC);

	clock_t t1 = clock();
	int sum1 = straight_array_sum(arr, n);
	t1 = clock() - t1;
	printf("Straight %f seconds\n", double(t1) / CLOCKS_PER_SEC);

	printf("sum parallel = %i, sum straight = %i - ", sum, sum1);
	cout << (sum == sum1 ? "true\n" : "false\n");
	return t < t1;

}

int find_M(int threads_num, int k)
{
	bool res = false;
	int n = 1;
	while (!res)
	{
		n *= k;
		int* arr = full_array(n);
		res = parallel_and_straight_time(arr, n, threads_num);
		cout << res << endl;
	}
	printf("M = %i\n", n);
	return n;
}

int hybrid_array_sum(int *arr, int n, int m, string filename)
{
	int M = stoi(read_from_file(filename));
	int sum = 0;
	if (n < M)
		sum = straight_array_sum(arr, n);
	else
		sum = parallel_array_sum(arr, n, m);
	return sum;
}

void Task1(int n, int threads_num)
{
	cout << "Task 1" << endl;
	int *arr = full_array(n);
	int sum = parallel_array_sum(arr, n, threads_num);
	int sum1 = straight_array_sum(arr, n);
	printf("sum parallel = %i, sum straight = %i - ", sum, sum1);
	cout << (sum == sum1 ? "true\n" : "false\n");
	cout << endl;
}

void Task2(int threads_num, int k, int r, string filename)
{
	cout << "Task 2" << endl;
	int M = find_M(threads_num, k);
	for (int i = 1, j = 0; i < r; i++)
	{
		j = find_M(threads_num, k);
		if (j > M)
			M = j;
	}
	cout << "minimum M = " << M << endl;
	write_to_file(to_string(M), filename);
	string text = read_from_file(filename);
	M = stoi(text);
	cout << "M from file - " << M << endl;
	cout << endl;
}

void Task3(int m, string filename)
{
	cout << "Task 3" << endl;
	int M = stoi(read_from_file(filename));
	int n = m * M;
	int *arr = full_array(n);
	cout << (n >= M ? "n >= M" : "n < M") << endl;

	clock_t t = clock();
	int sum = hybrid_array_sum(arr, n, 2, filename);
	t = clock() - t;
	printf("Hybrid %f seconds\n", double(t) / CLOCKS_PER_SEC);

	clock_t t1 = clock();
	int sum1 = straight_array_sum(arr, n);
	t1 = clock() - t1;
	printf("Straight %f seconds\n", double(t1) / CLOCKS_PER_SEC);

	printf("sum parallel = %i, sum straight = %i - ", sum, sum1);
	cout << (sum == sum1 ? "true\n\n" : "false\n\n");
}

void draw_histogram(float* x, int k, float min, float max)
{
	float r = 0.001;

	int equalizer = to_string(k).length();
	int equalizer2 = ceil((max - min) / r);
	
	cout << string((equalizer - 2) / 2, ' ') << "threads num" << string((equalizer - 2) / 2, ' ') << "  speed" << endl;
	for (int i = 1; i <= k; i++)
	{
		int m = 0;
		cout << "     " << i << string(equalizer, ' ') << "     ";
		for (float j = 0.; j <= x[i - 1] - min; j += r) 
		{
			cout << "|";
			m++;
		}
		cout << string(equalizer2 - m + 5, ' ') << " [" << x[i - 1] << " secs]\n";
		if (to_string(i + 1).length() > to_string(i).length())
			equalizer--;
	}
}

void Task4(int k, string filename)
{
	cout << "Task 4" << endl;
	int M = stoi(read_from_file(filename));
	int n = 2 * M;
	int* arr = full_array(n);
	float* x = new float[k];

	clock_t t = clock();
	int sum = straight_array_sum(arr, n);
	x[0] = float(clock() - t) / CLOCKS_PER_SEC;
	float min = x[0];
	float max = x[0];

	for (int i = 2; i <= k; i++)
	{
		t = clock();
		int sum1 = hybrid_array_sum(arr, n, i, filename);
		t = clock() - t;
		x[i - 1] = float(t) / CLOCKS_PER_SEC;
		if (sum1 != sum)
		{
			cout << "FALSE" << endl;
			return;
			//cout << n << ", count - " << count << " " << i << " " << sum << " " << sum1 << " " << sum - sum1;
		}
		if (x[i - 1] < min && x[i - 1] > 0.)
			min = x[i - 1];
		if (x[i - 1] > max)
			max = x[i - 1];
	}
	draw_histogram(x, k, min, max);
}

int main()
{
	string filename = "M.txt";
	Task1(100, 3);
	Task2(2, 10, 10, filename);
	Task3(8, filename);
	Task4(50, filename);

	return 0;
}

