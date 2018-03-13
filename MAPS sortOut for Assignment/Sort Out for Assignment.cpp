/***********************************************************************
MAPS Assignment 1 - Starting Point for OMP/SIMD Assignment

SortOut - integer sorting and file output program.

The data is held in 2000 rows of 1000 numbers each. Each row is sorted independently,
currently using a simple bubble sort.

Your version should also include sorting all of the data at once, i.e. all 2,000,000 numbers,
and this can be done after the rows have been sorted.

Outputting the strings should employ SIMD to exploit hardware parallelism.

This version includes basic timing information and uses strings to create the file output.

S.Andrews / A.Oram
Revised: A.Oram Feb 2018
************************************************************************
PLEASE ADD YOUR NAME AND STUDENT NUMBER HERE:
Bradley Ramsay
26004360
************************************************************************/
#include <omp.h>
#include <fstream>			//for file output
#include <iostream>			//for console output
#include <conio.h>			//for kbhit
#include "hr_time.h"		//for stopwatches
#include <stdio.h>			//for fputs
#include <vector>
#include "SampleSort.h"

using namespace std;

#define MAX_ROWS 2000
#define MAX_COLS 1000
#define MAX_CHARS 6			// numbers are in the range 1- 32,767, so 6 digits is enough.

#define SortedRows "SortedRows.txt"
#define SortedAll  "SortedAll.txt"		// for future use
#define SortedTimes "SortedTimes.txt"


int _data [MAX_ROWS][MAX_COLS];		// 2000 rows of 1000 numbers to sort!
const int rseed = 123;				// arbitrary seed for random number generator - PLEASE DON'T ALTER
									// After sorting data generated with seed 123 these results should be true:
const int	checkBeg = 87,			// at [0][0]
			checkMid = 16440,		// at [MAX_ROWS/2][MAX_COLS/2]
			checkEnd = 32760;		// at [MAX_ROWS-1][MAX_COLS-1]

CStopWatch s1, s2, s3, s4;
SampleSort sampleS;

int threads = 8;

void getData(void);
void sortEachRow(void);
//void SampleSort(void);
void MergeSort(void);
void MergeSortP(int*, int, int);
void merge(int A[], int start, int mid, int end);
void mergeP(int* d, int start, int mid, int end);
void mergeOuter(int start, int mid, int end);
void displayCheckData(void);
void outputDataAsString(void);
void outputAllDataAsString(void);
void outputTimes(void);
void outputDataAsString(void);
int sIndex = _data[0][0];
int mIndex = _data[MAX_ROWS / 2][MAX_COLS / 2];
int eIndex = _data[MAX_ROWS - 1][MAX_COLS - 1];
int start1 = _data[0][0];
int mid1 = (MAX_ROWS / 2)*(MAX_COLS / 2);
int end1 = (MAX_ROWS - 1)*(MAX_COLS - 1);
int* _dataP = &_data[0][0];



int main(void)
{
	
	omp_set_num_threads(threads);
	getData();

	thread_count = threads;
	//sample_size = ;
	list_size = (MAX_ROWS*MAX_COLS) - 1;

	// Allocate memory for variables
	//thread_handles = malloc(thread_count * sizeof(pthread_t));
	list = (int *) malloc(list_size * sizeof(int));
	tmp_list = (int *) malloc(list_size * sizeof(int));
	sorted_list = (int *) malloc(list_size * sizeof(int));
	sample_keys = (int *) malloc(sample_size * sizeof(int));
	sorted_keys = (int *) malloc(sample_size * sizeof(int));
	splitters =(int *) malloc(thread_count * sizeof(int));

	// One dimensional distribution arrays
	raw_dist = (int *) malloc(thread_count * thread_count * sizeof(int));
	col_dist = (int *) malloc(thread_count * sizeof(int));
	prefix_dist = (int *) malloc(thread_count * thread_count * sizeof(int));
	prefix_col_dist = (int *) malloc(thread_count * sizeof(int));

	list = &_data[0][0];

    #pragma omp parallel for
	for (long thread = 0; thread < thread_count; thread++)
	{
		sampleS.Thread_work(sample_keys);
		//pthread_create(&thread_handles[thread], NULL,
		//	sampleS.Thread_work, (void*)thread);
	}

	sampleS.Print_list(sorted_list, list_size, "Sorted list");
	system("pause");
	/*
	#pragma omp parallel for
	for (thread = 0; thread < thread_count; thread++)
	{
		pthread_join(thread_handles[thread], NULL);
	}
	*/
	s1.startTimer();
	sortEachRow();
//#pragma omp barrier

	/*
	cout << "1st elemdat: " << _data[0][0] << endl;
	cout << "1st elem: " << _dataP[0] << endl;
	cout << "2nd elemdat: " << _data[0][1] << endl;
	cout << "2nd elem: " << _dataP[1] << endl;
	cout << "Middle elemdat: " << _data[MAX_ROWS / 2][MAX_COLS/2] << endl;
	cout << "Middle elem: " << _dataP[1000500] << endl;
	cout << "Last elem: " << _data[MAX_ROWS - 1][MAX_COLS - 1] << endl;
	cout << "Last elem: " << _dataP[(MAX_ROWS*MAX_COLS)-1] << endl;
	*/
	
	//MergeSort();
	
	s1.stopTimer();

	displayCheckData();

	s2.startTimer();
	MergeSortP(_dataP, 0, ((MAX_ROWS*MAX_COLS) - 1));
	s2.stopTimer();

	s3.startTimer();
	outputDataAsString();
	s3.stopTimer();

	s4.startTimer();
	outputAllDataAsString();
	s4.stopTimer();

	outputTimes();

	while (!_kbhit());  //to hold console
}

//*********************************************************************************
void getData()		// Generate the same sequence of 'random' numbers.
{
	srand(123); //random number seed PLEASE DON'T CHANGE!
	for (int i = 0; i<MAX_ROWS; i++)
		for (int j = 0; j<MAX_COLS; j++)
			_data[i][j] = rand(); //RAND_MAX = 32767
}

/*
void SampleSort()
{
	///*To devise a samplesort implementation, one needs to decide on the number of buckets p. When this is done, the actual algorithm operates in three phases:[3]
	1.Sample p−1 elements from the input (the splitters). Sort these; each pair of adjacent splitters then defines a bucket.
	2.Loop over the data, placing each element in the appropriate bucket. (This may mean: send it to a processor, in a multiprocessor system.)
	3.Sort each of the buckets.

	The full sorted output is the concatenation of the buckets.

	//
	cout << "Sorting data...";
	
		for (int j = 0; j <threads; j++)
		{
			int endPos = (j + 1);
			if (endPos > MAX_ROWS) endPos = MAX_ROWS;
			
		}
	
}
*/

void MergeSort()
{
	//#pragma omp parallel for
	for (int i = 0; i < MAX_ROWS; i++)
	{
		MergeSortP(_data[i], 0, MAX_COLS-1);
	}

	/*
	int mid = (MAX_ROWS / 2)*(MAX_COLS / 2);
	int end = (MAX_ROWS - 1)*(MAX_COLS - 1);
	int p = 0, q = (MAX_ROWS / 2)*(MAX_COLS / 2) + 1;

	std::vector< std::vector<int> > Arr( ((MAX_ROWS - 1)*(MAX_COLS - 1 )) + 1);
	int k = 0;

	for (int i = 0; i <= ((MAX_ROWS - 1)*(MAX_COLS - 1)); ++i)
	{
		if (p > mid)
		{
		//	Arr.insert(k++, _data[q++]);
		}
	}
	*/
	/*void merge(int A[], int start, int mid, int end) {
		//stores the starting position of both parts in temporary variables.
		int p = start, q = mid + 1;

		int Arr[end - start + 1];
		int k = 0;

		for (int i = start;i <= end;i++) {
			if (p > mid)      //checks if first part comes to an end or not .
				Arr[k++] = A[q++];

			else if (q > end)   //checks if second part comes to an end or not
				Arr[k++] = A[p++];

			else if (A[p] < A[q])     //checks which part has smaller element.
				Arr[k++] = A[p++];

			else
				Arr[k++] = A[q++];
		}
		for (int p = 0; p< k;p++) {
			// Now the real array has elements in sorted manner including both parts.
			A[start++] = Arr[p];
		}
	}
	*/
	
}

void mergeOuter(int start, int mid, int end)
{
	//stores the starting position of both parts in temporary variables.

	int p = start;
	int q = mid + 1;
	int newIndex = end - start + 1;

	//std::vector< std::vector< int > > item ( 2, std::vector<int> ( 2, 0 ) );
	std::vector< std::vector<int> > Arr(end - start + 1);
	//std::vector<int> Arr(end - start + 1);
	//int Arr[newIndex];
	int k = 0;
	/*
	for (int i = start;i <= end;i++) {
		merge(_data[i], 0, (MAX_COLS / 2) - 1, MAX_COLS - 1);
		if (p > mid)      //checks if first part comes to an end or not .
			Arr[k++] = _data[q++];

		else if (q > end)   //checks if second part comes to an end or not
			Arr[k++] = _data[p++];

		else if (_data[p] < _data[q])     //checks which part has smaller element.
			Arr[k++] = _data[p++];

		else
			Arr[k++] = _data[q++];
	}
	for (int p = 0; p< k;p++) {
		// Now the real array has elements in sorted manner including both parts.
		_data[start++] = Arr[p];
	}
	*/
}

void mergeP(int* _d, int start, int mid, int end)
{
	int i, j, k;
	int p = start;
	int q = mid + 1;
	int newIndex = end - start + 1;
	int n1 = mid - start + 1;
	int n2 = end - mid;

	int Arr[(MAX_ROWS*MAX_COLS) - 1];
	vector<int> L(n1);
	vector<int> R(n2);
	//int L[n1], R[n2];
	
	// Copy data to temp arrays L[] and R[] 
	/*
#pragma omp parallel sections
	{
#pragma omp section
		for (i = 0; i < n1; i++)
			L[i] = _d[start + i];
#pragma omp section
		for (j = 0; j < n2; j++)
			R[j] = _d[mid + 1 + j];
	}
	*/

    //#pragma omp parallel for
	for (i = 0; i < n1; i++)
	{
		L[i] = _d[start + i];
	}
	//#pragma omp parallel for
	for (j = 0; j < n2; j++)
	{
		R[j] = _d[mid + 1 + j];
	}

	// Merge the temp arrays back into arr[l..r]
	i = 0; // Initial index of first subarray
	j = 0; // Initial index of second subarray
	k = start; // Initial index of merged subarray
	while (i < n1 && j < n2)
	{
		if (L[i] <= R[j])
		{
			_d[k] = L[i];
			i++;
		}
		else
		{
			_d[k] = R[j];
			j++;
		}
		k++;
	}

	// Copy the remaining elements of L[], if there are any
	while (i < n1)
	{
		_d[k] = L[i];
		i++;
		k++;
	}

	//Copy the remaining elements of R[], if there are any 
	while (j < n2)
	{
		_d[k] = R[j];
		j++;
		k++;
	}
//}
	/*
	int k = 0;

	for (int i = start;i <= end;i++) {
		if (p > mid)      //checks if first part comes to an end or not .
			Arr[k++] = _d[q++];

		else if (q > end)   //checks if second part comes to an end or not
			Arr[k++] = _d[p++];

		else if (_d[p] < _d[q])     //checks which part has smaller element.
			Arr[k++] = _d[p++];

		else
			Arr[k++] = _d[q++];
	}
	for (int p = 0; p< k;p++) {
		// Now the real array has elements in sorted manner including both parts.
		_d[start++] = Arr[p];
	}
	*/
}

// l is for left index and r is right index of the sub-array of arr to be sorted 
void MergeSortP(int* arr, int l, int r)
{
	if (l < r)
	{
		// Same as (l+r)/2, but avoids overflow for
		// large l and h
		//int m = (l + r) / 2;
		int m = l + (r - l) / 2;

		// Sort first and second halves
		MergeSortP(arr, l, m);
		MergeSortP(arr, m + 1, r);

		mergeP(arr, l, m, r);
	}
}

void merge(int _d[], int start, int mid, int end) {
	//stores the starting position of both parts in temporary variables.

	int p = start;
	int q = mid +1;
	int newIndex = end - start + 1;

	std::vector<int> Arr(MAX_ROWS, MAX_COLS);
	//int Arr[newIndex];
	int k = 0;
	/*
	for (int i = 0; i < MAX_ROWS; i++)
	{
		// at [0][0]
		// at [MAX_ROWS/2][MAX_COLS/2]
		// at [MAX_ROWS-1][MAX_COLS-1]

		for (int n = MAX_COLS - 1; n >= 0; n--)
		{
			if (i >= MAX_ROWS / 2 && n > MAX_COLS / 2)
			{
				Arr[i][k++] = _data[i][q++];
			}
		}
	}*/

	for (int i = start;i <= end;i++) {
		if (p > mid)      //checks if first part comes to an end or not .
			Arr[k++] = _d[q++];

		else if (q > end)   //checks if second part comes to an end or not
			Arr[k++] = _d[p++];

		else if (_d[p] < _d[q])     //checks which part has smaller element.
			Arr[k++] = _d[p++];

		else
			Arr[k++] = _d[q++];
	}
	for (int p = 0; p< k;p++) {
		// Now the real array has elements in sorted manner including both parts.
		_d[start++] = Arr[p];
	}

}

//*********************************************************************************
void sortEachRow()
{
	cout << "Sorting data...";
	for(int i=0; i<MAX_ROWS; i++)
	{	
		//Use a bubble sort on a row 
		for(int n=MAX_COLS-1; n>=0; n--)
		{   for(int j=0; j<n; j++)
			{
                if(_data[i][j] > _data[i][j+1])
				{
                    int temp = _data[i][j];
                    _data[i][j] = _data[i][j+1];
                    _data[i][j+1] = temp;
                }
            }
		}
	}
}


//*********************************************************************************
void displayCheckData()
{
	cout << "\n\ndata[0][0]                   = " << _data[0][0] << "\t" << (_data[0][0] == checkBeg ? "OK" : "BAD");
	cout << "\ndata[MAX_ROWS/2][MAX_COLS/2] = " << _data[MAX_ROWS / 2][MAX_COLS / 2] << "\t" << (_data[MAX_ROWS / 2][MAX_COLS / 2] == checkMid ? "OK" : "BAD");
	cout << "\ndata[MAX_ROWS-1][MAX_COLS-1] = " << _data[MAX_ROWS - 1][MAX_COLS - 1] << "\t" << (_data[MAX_ROWS - 1][MAX_COLS - 1] == checkEnd ? "OK" : "BAD");
}




//*********************************************************************************
void outputTimes()
{
	ofstream soTimes;
	soTimes.open(SortedTimes, ios::app);
	cout << "\n\nTime for sorting all rows   (s) : " << s1.getElapsedTime();
	cout << "\n\nTime for sorting all numbers   (s) : " << s2.getElapsedTime();
	cout << "\nTime for outputting to file SortedRows.txt (s) : " << s3.getElapsedTime();
	cout << "\nTime for outputting to file SortedAll.txt (s) : " << s4.getElapsedTime();
	cout << "\nCombined time               (s) : " << s1.getElapsedTime() + s2.getElapsedTime() + s3.getElapsedTime() + s4.getElapsedTime() << "\n\n\nPress a key to terminate.";

	soTimes << "\n\nTime for sorting all rows : " << s1.getElapsedTime() << " seconds";
	soTimes << "\n\nTime for sorting all numbers : " << s2.getElapsedTime() << " seconds";
	soTimes << "\nTime for outputting to file SortedRows.txt : " << s3.getElapsedTime() << " seconds";
	soTimes << "\nTime for outputting to file SortedAll.txt : " << s4.getElapsedTime() << " seconds";
	soTimes << "\nCombined time               (s) : " << s1.getElapsedTime() + s2.getElapsedTime() + s3.getElapsedTime() + s4.getElapsedTime();
	soTimes << "\n__________________________________________________\n\n";
	soTimes.close();
}


//*********************************************************************************
//Builds a sorted number list as a long string then outputs the whole thing in one big fputs!

void outputDataAsString()
{
	char numString[MAX_CHARS];
	string odata;
	cout << "\n\nOutputting data to " << SortedRows << "...";

	for (int i = 0; i<MAX_ROWS; i++){
		for (int j = 0; j<MAX_COLS; j++){
			_itoa_s<6>(_data[i][j], numString, 10);
			odata += numString;
			odata += "\t";
		}
		odata += "\n";
	}

	FILE * sodata;
	fopen_s(&sodata, SortedRows, "w");
	fputs(odata.c_str(), sodata);
	fclose(sodata);
}
//*********************************************************************************


//*********************************************************************************
//Builds a sorted number list as a long string then outputs the whole thing in one big fputs!

void outputAllDataAsString()
{
	char numString[MAX_CHARS];
	string odata;
	cout << "\n\nOutputting data to " << SortedAll << "...";

	for (int i = 0; i<MAX_ROWS; i++) {
		for (int j = 0; j<MAX_COLS; j++) {
			_itoa_s<6>(_data[i][j], numString, 10);
			odata += numString;
			odata += "\t";
		}
		odata += "\n";
	}

	FILE * sodata;
	fopen_s(&sodata, SortedAll, "w");
	fputs(odata.c_str(), sodata);
	fclose(sodata);
}
//*********************************************************************************
