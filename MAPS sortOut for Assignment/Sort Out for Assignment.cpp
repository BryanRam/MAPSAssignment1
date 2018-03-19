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
#include <emmintrin.h>	//for integer intrinsics
#include <fstream>			//for file output
#include <iostream>			//for console output
#include <conio.h>			//for kbhit
#include "hr_time.h"		//for stopwatches
#include <stdio.h>			//for fputs
#include <vector>
#include "SampleSort.h"
#include "RadixSort.h"

using namespace std;

#define MAX_ROWS 2000
#define MAX_COLS 1000
#define MAX_CHARS 6			// numbers are in the range 1- 32,767, so 6 digits is enough.

#define SortedRows "SortedRows.txt"
#define SortedAll  "SortedAll.txt"		// for future use
#define SortedTimes "SortedTimes.txt"


int _data [MAX_ROWS][MAX_COLS], _data2[MAX_ROWS][MAX_COLS];		// 2000 rows of 1000 numbers to sort!
const int rseed = 123;				// arbitrary seed for random number generator - PLEASE DON'T ALTER
									// After sorting data generated with seed 123 these results should be true:
const int	checkBeg = 87,			// at [0][0]
			checkMid = 16440,		// at [MAX_ROWS/2][MAX_COLS/2]
			checkEnd = 32760;		// at [MAX_ROWS-1][MAX_COLS-1]

CStopWatch s1, s2, s3, s4;
SampleSort sampleS;
RadixSort radix;

__declspec(align(16)) int sData[MAX_COLS];
__declspec(align(16)) float pointFour[4] = { 0.4, 0.4, 0.4, 0.4 };
__declspec(align(16)) float five[4] = { 5, 5, 5, 5 };
__declspec(align(16)) int simdAscii[4] = { 48, 48, 48, 48 };
__m128i magicNumber = _mm_setr_epi32(0x66666667, 0x66666667, 0x66666667, 0x66666667);
__m128i zeros = _mm_setr_epi32(0x0, 0x0, 0x0, 0x0);

int threads = 8;

void getData(void);
void sortEachRow(void);
void rSortEachRow(void);
void qSortEachRow(void);
void qSortAll(void);
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
void ItoA_SIMD(int, char *);
int sIndex = _data[0][0];
int mIndex = _data[MAX_ROWS / 2][MAX_COLS / 2];
int eIndex = _data[MAX_ROWS - 1][MAX_COLS - 1];
int start1 = _data[0][0];
int mid1 = (MAX_ROWS / 2)*(MAX_COLS / 2);
int end1 = (MAX_ROWS - 1)*(MAX_COLS - 1);
int* _dataP = &_data2[0][0];



int main(void)
{
	
	omp_set_num_threads(threads);
	getData();

	/*
	thread_count = threads;
	sample_size = 1000;
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

   // #pragma omp parallel for
	for (long thread = 1; thread < thread_count+1; thread++)
	{
		//cout << "Thread: " << omp_get_thread_num() << "\tMax: " << omp_get_max_threads() << endl;
		sampleS.Thread_work(thread);
		//pthread_create(&thread_handles[thread], NULL,
		//	sampleS.Thread_work, (void*)thread);
	}

	sampleS.Print_list(sorted_list, list_size, "Sorted list");
	system("pause");
	*/
	/*
	#pragma omp parallel for
	for (thread = 0; thread < thread_count; thread++)
	{
		pthread_join(thread_handles[thread], NULL);
	}
	*/



//#pragma omp parallel sections
//	{
//		#pragma omp section
//		{
			s1.startTimer();
			//	sortEachRow();
			rSortEachRow();
			//qSortEachRow();
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

			s1.stopTimer();

			displayCheckData();

			s3.startTimer();
			outputDataAsString();
			s3.stopTimer();
		/*}

		#pragma omp section
		{*/
			s2.startTimer();
			//MergeSortP(_dataP, 0, ((MAX_ROWS*MAX_COLS) - 1));
			int n = sizeof(_data2) / sizeof(_data2[0][0]);
			radix.radixsort(_dataP, n);
			//qSortAll();
			s2.stopTimer();

			s4.startTimer();
			outputAllDataAsString();
			s4.stopTimer();
		/*}
	}*/
	outputTimes();

	while (!_kbhit());  //to hold console
}

//*********************************************************************************
void getData()		// Generate the same sequence of 'random' numbers.
{
	srand(123); //random number seed PLEASE DON'T CHANGE!
	for (int i = 0; i<MAX_ROWS; i++)
		for (int j = 0; j < MAX_COLS; j++)
		{
			_data[i][j] = rand(); //RAND_MAX = 32767
			_data2[i][j] = _data[i][j];
		}
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
	
//#pragma omp parallel sections
//	{
//#pragma omp section
//		for (i = 0; i < n1; i++)
//			L[i] = _d[start + i];
//#pragma omp section
//		for (j = 0; j < n2; j++)
//			R[j] = _d[mid + 1 + j];
//	}
	

    #pragma omp parallel for
	for (i = 0; i < n1; i++)
	{
		L[i] = _d[start + i];
	}
	#pragma omp parallel for
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

void rSortEachRow()
{
	#pragma omp parallel for
	for (int i = 0; i<MAX_ROWS; i++)
	{
		int n = sizeof(_data[i]) / sizeof(_data[i][0]);
		radix.radixsort(_data[i], n);
	}
}

//*********************************************************************************
void qSortEachRow()
{
	//cout << "Sorting data...";
    #pragma omp parallel for
	for (int i = 0; i<MAX_ROWS; i++)
	{
		void* _pData = &_data[i][0];
		qsort(_pData, MAX_COLS, sizeof(int), Int_comp);
		/*
		//Use a bubble sort on a row 
		for (int n = MAX_COLS - 1; n >= 0; n--)
		{
			for (int j = 0; j<n; j++)
			{
				if (_data[i][j] > _data[i][j + 1])
				{
					int temp = _data[i][j];
					_data[i][j] = _data[i][j + 1];
					_data[i][j + 1] = temp;
				}
			}
		}
		*/
	}
}

void qSortAll()
{
	qsort(_dataP, MAX_ROWS * MAX_COLS, sizeof(int), Int_comp);
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
	char SIMD_String[MAX_CHARS * 4] = "heree\0";	// can accommodate 4 converted numbers
	string odata;
	cout << "\n\nOutputting data to " << SortedAll << "...";

	for (int i = 0; i<MAX_ROWS; i++) {
		for (int j = 0; j<MAX_COLS; j++) {
			_itoa_s<6>(_data2[i][j], numString, 10);
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

//***********************************************************
// This code should avoid division and use SIMD.

void ItoA_SIMD(int num, char * numStr)
{
/*
//__mi128i pointers to data
__m128i* pidata = (__m128i*) sData;
__m128i* ascii = (__m128i*) simdAscii;
__m128* pFour = (__m128*) pointFour;
__m128* pFive = (__m128*) five;
*/

// *** to be written ***
///*
__asm {
	//!mov		esi, pidata	//point esi to sData
	//!mov     eax, pFour //point eax to pFour
	mov		ebx, numStr		// point EBX to numStr
	//!mov		edi, ascii //point EDI to ascii	
	//!mov		esp, pFive //point ESP to pFive
	movdqa	xmm0, sData	//store first four elements of sData in xmm0
	movdqa	xmm2, xmm0	//make a copy of xmm0
	movdqa  xmm4, pointFour
	movdqa  xmm5, five
	movdqa	xmm6, simdAscii
	vpcmpeqb    xmm2, xmm0, xmm2    //Compare data against zeros
	je		endItoA

	nextSet :
	//jne	nextset
	//mov	[ebx], ASCII		// THEN simply set numStr to ASCII "0"
	//mov	[ebx + 1], NULL	//      add terminating null character
	//jmp	endItoA			//      and end

	// -------- divide num by 10 to get next digit, using a performance trick. ---------------
	//nextset:	mov		eax, 66666667h	// 66666667h = 2^34 / 10
	//mov		eax, 66666667h	// 66666667h = 2^34 / 10
	//pmulhw	xmm0, xmm4
	movdpa	xmm2, xmm0		//make a copy of xmm0
	mulps	xmm0, xmm4		//multiply xmm0 by 2^2/10, store in xmm0
	psrad	xmm0, 2			//xmm0 = num * 2^2/10 * 1/2^2, 
								// cancel out the 2^2, and you get num/10
	movdqa  xmm1, xmm0
	mulps   xmm1, xmm5
	//movdqa  xmm2, [xmm0 + xmm0 * 4]
	addps   xmm1, xmm1
	//lea		ecx,[edx+edx*4]	// ECX = EDX * 5
	//add		ecx,ecx			// ECX = EDX * 10 (could use "sal ecx,1")
	// therefore ECX = (number div 10)*10

	movdqa	xmm3, xmm2	//xmm3 is a copy of A (will become S)
	subps	xmm3, xmm1
	addps	xmm3, xmm6		// add 30h to make EAX the digit's ASCII code

	movdqa[oword ptr numStr + edx], xmm3      //Store xmm3 to numStr array
	 /*
	  x = _mm_packs_epi32(x, x);
	  x = _mm_packus_epi16(x, x);
	  *((int*)array) = _mm_cvtsi128_si32(x);
	 */
	 /*
	 packsswb	xmm4, xmm3
	 packuswb	xmm4, xmm4
	 movd		ebp, xmm4
	 mov			[ebx], ebp
	 */
	 //movdqa		[ebx], xmm4
	add		ebx, 16
	//mov		[ebx], al		// store digit character in "numStr"
	//inc		ebx				// move pointer along string

	//cmp		esi, 0			// if (number div 10) = 0, we've finished
	//jnz		nextset

	mov[ebx], NULL		// so add terminating null character

	//ASM Ver
	//reverse string
	mov		edx, numStr		// the number is in reverse order of digits
	nextChar : dec		ebx				// so we need to reverse the string
	cmp		ebx, edx
	jle		endItoA
	mov		eax, [edx]
	mov		ecx, [ebx]
	mov		[ebx], al
	mov		[edx], cl
	inc		edx
	jmp     nextChar

	ptest	xmm0, zeros
	jne	nextSet

	endItoA :
				pop esi
				pop ebx
		//End ASM ver	
}// end of ASM block
 //*/
}
//***********************************************************