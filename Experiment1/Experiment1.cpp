#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <iomanip>
using namespace std;

struct TimeRuns {
    double SelectionSort = 0.0;
    double InsertionSort = 0.0;
    double BinaryInsertionSort = 0.0;
    double BubbleSort = 0.0;
    double ShakerSort = 0.0;
    double ShellSort = 0.0;
    double HeapSort = 0.0;
    double MergeSort = 0.0;
    double NaturalMergeSort = 0.0;
    double QuickSort = 0.0;
    double FuncSort = 0.0;
    double RadixSort = 0.0;
    double CountingSort = 0.0;
};
void makeRandomArray(vector<int>& vec, int n, int k);
void makeAlreadySortedArray(vector<int>& vec, int n, int k);
void makeReverseSortedArray(vector<int>& vec, int n, int k);
void makeNearlySortedArray(vector<int>& vec, int n, int k);
bool isSorted(vector<int>& vec);
void swap(int& a, int& b);
void insertionSort(vector<int>& a);
void selectionSort(vector<int>& a);
int binarySearch(vector<int>& a, int left, int right, int key);
void binaryInsertionSort(vector<int>& a);
void bubbleSort(vector<int>& a);
void shakerSort(vector<int>& a);
void shellSort(vector<int>& a);
void heapify(vector<int>& a, int n, int i);
void buildHeap(vector<int>& a, int n);
void heapSort(vector<int>& a);
void merge(vector<int>& a, int left, int mid, int right);
void mergeSort(vector<int>& a, int left, int right);
void naturalMergeSort(vector<int>& a);
int medianOfThree(vector<int>& a, int low, int high);
int partition(vector<int>& a, int low, int high);
void quickSort(vector<int>& a, int low, int high);
int findMax(vector<int>& a);
int findMin(vector<int>& a);
int findDigits(int num);
void radixSort(vector<int>& a);
void countingSort(vector<int>& a);

int main()
{
    int n, k;
    cout << "Input size of array and max value of array: ";
    cin >> n >> k;
    vector<int>vec(n);
    TimeRuns Time;
    int runs;
    cout << "Input number runs for every algorithm: ";
    cin >> runs;
    cout << endl;
    
    //EXPERIMENT 1:
    for (int i = 1; i <= runs; i++) {
        //makeRandomArray(vec, n, k); //Tạo mảng
        makeAlreadySortedArray(vec, n, k);
        //makeReverseSortedArray(vec, n, k);
        //makeNearlySortedArray(vec, n, k);

        vector<int>temp(vec); //Tạo mảng bản sao để thực hiện sắp xếp (chắc chắn các mảng sử dụng cho mỗi thuật toán là giống nhau)
        auto start = chrono::high_resolution_clock::now();
        selectionSort(temp);
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        if (!isSorted(temp)) {
            cout << "Array is not sorted!!!" << endl;
            return 0;
        }
        Time.SelectionSort += duration.count();

        temp = vec; //Tạo lại bản sao cho mảng ban đầu
        start = chrono::high_resolution_clock::now();
        insertionSort(temp);
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        if (!isSorted(temp)) {
            cout << "Array is not sorted!!!" << endl;
            return 0;
        }
        Time.InsertionSort += duration.count();

        temp = vec; //Tạo lại bản sao cho mảng ban đầu
        start = chrono::high_resolution_clock::now();
        binaryInsertionSort(temp);
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        if (!isSorted(temp)) {
            cout << "Array is not sorted!!!" << endl;
            return 0;
        }
        Time.BinaryInsertionSort += duration.count();

        temp = vec; //Tạo lại bản sao cho mảng ban đầu
        start = chrono::high_resolution_clock::now();
        bubbleSort(temp);
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        if (!isSorted(temp)) {
            cout << "Array is not sorted!!!" << endl;
            return 0;
        }
        Time.BubbleSort += duration.count();

        temp = vec; //Tạo lại bản sao cho mảng ban đầu
        start = chrono::high_resolution_clock::now();
        shakerSort(temp);
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        if (!isSorted(temp)) {
            cout << "Array is not sorted!!!" << endl;
            return 0;
        }
        Time.ShakerSort += duration.count();

        temp = vec; //Tạo lại bản sao cho mảng ban đầu
        start = chrono::high_resolution_clock::now();
        shellSort(temp);
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        if (!isSorted(temp)) {
            cout << "Array is not sorted!!!" << endl;
            return 0;
        }
        Time.ShellSort += duration.count();

        temp = vec; //Tạo lại bản sao cho mảng ban đầu
        start = chrono::high_resolution_clock::now();
        heapSort(temp);
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        if (!isSorted(temp)) {
            cout << "Array is not sorted!!!" << endl;
            return 0;
        }
        Time.HeapSort += duration.count();

        temp = vec; //Tạo lại bản sao cho mảng ban đầu
        start = chrono::high_resolution_clock::now();
        mergeSort(temp, 0, n - 1);
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        if (!isSorted(temp)) {
            cout << "Array is not sorted!!!" << endl;
            return 0;
        }
        Time.MergeSort += duration.count();

        temp = vec; //Tạo lại bản sao cho mảng ban đầu
        start = chrono::high_resolution_clock::now();
        naturalMergeSort(temp);
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        if (!isSorted(temp)) {
            cout << "Array is not sorted!!!" << endl;
            return 0;
        }
        Time.NaturalMergeSort += duration.count();

        temp = vec; //Tạo lại bản sao cho mảng ban đầu
        start = chrono::high_resolution_clock::now();
        quickSort(temp, 0, n - 1);
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        if (!isSorted(temp)) {
            cout << "Array is not sorted!!!" << endl;
            return 0;
        }
        Time.QuickSort += duration.count();

        temp = vec; //Tạo lại bản sao cho mảng ban đầu
        start = chrono::high_resolution_clock::now();
        sort(temp.begin(), temp.end());
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        if (!isSorted(temp)) {
            cout << "Array is not sorted!!!" << endl;
            return 0;
        }
        Time.FuncSort += duration.count();

        temp = vec; //Tạo lại bản sao cho mảng ban đầu
        start = chrono::high_resolution_clock::now();
        radixSort(temp);
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        if (!isSorted(temp)) {
            cout << "Array is not sorted!!!" << endl;
            return 0;
        }
        Time.RadixSort += duration.count();

        temp = vec; //Tạo lại bản sao cho mảng ban đầu
        start = chrono::high_resolution_clock::now();
        countingSort(temp);
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        if (!isSorted(temp)) {
            cout << "Array is not sorted!!!" << endl;
            return 0;
        }
        Time.CountingSort += duration.count();
    }

    cout << "RESULTS: " << endl;
    cout << "Already Sorted Array: Size = " << n << endl;
    //cout << "Reverse Sorted Array: Size = " << n << endl;
    //cout << "Nearly Sorted Array: Size = " << n << endl;
    /*cout << "Random Array: Size = " << n << endl;*/

    cout << "Selection Sort: " << fixed << setprecision(2) << Time.SelectionSort / runs << "ms" << endl;
    cout << "Insertion Sort: " << fixed << setprecision(2) << Time.InsertionSort / runs << "ms" << endl;
    cout << "Binary Insertion Sort: " << fixed << setprecision(2) << Time.BinaryInsertionSort / runs << "ms" << endl;
    cout << "Bubble Sort: " << fixed << setprecision(2) << Time.BubbleSort / runs << "ms" << endl;
    cout << "Shaker Sort: " << fixed << setprecision(2) << Time.ShakerSort / runs << "ms" << endl;
    cout << "Shell Sort: " << fixed << setprecision(2) << Time.ShellSort / runs << "ms" << endl;
    cout << "Heap Sort: " << fixed << setprecision(2) << Time.HeapSort / runs << "ms" << endl;
    cout << "Merge Sort: " << fixed << setprecision(2) << Time.MergeSort / runs << "ms" << endl;
    cout << "Natural Merge Sort: " << fixed << setprecision(2) << Time.NaturalMergeSort / runs << "ms" << endl;
    cout << "Quick Sort: " << fixed << setprecision(2) << Time.QuickSort / runs << "ms" << endl;
    cout << "Function std::sort: " << fixed << setprecision(2) << Time.FuncSort / runs << "ms" << endl;
    cout << "Radix Sort: " << fixed << setprecision(2) << Time.RadixSort / runs << "ms" << endl;
    cout << "Counting Sort: " << fixed << setprecision(2) << Time.CountingSort / runs << "ms" << endl;
    

    //EXPERIMENT 2:
    //for (int i = 1; i <= runs; i++) {
    //    makeRandomArray(vec, n, k); //Tạo mảng
    //    //makeAlreadySortedArray(vec, n, k);
    //    //makeReverseSortedArray(vec, n, k);
    //    //makeNearlySortedArray(vec, n, k);

    //    vector<int>temp(vec); //Tạo mảng bản sao để thực hiện sắp xếp (chắc chắn các mảng sử dụng cho mỗi thuật toán là giống nhau)
    //    auto start = chrono::high_resolution_clock::now();
    //    shellSort(temp);
    //    auto end = chrono::high_resolution_clock::now();
    //    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //    if (!isSorted(temp)) {
    //        cout << "Array is not sorted!!!" << endl;
    //        return 0;
    //    }
    //    Time.ShellSort += duration.count();

    //    temp = vec; //Tạo lại bản sao cho mảng ban đầu
    //    start = chrono::high_resolution_clock::now();
    //    heapSort(temp);
    //    end = chrono::high_resolution_clock::now();
    //    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //    if (!isSorted(temp)) {
    //        cout << "Array is not sorted!!!" << endl;
    //        return 0;
    //    }
    //    Time.HeapSort += duration.count();

    //    temp = vec; //Tạo lại bản sao cho mảng ban đầu
    //    start = chrono::high_resolution_clock::now();
    //    mergeSort(temp, 0, n - 1);
    //    end = chrono::high_resolution_clock::now();
    //    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //    if (!isSorted(temp)) {
    //        cout << "Array is not sorted!!!" << endl;
    //        return 0;
    //    }
    //    Time.MergeSort += duration.count();

    //    temp = vec; //Tạo lại bản sao cho mảng ban đầu
    //    start = chrono::high_resolution_clock::now();
    //    naturalMergeSort(temp);
    //    end = chrono::high_resolution_clock::now();
    //    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //    if (!isSorted(temp)) {
    //        cout << "Array is not sorted!!!" << endl;
    //        return 0;
    //    }
    //    Time.NaturalMergeSort += duration.count();

    //    temp = vec; //Tạo lại bản sao cho mảng ban đầu
    //    start = chrono::high_resolution_clock::now();
    //    quickSort(temp, 0, n - 1);
    //    end = chrono::high_resolution_clock::now();
    //    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //    if (!isSorted(temp)) {
    //        cout << "Array is not sorted!!!" << endl;
    //        return 0;
    //    }
    //    Time.QuickSort += duration.count();

    //    temp = vec; //Tạo lại bản sao cho mảng ban đầu
    //    start = chrono::high_resolution_clock::now();
    //    sort(temp.begin(), temp.end());
    //    end = chrono::high_resolution_clock::now();
    //    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //    if (!isSorted(temp)) {
    //        cout << "Array is not sorted!!!" << endl;
    //        return 0;
    //    }
    //    Time.FuncSort += duration.count();

    //    temp = vec; //Tạo lại bản sao cho mảng ban đầu
    //    start = chrono::high_resolution_clock::now();
    //    radixSort(temp);
    //    end = chrono::high_resolution_clock::now();
    //    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //    if (!isSorted(temp)) {
    //        cout << "Array is not sorted!!!" << endl;
    //        return 0;
    //    }
    //    Time.RadixSort += duration.count();

    //    temp = vec; //Tạo lại bản sao cho mảng ban đầu
    //    start = chrono::high_resolution_clock::now();
    //    countingSort(temp);
    //    end = chrono::high_resolution_clock::now();
    //    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //    if (!isSorted(temp)) {
    //        cout << "Array is not sorted!!!" << endl;
    //        return 0;
    //    }
    //    Time.CountingSort += duration.count();
    //}

    //cout << "RESULTS: " << endl;
    ////cout << "Already Sorted Array: Size = " << n << endl;
    ////cout << "Reverse Sorted Array: Size = " << n << endl;
    ////cout << "Nearly Sorted Array: Size = " << n << endl;
    //cout << "Random Array: Size = " << n << endl;

    //cout << "Shell Sort: " << fixed << setprecision(2) << Time.ShellSort / runs << "ms" << endl;
    //cout << "Heap Sort: " << fixed << setprecision(2) << Time.HeapSort / runs << "ms" << endl;
    //cout << "Merge Sort: " << fixed << setprecision(2) << Time.MergeSort / runs << "ms" << endl;
    //cout << "Natural Merge Sort: " << fixed << setprecision(2) << Time.NaturalMergeSort / runs << "ms" << endl;
    //cout << "Quick Sort: " << fixed << setprecision(2) << Time.QuickSort / runs << "ms" << endl;
    //cout << "Function std::sort: " << fixed << setprecision(2) << Time.FuncSort / runs << "ms" << endl;
    //cout << "Radix Sort: " << fixed << setprecision(2) << Time.RadixSort / runs << "ms" << endl;
    //cout << "Counting Sort: " << fixed << setprecision(2) << Time.CountingSort / runs << "ms" << endl;
    //

    return 0;
}

void makeRandomArray(vector<int>& vec, int n, int k) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<>dis(0, k);
    for (int& num : vec) {
        num = dis(gen);
    }
}
void makeAlreadySortedArray(vector<int>& vec, int n, int k) {
    makeRandomArray(vec, n, k);
    sort(vec.begin(), vec.end());
}
void makeReverseSortedArray(vector<int>& vec, int n, int k) {
    makeRandomArray(vec, n, k);
    sort(vec.begin(), vec.end(), greater<int>());
}
void makeNearlySortedArray(vector<int>& vec, int n, int k) {
    makeRandomArray(vec, n, k);
    sort(vec.begin(), vec.end());
    double sorted;
    cout << "Input Ratio has been sorted: ";
    cin >> sorted;

    //Tạo các index ngẫu nhiên để hoán đổi vị trí
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<>indexDis(0, n - 1);

    int numSwaps = n * (1.0 - sorted);
    for (int i = 1; i <= numSwaps; i++) {
        swap(vec[indexDis(gen)], vec[indexDis(gen)]);
    }
}
bool isSorted(vector<int>& vec) {
    for (int i = 1; i < vec.size(); i++) {
        if (vec[i] < vec[i - 1]) return 0;
    }
    return 1;
}
void swap(int& a, int& b) {
    int temp = a;
    a = b;
    b = temp;
}
void insertionSort(vector<int>& a) {
    for (int i = 1; i < a.size(); i++) {
        int temp = a[i];
        int j = i - 1;
        while (j >= 0 && a[j] > temp) {
            a[j + 1] = a[j];
            j--;
        }
        a[j + 1] = temp;
    }
}
void selectionSort(vector<int>& a) {
    for (int i = 0; i < a.size() - 1; i++) {
        int pos = i;
        for (int j = i + 1; j < a.size(); j++) {
            if (a[j] < a[pos]) {
                pos = j;
            }
        }
        swap(a[pos], a[i]);
    }
}
int binarySearch(vector<int>& a, int left, int right, int key) {
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (a[mid] > key) {
            right = mid - 1;
        }
        else {
            left = mid + 1;
        }
    }
    return left;
}
void binaryInsertionSort(vector<int>& a) { 
    for (int i = 1; i < a.size(); i++) {
        int key = a[i];
        int pos = binarySearch(a, 0, i - 1, key);
        for (int j = i; j > pos; j--) {
            a[j] = a[j - 1];
        }
        a[pos] = key;
    }
}
void bubbleSort(vector<int>& a) {
    for (int i = 1; i < a.size(); i++) {
        for (int j = a.size() - 1; j >= i; j--) {
            if (a[j - 1] > a[j]) {
                swap(a[j - 1], a[j]);
            }
        }
    }
}
void shakerSort(vector<int>& a) {
    int left = 0, right = a.size() - 1;
    while (left < right) {
        for (int i = left; i < right; i++) {
            if (a[i + 1] < a[i]) {
                swap(a[i], a[i + 1]);
            }
        }
        right--;
        for (int i = right; i > left; i--) {
            if (a[i] < a[i - 1]) {
                swap(a[i], a[i - 1]);
            }
        }
        left++;
    }
}
void shellSort(vector<int>& a) {
    for (int gap = a.size() / 2; gap > 0; gap /= 2) {
        for (int i = gap; i < a.size(); i++) {
            int temp = a[i];
            int j = i - gap;
            while (j >= 0 && a[j] > temp) {
                a[j + gap] = a[j];
                j -= gap;
            }
            a[j + gap] = temp;
        }
    }
}
void heapify(vector<int>& a, int n, int i) {
    int saved = a[i];
    while (i < n / 2) {
        int child = i * 2 + 1;
        if (child + 1 < n) {
            if (a[child] < a[child + 1]) {
                child++;
            }
        }
        if (saved >= a[child]) break;
        a[i] = a[child];
        i = child;
    }
}
void buildHeap(vector<int>& a, int n) {
    for (int i = n / 2 - 1; i >= 0; i--) {
        heapify(a, n, i);
    }
}
void heapSort(vector<int>& a) {
    buildHeap(a, a.size());
    for (int i = a.size() - 1; i > 0;i--) {
        swap(a[i], a[0]);
        heapify(a, i, 0);
    }
}
void merge(vector<int>& a, int left, int mid, int right) {
    int n1 = mid - left + 1;
    int n2 = right - mid;
    vector<int>L1(n1), L2(n2);
    for (int i = 0; i < n1; i++) {
        L1[i] = a[i + left];
    }
    for (int i = 0; i < n2; i++) {
        L2[i] = a[i + mid + 1];
    }
    int i = 0, j = 0, k = left;
    while (i < n1 && j < n2) {
        if (L1[i] <= L2[j]) {
            a[k] = L1[i];
            i++;
        }
        else {
            a[k] = L2[j];
            j++;
        }
        k++;
    }
    while (i < n1) {
        a[k] = L1[i];
        k++;
        i++;
    }
    while (j < n2) {
        a[k] = L2[j];
        k++;
        j++;
    }
}
void mergeSort(vector<int>& a, int left, int right) {
    if (left < right) {
        int mid = left + (right - left) / 2;
        mergeSort(a, left, mid);
        mergeSort(a, mid + 1, right);
        merge(a, left, mid, right);
    }
}
void naturalMergeSort(vector<int>& a) {
    if (a.size() <= 1) return;
    while (true) {
        vector<pair<int, int>>runs;
        int start = 0;
        while (start < a.size()) {
            int end = start;
            while (end + 1 < a.size() && a[end] <= a[end + 1]) {
                end++;
            }
            runs.push_back({ start, end });
            start = end + 1;
        }
        if (runs.size() <= 1) return;
        for (int i = 0; i + 1 < runs.size(); i += 2) {
            merge(a, runs[i].first, runs[i].second, runs[i + 1].second);
        }
    }
}
int medianOfThree(vector<int>& a, int low, int high) {
    int mid = low + (high - low) / 2;

    if (a[low] > a[mid]) swap(a[low], a[mid]);
    if (a[low] > a[high]) swap(a[low], a[high]);
    if (a[mid] > a[high]) swap(a[mid], a[high]);
    swap(a[mid], a[high]);
    return a[high];
}
int partition(vector<int>& a, int low, int high) {
    int pivot = medianOfThree(a, low, high);
    int i = low - 1;
    for (int j = low; j < high; j++) {
        if (a[j] < pivot) {
            i++;
            swap(a[i], a[j]);
        }
    }
    swap(a[i + 1], a[high]);
    return i + 1;
}
void quickSort(vector<int>& a, int low, int high) {
    if (low < high) {
        int pivot = partition(a, low, high);
        quickSort(a, low, pivot - 1);
        quickSort(a, pivot + 1, high);
    }
}
int findMax(vector<int>& a) {
    int ans = INT_MIN;
    for (int num : a) {
        if (num > ans) {
            ans = num;
        }
    }
    return ans;
}
int findMin(vector<int>& a) {
    int ans = INT_MAX;
    for (int num : a) {
        if (ans > num) ans = num;
    }
    return ans;
}
int findDigits(int num) {
    if (num == 0) return 1;
    int cnt = 0;
    while (num != 0) {
        num /= 10;
        cnt++;
    }
    return cnt;
}
void radixSort(vector<int>& a) {
    int  n = findDigits(findMax(a));
    for (int i = 0; i < n; i++) {
        vector<vector<int>>bin(10);
        for (int num : a) {
            int digit = (num / (int)pow(10, i)) % 10;
            bin[digit].push_back(num);
        }
        int idx = 0;
        for (int j = 0; j <= 9; j++) {
            for (int num : bin[j]) {
                a[idx++] = num;
            }
        }
    }
}
void countingSort(vector<int>& a) {
    int min = findMin(a), max = findMax(a);
    vector<int>freq(max - min + 1, 0);
    for (int num : a) {
        freq[num - min]++;
    }
    int idx = 0;
    for (int i = 0; i < freq.size(); i++) {
        while (freq[i] > 0) {
            a[idx] = min + i;
            freq[i]--;
            idx++;
        }
    }
}