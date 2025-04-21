#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <string>
#include <fstream>
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
bool isSpaceOnly(string& s);
void makeRandomArray(const string path, vector<string>& vec);
void makeAlreadySortedArray(const string path, vector<string>& vec);
void makeReverseSortedArray(const string path, vector<string>& vec);
void makeNearlySortedArray(const string path, vector<string>& vec);
bool isSorted(const vector<string>& vec);
void insertionSort(vector<string>& a);
void selectionSort(vector<string>& a);
int binarySearch(vector<string>& a, int left, int right, string key);
void binaryInsertionSort(vector<string>& a);
void bubbleSort(vector<string>& a);
void shakerSort(vector<string>& a);
void shellSort(vector<string>& a);
void heapify(vector<string>& a, int n, int i);
void buildHeap(vector<string>& a, int n);
void heapSort(vector<string>& a);
void merge(vector<string>& a, int left, int mid, int right);
void mergeSort(vector<string>& a, int left, int right);
void naturalMergeSort(vector<string>& a);
int medianOfThree(vector<string>& a, int low, int high);
int partition(vector<string>& a, int low, int high);
void quickSort(vector<string>& a, int low, int high);
void radixSort(vector<string>& a);
void counting(vector<string>& a, int pos);
void countingSort(vector<string>& a);
int main()
{
    cout << "Input address of file: ";
    string path;
    cin >> path;
    TimeRuns Time;
    int runs;
    cout << "Input number runs for every algorithm: ";
    cin >> runs;
    cout << endl;

    for (int i = 1; i <= runs; i++) {
        vector<string>vec;
        //makeRandomArray(path, vec); //Tạo mảng
        //makeAlreadySortedArray(path, vec);
        //makeReverseSortedArray(path, vec);
        makeNearlySortedArray(path, vec);

        int n = vec.size();

        vector<string>temp(vec); //Tạo mảng bản sao để thực hiện sắp xếp (chắc chắn các mảng sử dụng cho mỗi thuật toán là giống nhau)
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
    
    ofstream ofs("Result.txt");

    ofs << "RESULTS: " << endl;
    //ofs << "Already Sorted Array:" << endl;
    //ofs << "Reverse Sorted Array:" << endl;
    ofs << "Nearly Sorted Array:" << endl;
    //ofs << "Random Array: " << endl;

    ofs << "Selection Sort: " << fixed << setprecision(2) << Time.SelectionSort / runs << "ms" << endl;
    ofs << "Insertion Sort: " << fixed << setprecision(2) << Time.InsertionSort / runs << "ms" << endl;
    ofs << "Binary Insertion Sort: " << fixed << setprecision(2) << Time.BinaryInsertionSort / runs << "ms" << endl;
    ofs << "Bubble Sort: " << fixed << setprecision(2) << Time.BubbleSort / runs << "ms" << endl;
    ofs << "Shaker Sort: " << fixed << setprecision(2) << Time.ShakerSort / runs << "ms" << endl;
    ofs << "Shell Sort: " << fixed << setprecision(2) << Time.ShellSort / runs << "ms" << endl;
    ofs << "Heap Sort: " << fixed << setprecision(2) << Time.HeapSort / runs << "ms" << endl;
    ofs << "Merge Sort: " << fixed << setprecision(2) << Time.MergeSort / runs << "ms" << endl;
    ofs << "Natural Merge Sort: " << fixed << setprecision(2) << Time.NaturalMergeSort / runs << "ms" << endl;
    ofs << "Quick Sort: " << fixed << setprecision(2) << Time.QuickSort / runs << "ms" << endl;
    ofs << "Function std::sort: " << fixed << setprecision(2) << Time.FuncSort / runs << "ms" << endl;
    ofs << "Radix Sort: " << fixed << setprecision(2) << Time.RadixSort / runs << "ms" << endl;
    ofs << "Counting Sort: " << fixed << setprecision(2) << Time.CountingSort / runs << "ms" << endl;

    ofs.close();

    //cout << "RESULTS: " << endl;
    //cout << "Already Sorted Array:" << endl;
    ////cout << "Reverse Sorted Array:"<< endl;
    ////cout << "Nearly Sorted Array: "<< endl;
    ///*cout << "Random Array: "<< endl;*/

    //cout << "Selection Sort: " << fixed << setprecision(2) << Time.SelectionSort / runs << "ms" << endl;
    //cout << "Insertion Sort: " << fixed << setprecision(2) << Time.InsertionSort / runs << "ms" << endl;
    //cout << "Binary Insertion Sort: " << fixed << setprecision(2) << Time.BinaryInsertionSort / runs << "ms" << endl;
    //cout << "Bubble Sort: " << fixed << setprecision(2) << Time.BubbleSort / runs << "ms" << endl;
    //cout << "Shaker Sort: " << fixed << setprecision(2) << Time.ShakerSort / runs << "ms" << endl;
    //cout << "Shell Sort: " << fixed << setprecision(2) << Time.ShellSort / runs << "ms" << endl;
    //cout << "Heap Sort: " << fixed << setprecision(2) << Time.HeapSort / runs << "ms" << endl;
    //cout << "Merge Sort: " << fixed << setprecision(2) << Time.MergeSort / runs << "ms" << endl;
    //cout << "Natural Merge Sort: " << fixed << setprecision(2) << Time.NaturalMergeSort / runs << "ms" << endl;
    //cout << "Quick Sort: " << fixed << setprecision(2) << Time.QuickSort / runs << "ms" << endl;
    //cout << "Function std::sort: " << fixed << setprecision(2) << Time.FuncSort / runs << "ms" << endl;
    //cout << "Radix Sort: " << fixed << setprecision(2) << Time.RadixSort / runs << "ms" << endl;
    //cout << "Counting Sort: " << fixed << setprecision(2) << Time.CountingSort / runs << "ms" << endl;

    return 0;
}

bool isSpaceOnly(string& s) {
    return all_of(s.begin(), s.end(), [](char ch) {
        return isspace(static_cast<unsigned char>(ch));
        });
}
void makeRandomArray(const string path, vector<string>& vec) {
    makeAlreadySortedArray(path, vec);

    //Tạo các index ngẫu nhiên để hoán đổi vị trí
    int n = vec.size();
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<>indexDis(0, n - 1);

    for (int i = 1; i <= n/2; i++) {
        int idx1 = indexDis(gen), idx2;
        do {
            idx2 = indexDis(gen);
        } while (idx1 == idx2);
        swap(vec[idx1], vec[idx2]);
    }
}
void makeAlreadySortedArray(const string path, vector<string>& vec) {
    ifstream ifs(path);
    if (!ifs.is_open()) {
        cout << "Can't open file" << endl;
        return;
    }
    string temp;
    while (getline(ifs, temp)) {
        if (temp.empty() || isSpaceOnly(temp)) {
            continue;
        }
        vec.push_back(temp);
    }
    ifs.close();
}
void makeReverseSortedArray(const string path, vector<string>& vec) {
    makeAlreadySortedArray(path, vec);
    reverse(vec.begin(), vec.end());
}
void makeNearlySortedArray(const string path, vector<string>& vec) {
    makeAlreadySortedArray(path, vec);

    double sorted;
    cout << "Input Ratio has been sorted: ";
    cin >> sorted;

    //Tạo các index ngẫu nhiên để hoán đổi vị trí
    int n = vec.size();
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<>indexDis(0, n - 1);

    int numSwaps = n * (1.0 - sorted) / 2;
    for (int i = 1; i <= numSwaps; i++) {
        int idx1 = indexDis(gen), idx2;
        do {
            idx2 = indexDis(gen);
        } while (idx1 == idx2);
        swap(vec[idx1], vec[idx2]);
    }
}
bool isSorted(const vector<string>& vec) {
    for (size_t i = 1; i < vec.size(); ++i) {
        if (vec[i] < vec[i - 1]) return false;
    }
    return true;
}
void insertionSort(vector<string>& a) {
    for (int i = 1; i < a.size(); i++) {
        string temp = a[i];
        int j = i - 1;
        while (j >= 0 && a[j] > temp) {
            a[j + 1] = a[j];
            j--;
        }
        a[j + 1] = temp;
    }
}
void selectionSort(vector<string>& a) {
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
int binarySearch(vector<string>& a, int left, int right, string key) {
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
void binaryInsertionSort(vector<string>& a) {
    for (int i = 1; i < a.size(); i++) {
        string key = a[i];
        int pos = binarySearch(a, 0, i - 1, key);
        for (int j = i; j > pos; j--) {
            a[j] = a[j - 1];
        }
        a[pos] = key;
    }
}
void bubbleSort(vector<string>& a) {
    for (int i = 1; i < a.size(); i++) {
        for (int j = a.size() - 1; j >= i; j--) {
            if (a[j - 1] > a[j]) {
                swap(a[j - 1], a[j]);
            }
        }
    }
}
void shakerSort(vector<string>& a) {
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
void shellSort(vector<string>& a) {
    for (int gap = a.size() / 2; gap > 0; gap /= 2) {
        for (int i = gap; i < a.size(); i++) {
            string temp = a[i];
            int j = i - gap;
            while (j >= 0 && a[j] > temp) {
                a[j + gap] = a[j];
                j -= gap;
            }
            a[j + gap] = temp;
        }
    }
}
void heapify(vector<string>& a, int n, int i) {
    string saved = a[i];
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
void buildHeap(vector<string>& a, int n) {
    for (int i = n / 2 - 1; i >= 0; i--) {
        heapify(a, n, i);
    }
}
void heapSort(vector<string>& a) {
    buildHeap(a, a.size());
    for (int i = a.size() - 1; i > 0;i--) {
        swap(a[i], a[0]);
        heapify(a, i, 0);
    }
}
void merge(vector<string>& a, int left, int mid, int right) {
    int n1 = mid - left + 1;
    int n2 = right - mid;
    vector<string>L1(n1), L2(n2);
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
void mergeSort(vector<string>& a, int left, int right) {
    if (left < right) {
        int mid = left + (right - left) / 2;
        mergeSort(a, left, mid);
        mergeSort(a, mid + 1, right);
        merge(a, left, mid, right);
    }
}
void naturalMergeSort(vector<string>& a) {
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
int medianOfThree(vector<string>& a, int low, int high) {
    int mid = low + (high - low) / 2;
    if (a[low] <= a[mid] && a[mid] <= a[high]) return mid;
    else if (a[low] <= a[high] && a[high] <= a[mid]) return high;
    else return low;
}
int partition(vector<string>& a, int low, int high) {
    int pi = medianOfThree(a, low, high);
    string pivot = a[pi];
    swap(a[pi], a[high]);
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
void quickSort(vector<string>& a, int low, int high) {
    if (low < high) {
        int pivot = partition(a, low, high);
        quickSort(a, low, pivot - 1);
        quickSort(a, pivot + 1, high);
    }
}
void radixSort(vector<string>& a) {
    if (a.empty()) return;
    size_t max_len = 0;
    for (const string& s : a) {
        max_len = max(max_len, s.length());
    }
    for (string& s : a) {
        if (s.length() < max_len) {
            s.append(max_len - s.length(), '\0');
        }
    }

    for (int pos = max_len - 1; pos >= 0; --pos) {
        vector<vector<string>> buckets(256);

        for (const string& s : a) {
            unsigned char ch = static_cast<unsigned char>(s[pos]);
            buckets[ch].push_back(s);
        }
        a.clear();
        for (auto& bucket : buckets) {
            for (const string& s : bucket) {
                a.push_back(s);
            }
        }
    }
    for (string& s : a) {
        size_t end = s.find_last_not_of('\0');
        if (end != std::string::npos) {
            s = s.substr(0, end + 1);
        }
        else {
            s.clear();
        }
    }
}
void counting(vector<string>& a, int pos) {
    vector<vector<string>> count(256);
    for (string s : a) {
        unsigned char key = (pos < s.size()) ? s[pos] : 0;
        count[key].push_back(s);
    }
    a.clear();
    for (int i = 0; i < 256; ++i) {
        for (string s : count[i]) {
            a.push_back(s);
        }
    }
}
void countingSort(vector<string>& a) {
    size_t maxLen = 0;
    for (const string& s : a)
        maxLen = max(maxLen, s.length());

    for (int pos = maxLen - 1; pos >= 0; --pos) {
        counting(a, pos);
    }
}