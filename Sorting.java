import java.util.Random;
import Plotter.Plotter;

public class Sorting {

	final static int SELECTION_VS_QUICK_LENGTH = 12;
	final static int MERGE_VS_QUICK_LENGTH = 15;
	final static int MERGE_VS_QUICK_SORTED_LENGTH = 12;
	final static int SELECT_VS_MERGE_LENGTH = 16;
	final static double T = 600.0;
	/**
	 * Sorts a given array using the quick sort algorithm.
	 * At each stage the pivot is chosen to be the rightmost element of the subarray.
	 * 
	 * Should run in average complexity of O(nlog(n)), and worst case complexity of O(n^2)
	 * 
	 * @param arr - the array to be sorted
	 */
	
	public static void quickSort(double[] arr){
		qSorter(arr, 0, arr.length - 1);
	}
	
	private static void qSorter(double[] a, int start, int end){
		if (start < end){
			int q = partition(a, start, end);
			qSorter(a, start, q - 1);
			qSorter(a, q + 1, end);
		}
	}
	
	private static int partition(double[] a, int start, int end){
		double pivot = a[end];
		int i = start;
		double temp = 0;
		for (int j = start; j < end; j++){
			if (a[j] <= pivot){
				temp = a[i];
				a[i] = a[j];
				a[j] = temp;
				i = i + 1;
			}			
		}
		temp = a[i];
		a[i] = a[end];
		a[end] = temp;
		return i;
	}
	
	/**
	 * Sorts a given array using the merge sort algorithm.
	 * 
	 * Should run in complexity O(nlog(n)) in the worst case.
	 * 
	 * @param arr - the array to be sorted
	 */
	public static void mergeSort(double[] arr){
		int n = arr.length;
		if (n < 2) return;
		int mid = n/2;
		double[] left = new double[mid];
		double[] right = new double[n - mid];
		System.arraycopy(arr, 0, left, 0, left.length);
		System.arraycopy(arr, left.length, right, 0, right.length);
		mergeSort(left);
		mergeSort(right);
		merge(left, right, arr);
	}
	
	private static void merge(double[] left, double[] right, double[] arr){
		int li = 0;
		int ri = 0;
		int arri = 0;		
		while (li < left.length && ri < right.length){
			if (left[li] <= right[ri]){
				arr[arri++] = left[li++];
			}
			else {
				arr[arri++] = right[ri++];
			}
		}
		while (li < left.length) {
			arr[arri++] = left[li++];
		}
		while (ri < right.length) {
			arr[arri++] = right[ri++];
		}
	}

	/**
	 * finds the i'th order statistic of a given array.
	 * 
	 * Should run in complexity O(n) in the worst case.
	 * 
	 * @param arr - the array.
	 * @param i - a number between 0 and arr.length - 1.
	 * @return the number which would be at index i, if the array was to be sorted
	 */
	public static double select (double[] arr, int i){
		return select(arr, 0, arr.length-1, i);
	}

	private static double median(double[] a, int start, int size)
	{
		double[] arr = new double[size];
		System.arraycopy(a, start, arr, 0, size);
		quickSort(arr);
		return arr[(arr.length-1) / 2];
	}

	private static int IndexOf(double[] a, double val)
	{
		for(int i = 0; i < a.length; i++)
		{
			if (a[i] == val){
				return i;
			}
		}
		return -1;
	}

	private static double select(double[] a, int start, int end, int i)
	{
		// if one element
		if (start == end) {
			return a[start];
		}
		int j = 0;
		int medSize = 0;
		if ((end-start)%5 != 0){
			medSize = ((end - start) / 5) + 1;
		}
		else {
			medSize = ((end - start) / 5);
		}
		double[] med = new double[medSize];
		for(int k = 0; k < med.length; k++)
		{
			j = k * 5 + start;
			if ((end - j + 1)  >= 5) {
				med[k] = median(a, j, 5);
			}
			else {
				med[k] = median(a, j, end-j+1);
			}
		}

		double pivot = select(med, 0, (med.length - 1), ((med.length-1)/2));
		int pivotIndex = IndexOf(a, pivot);
		double temp = a[pivotIndex];
		a[pivotIndex] = a[end];
		a[end] = temp;
		int q = partition(a, start, end);
		if (q < i)
			return select(a, q+1, end, i);
		if (q > i)
			return select(a, start, q-1, i);
		return a[i];
	}
	
	/**
	 * Sorts a given array using the selection sort algorithm.
	 * 
	 * Should run in complexity O(n^2) in the worst case.
	 * 
	 * @param arr - the array to be sorted
	 */
	public static void selectionSort(double[] arr){
		int L = arr.length;
		for (int i = 0; i < L; i++) {
			for (int j = i + 1; j < L; j++){
				if (arr[i] > arr[j]){
					double temp = arr[i];
					arr[i] = arr[j];
					arr[j] = temp;
				}
			}
		}
	}
	
	public static void main(String[] args) {
		selectionVsQuick();
		mergeVsQuick();
		mergeVsQuickOnSortedArray();
		selectVsMerge();
	}
	
	/**
	 * Compares the selection sort algorithm against quick sort on random arrays
	 */
	public static void selectionVsQuick(){
		double[] quickTimes = new double[SELECTION_VS_QUICK_LENGTH];
		double[] selectionTimes = new double[SELECTION_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < SELECTION_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumSelection = 0;
			for(int k = 0; k < T; k++){
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				selectionSort(b);
				endTime = System.currentTimeMillis();
				sumSelection += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			selectionTimes[i] = sumSelection/T;
		}
		Plotter.plot("quick sort", quickTimes, "selection sort", selectionTimes);
	}
	
	/**
	 * Compares the merge sort algorithm against quick sort on random arrays
	 */
	public static void mergeVsQuick(){
		double[] quickTimes = new double[MERGE_VS_QUICK_LENGTH];
		double[] mergeTimes = new double[MERGE_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < MERGE_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumMerge = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				mergeSort(b);
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			mergeTimes[i] = sumMerge/T;
		}
		Plotter.plot("quick sort", quickTimes, "merge sort", mergeTimes);
	}

	/**
	 * Compares the merge sort algorithm against quick sort on pre-sorted arrays
	 */
	public static void mergeVsQuickOnSortedArray(){
		double[] quickTimes = new double[MERGE_VS_QUICK_SORTED_LENGTH];
		double[] mergeTimes = new double[MERGE_VS_QUICK_SORTED_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < MERGE_VS_QUICK_SORTED_LENGTH; i++) {
			long sumQuick = 0;
			long sumMerge = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = j;
					b[j] = j;
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				mergeSort(b);
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			mergeTimes[i] = sumMerge/T;
		}
		Plotter.plot("quick sort on sorted array", quickTimes, "merge sort on sorted array", mergeTimes);
	}

	/**
	 * Compares the select algorithm against sorting an array.
	 */
	public static void selectVsMerge(){
		double[] mergeTimes = new double[MERGE_VS_QUICK_LENGTH];
		double[] selectTimes = new double[MERGE_VS_QUICK_LENGTH];
		long startTime, endTime;
		double x;
		Random r = new Random();
		for (int i = 0; i < MERGE_VS_QUICK_LENGTH; i++) {
			long sumMerge = 0;
			long sumSelect = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				int index = (int)(Math.random() * size);
				startTime = System.currentTimeMillis();
				mergeSort(a);
				x = a[index];
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
				startTime = System.currentTimeMillis();
				x = select(b, index);
				endTime = System.currentTimeMillis();
				sumSelect += endTime - startTime;
			}
			mergeTimes[i] = sumMerge/T;
			selectTimes[i] = sumSelect/T;
		}
		Plotter.plot("merge sort and select", mergeTimes, "select", selectTimes);
	}
}
