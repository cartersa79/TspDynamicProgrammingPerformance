import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.io.*;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.function.Function;

public class TspDynamicProgrammingPerformance {

    static ThreadMXBean bean = ManagementFactory.getThreadMXBean();

    /* define constants */
    static long MAXVALUE = 2000000000;
    static long MINVALUE = -2000000000;
    static int numberOfTrials = 1;                   // adjust numberOfTrials and MAXINPUTSIZE based on available
    static int MAXINPUTSIZE = 22;  // heap space limitation
    static int MININPUTSIZE = 3;

    static String ResultsFolderPath = "/home/steve/Results/"; // pathname to results folder
    static FileWriter resultsFile;
    static PrintWriter resultsWriter;

    public static void main(String[] args) {
        // run the whole experiment at least twice, and expect to throw away the data from the earlier runs, before java has fully optimized
        System.out.println("Running first full experiment...");
        runFullExperiment("TspDP-Exp1-ThrowAway.txt");
        System.out.println("Running second full experiment...");
        runFullExperiment("TspDP-Exp2.txt");
        System.out.println("Running third full experiment...");
        runFullExperiment("TspDP-Exp3.txt");

    }

    static void runFullExperiment(String resultsFileName) {
        try {
            resultsFile = new FileWriter(ResultsFolderPath + resultsFileName);
            resultsWriter = new PrintWriter(resultsFile);
        } catch (Exception e) {
            System.out.println("*****!!!!!  Had a problem opening the results file " + ResultsFolderPath + resultsFileName);
            return; // not very foolproof... but we do expect to be able to create/open the file...
        }

        ThreadCpuStopWatch BatchStopwatch = new ThreadCpuStopWatch(); // for timing an entire set of trials
        ThreadCpuStopWatch TrialStopwatch = new ThreadCpuStopWatch(); // for timing an individual trial

        resultsWriter.println("#InputSize    AverageTime"); // # marks a comment in gnuplot data
        resultsWriter.flush();

        /* for each size of input we want to test: in this case starting small and doubling the size each time */
        for (int inputSize = MININPUTSIZE; inputSize <= MAXINPUTSIZE; inputSize++) {
            // progress message...
            System.out.println("Running test for input size " + inputSize + " ... ");

            /* repeat for desired number of trials (for a specific size of input)... */
            long batchElapsedTime = 0;
            // generate a list of random integers in random order to use as test input
            // In this case we're generating one list to use for the entire set of trials (of a given input size)
            //System.out.print("    Generating test data...");
            //long[] testList = createRandomIntegerList(inputSize);
            //System.out.println("...done.");
            //System.out.print("    Running trial batch...");

            /* force garbage collection before each batch of trials run so it is not included in the time */
            System.gc();

            // instead of timing each individual trial, we will time the entire set of trials (for a given input size)
            // and divide by the number of trials -- this reduces the impact of the amount of time it takes to call the
            // stopWatch methods themselves
            //BatchStopwatch.start(); // comment this line if timing trials individually

            // run the trials
            for (long trial = 0; trial < numberOfTrials; trial++) {
                // generate a random list of integers each trial
                double[][] testMatrix = generateRandomEuclideanCostMatrix(inputSize, 100, 100);

                // generate a random key to search in the range of a the min/max numbers in the list
                // long testSearchKey = (long) (0 + Math.random() * (testList[testList.length - 1]));
                /* force garbage collection before each trial run so it is not included in the time */
                // System.gc();

                TrialStopwatch.start(); // *** uncomment this line if timing trials individually
                /* run the function we're testing on the trial input */
                TspDynamicProgrammingRecursive expDP = new TspDynamicProgrammingRecursive(testMatrix);
                expDP.solve();
                batchElapsedTime = batchElapsedTime + TrialStopwatch.elapsedTime(); // *** uncomment this line if timing trials individually
            }
            //batchElapsedTime = BatchStopwatch.elapsedTime(); // *** comment this line if timing trials individually
            double averageTimePerTrialInBatch = (double) batchElapsedTime / (double) numberOfTrials; // calculate the average time per trial in this batch

            /* print data for this size of input */
            resultsWriter.printf("%12d  %15.2f \n", inputSize, averageTimePerTrialInBatch); // might as well make the columns look nice
            resultsWriter.flush();
            System.out.println(" ....done.");
        }
    }

    public static double[][] generateRandomCostMatrix(int numVertices, double maxEdgeLength) {
        double[][] costMatrix = new double[numVertices][numVertices];
        for (int i = 0; i < numVertices; i++) {
            costMatrix[i][i] = 0;
            for (int j = i+1; j < numVertices; j++) {
                costMatrix[i][j] = (Math.random() * ((maxEdgeLength - 1) + 1)) + 1;
                costMatrix[j][i] = costMatrix[i][j];
            }
        }
        return costMatrix;
    }

    public static double[][] generateRandomEuclideanCostMatrix(int numVertices, int xMax, int yMax){
        double [][] costMatrix = new double[numVertices][numVertices];
        int[] x = new int[numVertices];
        int[] y = new int[numVertices];
        // generate random coordinates
        for (int i = 0; i < numVertices; i++){
            x[i] = (int) (Math.random() * ((xMax - 1) + 1)) + 1;
            y[i] = (int) (Math.random() * ((yMax - 1) + 1)) + 1;
        }
        // build costMatrix
        for (int i = 0; i < numVertices; i++){
            costMatrix[i][i] = 0;
            for (int j = i+1; j < numVertices; j++){
                costMatrix[i][j] = Math.sqrt((Math.pow((x[i]-x[j]),2) + Math.pow((y[i]-y[j]),2)));
                costMatrix[j][i] = costMatrix[i][j];
            }
        }
        return costMatrix;
    }

    public static double[][] generateRandomCircularGraphCostMatrix(int numVertices, int radius){
        double [][] costMatrix = new double[numVertices][numVertices];
        double[] x = new double[numVertices];
        double[] y = new double[numVertices];
        double stepAngle = 2*Math.PI/numVertices;

        // generate circular coordinates
        for (int i = 0; i < numVertices; i++){
            x[i] = radius * Math.sin(i * stepAngle);
            y[i] = radius * Math.cos(i * stepAngle);
        }
        shuffleArray(x, y);

        // build costMatrix
        for (int i = 0; i < numVertices; i++){
            costMatrix[i][i] = 0;
            for (int j = i+1; j < numVertices; j++){
                costMatrix[i][j] = Math.sqrt((Math.pow((x[i]-x[j]),2) + Math.pow((y[i]-y[j]),2)));
                costMatrix[j][i] = costMatrix[i][j];
            }
        }

        // for verification testing
        //xCoord = x.clone();
        //yCoord = y.clone();
        //expectedPathLength =

        return costMatrix;
    }

    // shuffle array
    // code modified from: https://www.vogella.com/tutorials/JavaAlgorithmsShuffle/article.html
    public static void shuffleArray(double[] x, double[] y) {
        int n = x.length;
        Random random = new Random();
        random.nextInt();
        for (int i = 0; i < n; i++) {
            int change = i + random.nextInt(n - i);
            swapDoubles(x, i, change);
            swapDoubles(y, i, change);
        }
    }
    // swap i and j
    public static double[] swapDoubles(double[] arr, int i, int j)
    {
        double temp = arr[i] ;
        arr[i] = arr[j];
        arr[j] = temp;
        return arr;
    }

    public static void print2dIntArray(int[][] mat){
        for (int[] row : mat)
            System.out.println(Arrays.toString(row));
    }

    public static void print2dDoubleArray(double[][] arr){
        for (double[] row : arr)
            System.out.println(Arrays.toString(row));
    }

    public static void print2dDoubleArrayPretty(double[][] arr){
        for (double[] row : arr) {
            for (double x : row)
                System.out.printf("%.2f\t", x);
            System.out.println();
        }
    }

    public static void printDoubleArrayPretty(double[] arr){
        for (int i = 0; i < arr.length; i++)
            System.out.printf("%8.2f", arr[i]);
        System.out.println();
    }

    public static void printIntArrayPretty(int[] arr){
        for (int i = 0; i < arr.length; i++)
            System.out.printf("%8d", arr[i]);
        System.out.println();
    }
}