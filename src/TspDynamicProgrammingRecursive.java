import java.util.*;

public class TspDynamicProgrammingRecursive {

    public static double[][] costMatrix;
    public static double[] xCoord;
    public static double[] yCoord;

    private final int N;
    private final int START_NODE;
    private final int FINISHED_STATE;

    private double[][] distance;
    private double minTourCost = Double.POSITIVE_INFINITY;

    private List<Integer> tour = new ArrayList<>();
    private boolean ranSolver = false;

    public static void main(String[] args) {

        // Create adjacency matrix
        int n = 5;
        int[] verticesArray = new int[n];
        for (int i = 0; i < n; i++)
            verticesArray[i] = i;
        costMatrix = generateRandomCircularGraphCostMatrix(5, 100);
        TspDynamicProgrammingRecursive solver = new TspDynamicProgrammingRecursive(costMatrix);

        System.out.println("\nCost Matrix:");
        print2dDoubleArrayPretty(costMatrix);
        System.out.println("\nTesting using circular matrix:");
        System.out.println("The shortest path is: " + solver.getTour());
        System.out.printf("The shortest length is: %.2f\n\n", solver.getTourCost());
        System.out.print("Vertex Coords:  ");
        printIntArrayPretty(verticesArray);
        System.out.print("The x coords are: ");
        printDoubleArrayPretty(xCoord);
        System.out.print("The y coords are: ");
        printDoubleArrayPretty(yCoord);


    }

    public TspDynamicProgrammingRecursive(double[][] distance) {
        this(0, distance);
    }

    public TspDynamicProgrammingRecursive(int startNode, double[][] distance) {

        this.distance = distance;
        N = distance.length;
        START_NODE = startNode;

        // Validate inputs.
        if (N <= 2) throw new IllegalStateException("TSP on 0, 1 or 2 nodes doesn't make sense.");
        if (N != distance[0].length)
            throw new IllegalArgumentException("Matrix must be square (N x N)");
        if (START_NODE < 0 || START_NODE >= N)
            throw new IllegalArgumentException("Starting node must be: 0 <= startNode < N");
        if (N > 32)
            throw new IllegalArgumentException(
                    "Matrix too large! A matrix that size for the DP TSP problem with a time complexity of"
                            + "O(n^2*2^n) requires way too much computation for any modern home computer to handle");

        // The finished state is when the finished state mask has all bits are set to
        // one (meaning all the nodes have been visited).
        FINISHED_STATE = (1 << N) - 1;
    }

    // Returns the optimal tour for the traveling salesman problem.
    public List<Integer> getTour() {
        if (!ranSolver) solve();
        return tour;
    }

    // Returns the minimal tour cost.
    public double getTourCost() {
        if (!ranSolver) solve();
        return minTourCost;
    }

    public void solve() {

        // Run the solver
        int state = 1 << START_NODE;
        Double[][] memo = new Double[N][1 << N];
        Integer[][] prev = new Integer[N][1 << N];
        minTourCost = tsp(START_NODE, state, memo, prev);

        // Regenerate path
        int index = START_NODE;
        while (true) {
            tour.add(index);
            Integer nextIndex = prev[index][state];
            if (nextIndex == null) break;
            int nextState = state | (1 << nextIndex);
            state = nextState;
            index = nextIndex;
        }
        tour.add(START_NODE);
        ranSolver = true;
    }

    private double tsp(int i, int state, Double[][] memo, Integer[][] prev) {

        // Done this tour. Return cost of going back to start node.
        if (state == FINISHED_STATE) return distance[i][START_NODE];

        // Return cached answer if already computed.
        if (memo[i][state] != null) return memo[i][state];

        double minCost = Double.POSITIVE_INFINITY;
        int index = -1;
        for (int next = 0; next < N; next++) {

            // Skip if the next node has already been visited.
            if ((state & (1 << next)) != 0) continue;

            int nextState = state | (1 << next);
            double newCost = distance[i][next] + tsp(next, nextState, memo, prev);
            if (newCost < minCost) {
                minCost = newCost;
                index = next;
            }
        }

        prev[i][state] = index;
        return memo[i][state] = minCost;
    }

    public static double[][] generateRandomCircularGraphCostMatrix(int numVertices, int radius) {
        double[][] costMatrix = new double[numVertices][numVertices];
        double[] x = new double[numVertices];
        double[] y = new double[numVertices];
        double stepAngle = 2 * Math.PI / numVertices;

        // generate circular coordinates
        for (int i = 0; i < numVertices; i++) {
            x[i] = radius * Math.sin(i * stepAngle);
            y[i] = radius * Math.cos(i * stepAngle);
        }
        shuffleArray(x, y);

        // build costMatrix
        for (int i = 0; i < numVertices; i++) {
            costMatrix[i][i] = 0;
            for (int j = i + 1; j < numVertices; j++) {
                costMatrix[i][j] = Math.sqrt((Math.pow((x[i] - x[j]), 2) + Math.pow((y[i] - y[j]), 2)));
                costMatrix[j][i] = costMatrix[i][j];
            }
        }

        // for verification testing
        xCoord = x.clone();
        yCoord = y.clone();
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
    public static double[] swapDoubles(double[] arr, int i, int j) {
        double temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
        return arr;
    }

    public static void printDoubleArrayPretty(double[] arr) {
        for (int i = 0; i < arr.length; i++)
            System.out.printf("%8.2f", arr[i]);
        System.out.println();
    }

    public static void printIntArrayPretty(int[] arr) {
        for (int i = 0; i < arr.length; i++)
            System.out.printf("%8d", arr[i]);
        System.out.println();
    }

    public static void print2dDoubleArrayPretty(double[][] arr) {
        for (double[] row : arr) {
            for (double x : row)
                System.out.printf("%.2f\t", x);
            System.out.println();
        }
    }
}