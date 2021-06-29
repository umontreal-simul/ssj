package umontreal.ssj.hups;

import org.junit.jupiter.api.*;
import static org.junit.jupiter.api.Assertions.*;

public class DigitalNetBase2Test {

    protected class DigitalNetBase2ForTest extends DigitalNetBase2 {
		public DigitalNetBase2ForTest (int numRows, int numCols, int dim, int interlacing, int numPoints) {
            this.numRows = numRows;
            this.numCols = numCols;
            this.dim = dim;
            this.interlacing = interlacing;
            this.numPoints = numPoints;
            this.outDigits = MAXBITS;
		}
	}
    
    @Test
    public void testGeneratorMatricesStandardFormat(){
        int numRows = 2;
        int numCols = 2;
        int dim = 2;
        DigitalNetBase2 net = new DigitalNetBase2ForTest(numRows, numCols, dim, 0, 0);
        int[][][] matrices = {{{1, 0}, {0, 1}}, {{0, 1}, {1, 0}}};
        net.generatorMatricesFromStandardFormat(matrices);

        int[] genMatrices = net.getGeneratorMatricesTrans();
        assertEquals(4, genMatrices.length);
        assertArrayEquals(new int[]{1073741824, 536870912, 536870912, 1073741824}, genMatrices);

        int[][][] outputMatrices = net.generatorMatricesToStandardFormat();
        assertArrayEquals(matrices, outputMatrices);
    }

    @Test
	public void testOutputInterlace(){
        int dim = 2;
        int interlacing = 2;
        int numPoints = 1;
        DigitalNetBase2 net = new DigitalNetBase2ForTest(0, 0, dim, interlacing, numPoints);
        int[][] points = {{1073741824, 536870912}};   // 2^30, 2^29
        double[][] result = new double[1][1];
        net.outputInterlace(points, result);
        assertEquals(0.5625, result[0][0], 0.000001);   // 1/2 + 1/2^4
    }

    @Test
    public void testMatrixInterlace(){
        int numRows = 2;
        int numCols = 2;
        int dim = 6;
        int interlacing = 3;
        int numPoints = 1 << numCols;
        DigitalNetBase2 net = new DigitalNetBase2ForTest(numRows, numCols, dim, interlacing, numPoints);
        // Shape 6 x 2 x 2
        int[][][] matrices = {{{1, 0}, {0, 1}}, {{1, 1}, {1, 0}}, {{0, 1}, {1, 1}}, {{0, 1}, {1, 0}}, {{1, 1}, {1, 0}}, {{0, 1}, {1, 1}}};
        net.generatorMatricesFromStandardFormat(matrices);

        DigitalNetBase2 interlacedNet = net.matrixInterlace();

        int[][][] interlacedMatrices = interlacedNet.generatorMatricesToStandardFormat();
        // Shape 2 x 6 x 2
        int[][][] objectiveMatrices = {{{1, 0}, {1, 1}, {0, 1}, {0, 1}, {1, 0}, {1, 1}}, {{0, 1}, {1, 1}, {0, 1}, {1, 0}, {1, 0}, {1, 1}}};
        assertArrayEquals(objectiveMatrices, interlacedMatrices);
    }

}
