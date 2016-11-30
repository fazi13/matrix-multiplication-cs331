import java.util.Random;

public class Multiplication {

	public static void main(String[] args) {
		int n = 1024;
		double iterations = 1.0;
		//start with n x n matrix and increase n by 2 fold each iteration
		
		long startTime = 0;
		long stopTime = 0;
		long time1 = 0;
		long time2 = 0;
		long time3 = 0;
		
		int f = 0;
		while(f < iterations){			
			int[][] a = new int[n][n];
			int[][] b = new int[n][n];
			Random rand = new Random();
			//fill matrix with random int from 0-99
			for(int i = 0; i < a.length; i++){
				for(int j = 0; j < a.length; j++){
					int x = rand.nextInt(100);
					a[i][j] = x;
					int y = rand.nextInt(100);
					b[i][j] = y;
				}			
			}
			int[][] c = new int[0][0]; 
			
			//calculate time for each iteration of each multiplication method
			startTime = System.currentTimeMillis();
			c = multiply(a, b);
			stopTime = System.currentTimeMillis();
			time1 += (stopTime - startTime);

			startTime = System.currentTimeMillis();
			c = dcMultiply(a, b);
			stopTime = System.currentTimeMillis();
			time2 += (stopTime - startTime);

			startTime = System.currentTimeMillis();
			c = sMultiply(a, b);
			stopTime = System.currentTimeMillis();
			time3 += (stopTime - startTime);

			
			f++;
		}
		
		System.out.println("Standard Multiplication time: " + Long.toString(time1));
		double avg1 = time1 / iterations; //div by total iterations
		System.out.println("Average time per instance: " + Double.toString(avg1));
		System.out.println("Divide & Conquer Multiplication time: " + Long.toString(time2));
		double avg2 = time2 / iterations; //div by total iterations
		System.out.println("Average time per instance: " + Double.toString(avg2));
		System.out.println("Strassen Multiplication time: " + Long.toString(time3));
		double avg3 = time3 / iterations; //div by total iterations
		System.out.println("Average time per instance: " + Double.toString(avg3));
	}
	
	//standard matrix multiplication
	private static int[][] multiply(int[][] x, int[][] y){
		int[][] z = new int[x.length][x.length];
		for(int i = 0; i < x.length; i++){
			for(int j = 0; j < x.length; j++){
				int c = 0;
				for(int k = 0; k < x.length; k++){
					c += x[i][k] * y[k][j];
				}
				z[i][j] = c;
			}
		}
		return z;
	}
	
	//divide and conquer
	private static int[][] dcMultiply(int[][] x, int[][] y){
	    return matrixMult(x, y, 0, 0, 0, 0, x.length);
	}

	private static int[][] matrixMult(int[][] x, int[][] y, int rowX, int colX, int rowY, int colY, int n){
	    int[][] z = new int[n][n];
	    //base case
	    if(n == 1)
	        z[0][0] = x[rowX][colX] * y[rowY][colY];
	    else {
	    	//split into n/2 sections
	        int m = n/2;
	        //C calculation
	         matrixSum(z, matrixMult(x, y, rowX, colX, rowY, colY, m), matrixMult(x, y, rowX, colX + m, rowY + m, colY, m), 0, 0);
	         matrixSum(z, matrixMult(x, y, rowX, colX, rowY, colY + m, m), matrixMult(x, y, rowX, colX + m, rowY + m, colY + m, m), 0, m);
	         matrixSum(z, matrixMult(x, y, rowX + m, colX, rowY, colY, m), matrixMult(x, y, rowX + m, colX + m, rowY + m, colY, m), m, 0);
	         matrixSum(z, matrixMult(x, y, rowX + m, colX, rowY, colY + m, m), matrixMult(x, y, rowX + m, colX + m, rowY + m, colY + m, m), m, m);
	    }
	    return z;
	}

	private static void matrixSum(int[][] z, int[][] x, int[][] y, int rowZ, int colZ){
	    int n = x.length;
	    for(int i =0; i<n; i++){
	        for(int j=0; j<n; j++)  
	            z[i+rowZ][j+colZ] = x[i][j] + y[i][j];
	    }
	}
	
	//strassen method
	private static int[][] sMultiply(int[][] x, int[][] y){
		int[][] z = new int[x.length][x.length];
		int n = x.length;
		if(n == 1)
			z[0][0] = x[0][0] * y[0][0];
		else{
			int[][] x11 = new int[n/2][n/2];
			int[][] x12 = new int[n/2][n/2];
			int[][] x21 = new int[n/2][n/2];
			int[][] x22 = new int[n/2][n/2];
			int[][] y11 = new int[n/2][n/2];
			int[][] y12 = new int[n/2][n/2];
			int[][] y21 = new int[n/2][n/2];
			int[][] y22 = new int[n/2][n/2];
			
			//split x and y into 4 halves
			split(x, x11, 0, 0);
			split(x, x12, 0, n/2);
			split(x, x21, n/2, 0);
			split(x, x22, n/2, n/2);
			
			split(y, y11, 0, 0);
			split(y, y12, 0, n/2);
			split(y, y21, n/2, 0);
			split(y, y22, n/2, n/2);
			//P
			int[][] p = sMultiply(sAdd(x11, x22), sAdd(y11, y22));
			//Q
			int[][] q = sMultiply(sAdd(x21, x22), y11);
			//R
			int[][] r = sMultiply(x11, sSub(y12, y22));
			//S
			int[][] s = sMultiply(x22, sSub(y21, y11));
			//T
			int[][] t = sMultiply(sAdd(x11, x12), y22);
			//U
			int[][] u = sMultiply(sSub(x21, x11), sAdd(y11, y12));
			//V
			int[][] v = sMultiply(sSub(x12, x22), sAdd(y21, y22));
			//calculate C's
			int[][] z11 = sAdd(sSub(sAdd(p, s), t), v);
			int[][] z12 = sAdd(r, t);
			int[][] z21 = sAdd(q, s);
			int[][] z22 = sAdd(sSub(sAdd(p, r), q), u);
			//combine C's into resulting C matrix
			comb(z11, z, 0, 0);
			comb(z12, z, 0, n/2);
			comb(z21, z, n/2, 0);
			comb(z22, z, n/2, n/2);
		}
		return z;
	}
	
	//matrix addition
	private static int[][] sAdd(int[][] x, int[][] y){
		int n = x.length;
		int[][] z = new int[n][n];
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				z[i][j] = x[i][j] + y[i][j];
			}
		}
		return z;
	}
	
	//matrix subtraction
	private static int[][] sSub(int[][] x, int[][] y){
		int n = x.length;
		int[][] z = new int[n][n];
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				z[i][j] = x[i][j] - y[i][j];
			}
		}
		return z;
	}
	
	//splits matrix into corresponding c11, c12, c21, c22
	private static void split(int[][] x, int[][] y, int n, int m){
		for(int i = 0, i2 = n; i < y.length; i++, i2++){
			for(int j = 0, j2 = m; j < y.length; j++, j2++){
				y[i][j] = x[i2][j2];
			}
		}
	}
	
	//combines back into C matrix
	private static void comb(int[][] x, int[][] y, int n, int m){
		for(int i = 0, i2 = n; i < x.length; i++, i2++){
			for(int j = 0, j2 = m; j < x.length; j++, j2++){
				y[i2][j2] = x[i][j];
			}
		}
	}
	
	//used to print a matrix for testing
	private static void printMatrix(int[][] x) {
	    for (int i = 0; i < x.length; i++) {
	        for (int j = 0; j < x[i].length; j++) {
	            System.out.printf("%12d", x[i][j]);
	        }
	        System.out.println();
	    }
	}
}
