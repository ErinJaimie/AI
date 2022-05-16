// https://pymoo.org/problems/multi/zdt.html
// https://esa.github.io/pagmo2/docs/cpp/problems/zdt.html
	
public class mopBenchmark {
	
	// zdt1: box constrained n-dimensional multi-objective problem
	// n > 1, m = 2, xE[0, 1], dim = 30
	public static double zdt1(int obj, double x[]) {
		if (obj == 0) {
			return x[0];
		} else {
			double sum = 0;
			for (int i = 1; i < x.length; i++) {
				sum += x[i];
			}
			double g = 1 + (9/(x.length-1))*sum;
			return g*(1 - Math.sqrt(x[0]/g));
		}
	}
	
	public static double[][] r() {
		double r[][] = new double[2][30];
		for (int i = 0; i < 30; i++) {
			r[0][i] = 0;
			r[1][i] = 1;
		}
		return r;
	}
	
	// zdt2: box constrained n-dimensional multi-objective problem
	// n > 1, m = 2, xE[0, 1], dim = 30
	public static double zdt2(int obj, double x[]) {
		if (obj == 0) {
			return x[0];
		} else {
			double sum = 0;
			for (int i = 1; i < x.length; i++) {
				sum += x[i];
			}
			double g = 1 + (9/(x.length-1))*sum;
			return g*(1 - (x[0]/g)*(x[0]/g));
		}
	}
	
	// zdt3: box constrained n-dimensional multi-objective problem
	// n > 1, m = 2, xE[0, 1], dim = 30
	public static double zdt3(int obj, double x[]) {
		if (obj == 0) {
			return x[0];
		} else {
			double sum = 0;
			for (int i = 1; i < x.length; i++) {
				sum += x[i];
			}
			double g = 1 + (9/(x.length-1))*sum;
			return g*(1 - Math.sqrt(x[0]/g) - (x[0]/g)*(Math.sin(10*Math.PI*x[0])));
		}
	}
	
	// zdt4: box constrained n-dimensional multi-objective problem
	// n > 1, m = 2, xE[0, 1], dim = 30
	public static double zdt4(int obj, double x[]) {
		if (obj == 0) {
			return x[0];
		} else {
			double sum = 0;
			for (int i = 1; i < x.length; i++) {
				sum += (x[i]*x[i] - 10*Math.cos(4*Math.PI*x[i]));
			}
			double g = 91 + sum;
			return g*(1 - Math.sqrt(x[0]/g))*x[0];
		}
	}
	
	public static double[][] r2() {
		double r[][] = new double[2][10];
		r[0][0] = 0;
		r[1][0] = 1;
		for (int i = 1; i < 10; i++) {
			r[0][i] = -5;
			r[1][i] = 5;
		}
		return r;
	}
	
	// zdt6: box constrained n-dimensional multi-objective problem
	// n > 1, m = 2, xE[0, 1], dim = 10
	public static double zdt6(int obj, double x[]) {
		if (obj == 0) {
			return 1 - Math.exp(-4*x[0])*Math.pow(Math.sin(6*Math.PI*x[0]), 6);
		} else {
			double sum = 0;
			for (int i = 1; i < x.length; i++) {
				sum += x[i]/9;
			}
			double g = 1 + 9*(Math.pow(sum, 0.25));
			return g*(1 - (zdt6(0, x)/g)*(zdt6(0, x)/g));
		}
	}
	
	public static double[][] r3() {
		double r[][] = new double[2][10];
		for (int i = 0; i < 10; i++) {
			r[0][i] = 0;
			r[1][i] = 1;
		}
		return r;
	}

}
