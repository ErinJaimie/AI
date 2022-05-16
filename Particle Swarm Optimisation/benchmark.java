import java.util.Arrays;

// https://towardsdatascience.com/optimization-eye-pleasure-78-benchmark-test-functions-for-single-objective-optimization-92e7ed1d1f12
public class benchmark {
	
	// Egg Crate Function
	// min = -5 | max = 5
	// known minimum: 0
	// ERROR: Limited to 2D
	public static double eggCrate(int dim, double x[]) {

		double sum = 0.0;
		double sinSum = 0.0;
		for (int i = 0; i < dim; i++) {
			sum += x[i]*x[i];
			sinSum += Math.sin(x[i])*Math.sin(x[i]);
		}
		return sum + 24*sinSum;

	}
	
	// Simpleton-n Function
	// min = 0 | max = 10
	// known minimum: -10*dim
	public static double simpleton(int dim, double x[]) {
		double sum = 0.0;
		for (int i = 0; i < dim; i++) {
			sum += x[i];
		}
		return -sum;
	}
	
	// Trigonometric 2 Function
	// min =  -500 | max = 500
	// known minimum: 1
	public static double trigonometric(int dim, double x[]) {
		double x_star[] = new double[dim];
		Arrays.fill(x_star, 0.9);
		double sum = 1.0;
		for (int i = 0; i < dim; i++) {
			sum += 8*(Math.sin(7*(x[i]-x_star[i])*(x[i]-x_star[i])))*(Math.sin(7*(x[i]-x_star[i])*(x[i]-x_star[i]))) + 6*Math.sin(14*(x[1]-x_star[1])*(x[1]-x_star[1]))*Math.sin(14*(x[1]-x_star[1])*(x[1]-x_star[1])) + (x[i] - x_star[i])*(x[i] - x_star[i]);
		}
		return sum;
	}

	// Elliptic Function
	// min = -100 | max = 100
	// known minimum: 0
	public static double elliptic(int dim, double x[]) {
		int m = dim -1;
		int i[] = new int[dim];
		for (int j = 0; j < dim; j++) {
			i[j] = j;
		}
		double sum = 0.0;
		for (int j = 0; j < dim; j++) {
			sum += x[j]*x[j]*Math.pow(10, 6*(i[j]/m));
		}
		return sum;
	}
	
	// Generalised Matyas Function 
	// min = -10 | max = 10
	// known minimum: 0
	// ERROR: limited to 2d
	public static double matyas(int dim, double x[]) {
		double sum = 0.0;
		double prod = 1.0;
		for (int i = 0; i < dim; i++) {
			sum += x[i]*x[i];
			prod *= x[i];
		}
		return 0.26*sum + 0.48*prod;
	}
	
	// Styblinski Tank Function
	// min = -5 | max = 5
	// known minimum: -39.16599*dim f(-2.903534,...)
	// multimodal & non-random
	public static double styblinski(int dim, double x[]) {
		double sum = 0.0;
		for (int i = 0; i < dim; i++) {
			sum += (Math.pow(x[i], 4)-16*x[i]*x[i]+5*x[i]);
		}
		return 0.5*sum;
	}
	
	// Schwefel Function
	// min = -500 | max = 500
	// known minimum: 0 f(420.9687,...)
	// multimodal & non-random
	public static double schwefel(int dim, double x[]) {
		double sum = 0.0;
		for (int i = 0; i < dim; i++) {
			sum += (x[i]*Math.sin(Math.sqrt(Math.abs(x[i]))));
		}
		return 418.9829*dim - sum;
	}
	
	// Sphere Function
	// min = -5.12 | max = 5.12
	// known minimum: 0 f(0,...)
	// non-multimodal & non-random
	public static double sphere(int dim, double x[]) {
		double sum = 0.0;
		for (int i = 0; i < dim; i++) {
			sum += (x[i]*x[i]);
		}
		return sum;
	}
	
	// Powell Function
	// min = -1 | max = 1
	// known minimum: 0 f(0,...)
	// non-multimodal & non-random
	public static double powell(int dim, double x[]) {
		double sum = 0.0;
		for (int i = 0; i < dim; i++) {
			sum += (Math.pow(Math.abs(x[i]), i+1));
		}
		return sum;
	}
	
	// Rastrigin Function
	// min = -5.12 | max = 5.12
	// known minimum: 0 f(0,...)
	// multimodal & non-random
	public static double rastrigin(int dim, double x[]) {
		double sum = 0.0;
		for (int i = 0; i < dim; i++) {
			sum += (x[i]*x[i] -10*Math.cos(2*Math.PI*x[i]));
		}
		return sum;
	}
	
	// Brown 
	// min = -1 | max = 4
	// known minimum: 0 f(0, ..., 0)
	// non-multimodal, non-random
	public static double brown(int dim, double x[]) {
		double sum = 0.0;
		for (int i = 0; i < dim - 1; i++) {
			sum += Math.pow(Math.pow(x[i], 2), Math.pow(x[i+1], 2)+1) + Math.pow(Math.pow(x[i+1], 2), Math.pow(x[i], 2)+1);
		}
		return sum;
	}

}
