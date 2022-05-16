import java.util.Arrays;
import java.io.PrintWriter;
import java.io.FileNotFoundException;

public class assignment1 {
	public static double w = 1.0;
	public static double c1 = 2.0;
	public static double c2 = 2.0;
	public static double K[] = {0.0, 0.1, 0.3, 0.5};
	public static int n = 30;
	public static int iterations = 5000;
	public static int dim = 30;
	public static Particle swarm[] = new Particle[n];
	public static double gBest[] = new double[dim];
	public static boolean left[] = new boolean[n];
	public static double knownMinimum[] = new double[dim];
	
	// benchmark function: Styblinski Tank
//	public static double min = -5.0;
//	public static double max = 5.0;
//	static int bm = 0;
	
	// benchmark function: Schwefel
//	public static double min = -500.0;
//	public static double max = 500.0;
//	static int bm = 1;
	
	// benchmark function: Sphere
//	public static double min = -5.12;
//	public static double max = 5.12;
//	static int bm = 2;
	
	// benchmark function: Powell
//	public static double min = -1.0;
//	public static double max = 1.0;
//	static int bm = 3;
	
	// benchmark function: Brown
	public static double min = -1.0;
	public static double max = 4.0;
	static int bm = 4;


	public static void main(String[] args) throws FileNotFoundException {
		double div[][][] = new double[iterations][20][4];
		double qual[][][] = new double[iterations][20][4];
		double vm[][][] = new double[iterations][20][4];
		double perc[][][] = new double[iterations][20][4];
		
		if (bm == 0) {
			Arrays.fill(knownMinimum, -2.903534);
		} else if (bm == 1) {
			Arrays.fill(knownMinimum, 420.9687);
		} else {
			Arrays.fill(knownMinimum, 0.0);
		}
		
		for (int k = 0; k < K.length; k++) {
			for (int run = 0; run < 20; run++) {
				
				// initialize 'left' boolean
				Arrays.fill(left, false);
				
				// Initialize swarm
				double defaultVel[] = new double[dim];
				Arrays.fill(defaultVel, 0.0);
				for (int i = 0; i < n; i++) {
					double pos[] = generatePosition();
					swarm[i] = new Particle(pos, pos, defaultVel);
				}
				
				// Initialize global best position
				double best = eval(swarm[0].position);
				for (int i = 0; i < n; i++) {
					if (eval(swarm[i].position) < best) {
						gBest = swarm[i].position;
					}
				}
				
				// run simulator for x iterations
				for (int t = 0; t < iterations; t++) {
					Arrays.fill(left, false);
					// update particle positions
					for (int j = 0; j < n; j++) {
						updatePosition(swarm[j], K[k]);
						if (swarm[j].dead) left[j] = true; 
					}
					setGlobalBest();
					div[t][run][k] = diversity(swarm);
					qual[t][run][k] = quality(gBest);
					vm[t][run][k] = velocityMagnitude(swarm);
					perc[t][run][k] = escaped();
				}
				System.out.println("globalBest: " + eval(gBest));
			}
		}
		
		PrintWriter out = new PrintWriter("Results/" + bm + "Diversity_mean.csv");
		PrintWriter sout = new PrintWriter("Results/" + bm + "Diversity_std.csv");
		out.println("t, no clamping, k=0.1, k=0.3, k=0.5");
		sout.println("t, no clamping, k=0.1, k=0.3, k=0.5");
		for (int i = 0; i < iterations; i++) {
			out.print(i + ", ");
			sout.print(i + ", ");
			for (int k = 0; k < 4; k++) {
				double sum = 0.0;
				double sdev = 0.0;
				for (int j = 0; j < 20; j++) {
					sum += div[i][j][k];
				}
				out.print(sum/20.0 + ", ");
				for (int j = 0; j < 20; j++) {
					sdev += (div[i][j][k]-(sum/20.0))*(div[i][j][k]-(sum/20.0));
				}
				sout.print(Math.sqrt(sdev/20.0) + ", ");
			}
			out.println();
			sout.println();
		}
		out.close();
		sout.close();
		
		out = new PrintWriter("Results/" + bm + "Quality_mean.csv");
		sout = new PrintWriter("Results/" + bm + "Quality_std.csv");
		out.println("t, no clamping, k=0.1, k=0.3, k=0.5");
		sout.println("t, no clamping, k=0.1, k=0.3, k=0.5");
		for (int i = 0; i < iterations; i++) {
			out.print(i + ", ");
			sout.print(i + ", ");
			for (int k = 0; k < 4; k++) {
				double sum = 0.0;
				double sdev = 0.0;
				for (int j = 0; j < 20; j++) {
					sum += qual[i][j][k];
				}
				out.print(sum/20.0 + ", ");
				for (int j = 0; j < 20; j++) {
					sdev += (qual[i][j][k]-(sum/20.0))*(qual[i][j][k]-(sum/20.0));
				}
				sout.print(Math.sqrt(sdev/20.0) + ", ");
			}
			out.println();
			sout.println();
		}
		out.close();
		sout.close();
		
		out = new PrintWriter("Results/" + bm + "Escaped_mean.csv");
		sout = new PrintWriter("Results/" + bm + "Escaped_std.csv");
		out.println("t, no clamping, k=0.1, k=0.3, k=0.5");
		sout.println("t, no clamping, k=0.1, k=0.3, k=0.5");
		for (int i = 0; i < iterations; i++) {
			out.print(i + ", ");
			sout.print(i + ", ");
			for (int k = 0; k < 4; k++) {
				double sum = 0.0;
				double sdev = 0.0;
				for (int j = 0; j < 20; j++) {
					sum += perc[i][j][k];
				}
				out.print(sum/20.0 + ", ");
				for (int j = 0; j < 20; j++) {
					sdev += (perc[i][j][k]-(sum/20.0))*(perc[i][j][k]-(sum/20.0));
				}
				sout.print(Math.sqrt(sdev/20.0) + ", ");
			}
			out.println();
			sout.println();
		}
		out.close();
		sout.close();

		out = new PrintWriter("Results/" + bm + "VM_mean.csv");
		sout = new PrintWriter("Results/" + bm + "VM_std.csv");
		out.println("t, no clamping, k=0.1, k=0.3, k=0.5");
		sout.println("t, no clamping, k=0.1, k=0.3, k=0.5");
		for (int i = 0; i < iterations; i++) {
			out.print(i + ", ");
			sout.print(i + ", ");
			for (int k = 0; k < 4; k++) {
				double sum = 0.0;
				double sdev = 0.0;
				for (int j = 0; j < 20; j++) {
					sum += vm[i][j][k];
				}
				out.print(sum/20.0 + ", ");
				for (int j = 0; j < 20; j++) {
					sdev += (vm[i][j][k]-(sum/20.0))*(vm[i][j][k]-(sum/20.0));
				}
				sout.print(Math.sqrt(sdev/20.0) + ", ");
			}
			out.println();
			sout.println();
		}
		out.close();
		sout.close();
	}
	
	// this function must be called PER PARTICLE
	// update velocity -> v(t+1)
	public static void updateVelocity(Particle p, double k) {
		double newVel[] = new double[dim];
		for (int i = 0; i < dim; i++) {
			double r1 = Math.random();
			double r2 = Math.random();
			// need to check that vector math is correct.
			newVel[i] = w*p.velocity[i] + c1*r1*(p.pBest[i] - p.position[i]) + c2*r2*(gBest[i] - p.position[i]);
		}
		
		// if clamping is set to true, check that new velocity does not exceed V_max
		if (k >= 0.1) {
			double Vmax = k*(max - min);
			for (int i = 0; i < dim; i++) {
				if (Math.abs(newVel[i]) > Vmax) {
					if (newVel[i] > Vmax) {
						newVel[i] = Vmax;
					} else {
						newVel[i] = -1*Vmax;
					}
				}
			}
		}
		
		p.updateVelocity(newVel);
	}
	
	// this function must be called PER PARTICLE
	// update position -> x(t+1)
	public static void updatePosition(Particle p, double k) {
		double newPos[] = new double[dim];
		updateVelocity(p, k);
		for (int i = 0; i < dim; i++) {
			newPos[i] = p.position[i] + p.velocity[i];
		}
		p.updatePosition(newPos);
		p.alive();
		// kill particle if position update caused it to leave the search space
		for (int i = 0; i < dim; i++) {
			if (p.position[i] > max || p.position[i] < min) {
				p.kill();
				break;
			}
		}
		if (((eval(p.position) - eval(p.pBest)) < 0.0) && !p.dead) {
			p.pBest = p.position;
		}
	}
	
	// generates an initial position ~U(x_min, x_max)
	public static double[] generatePosition() {
		double pos[] = new double[dim];
		for (int i = 0; i < dim; i++) {
			pos[i] = Math.random()*(max-min) + min;
		}
		return pos;
	}
	
	// determine the global best position
	public static void setGlobalBest() {
		double best = eval(gBest);
		for (int i = 0; i < n; i++) {
			if (((eval(swarm[i].position) - best) < 0.0) && !swarm[i].dead) {
				gBest = swarm[i].position;
				best = eval(swarm[i].position);
			}
		}
	}
	
	// evaluate position according to benchmark function
	public static double eval(double x[]) {
		if (bm == 0) {
			return benchmark.styblinski(dim, x);
		} else if (bm == 1) {
			return benchmark.schwefel(dim, x);
		} else if (bm == 2) {
			return benchmark.sphere(dim, x);
		} else if (bm == 3) {
			return benchmark.powell(dim, x);
		} else {
			return benchmark.brown(dim, x);
		}
	}
	
	// measures the quality of the global best position found by the swarm
	public static double quality(double gBest[]) {
		return eval(gBest);
	}
	
	// measures the diversity of the swarm over time
	// as the average Euclidean distance that particles are from the center of mass of the swarm
	public static double diversity(Particle swarm[]) {
		// get average position (assuming this is the particle centroid)
		double center[] = new double[dim];
		for (int i = 0;  i < n; i++) {
			for (int j = 0; j < dim; j++) {
				center[i] += swarm[i].position[j]/(double) n;
			}
		}
		
		// now find average Euclidean distance 
		double totalDist = 0.0;
		for (int i = 0; i < n; i++) {
			double dist = 0.0;
			for (int j = 0; j < dim; j++) {
				dist += (swarm[i].position[j] - center[j])*(swarm[i].position[j] - center[j]);
			}
			totalDist += Math.sqrt(dist);
		}
		return totalDist/(double) n;
	}
	
	// measures the percentage of particles that leave the search space boundaries
	public static double escaped() {
		int cnt = 0;
		for (int i = 0; i < n; i++) {
			if (left[i]) {
				cnt++;
			}
		}
		return (double) cnt/(double) n;
	}
	
	// measures the average velocity magnitude over time
	// needs to be called at every point in time
	public static double velocityMagnitude(Particle swarm[]) {
		double magnitudeSum = 0.0;
		for (int i = 0; i < n; i++) {
			double mag = 0.0;
			for (int j = 0; j < dim; j++) {
				mag += swarm[i].velocity[j]*swarm[i].velocity[j];
			}
			magnitudeSum += Math.sqrt(mag);
		}
		return magnitudeSum/(double) n;
	}
	

}
