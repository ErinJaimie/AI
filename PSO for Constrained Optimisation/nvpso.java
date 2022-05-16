import java.util.Arrays;
import java.io.PrintWriter;
import java.io.FileNotFoundException;

public class nvpso {
	public static double w = 0.7;
	public static double c1 = 1.0;
	public static double c2 = 1.0;
	public static int n = 40;
	public static int iterations = 5000;
	public static Particle swarm[] = new Particle[n];
	
	
	// benchmark function: g01
	public static double range[][] = cop.g01Range();
	public static int dim = range[0].length;
	static int bm = 0;
	public static double gBest[] = new double[dim];
	public static double min = cop.g01Min;
	public static String name = "g01";
	
	// benchmark function: g04
//	public static double range[][] = cop.g04Range();
//	public static int dim = range[0].length;
//	static int bm = 6;
//	public static double gBest[] = range[1];
//	public static double min = cop.g04Min;
//	public static String name = "g04";
	
	// benchmark function: g07
//	public static double range[][] = cop.g07Range();
//	public static int dim = range[0].length;
//	static int bm = 3;
//	public static double gBest[] = range[1];
//	public static double min = cop.g07Min;
//	public static String name = "g07";
	
	// benchmark function: g08
//	public static double range[][] = cop.g08Range();
//	public static int dim = range[0].length;
//	static int bm = 10;
//	public static double gBest[] = range[1];
//	public static double min = cop.g08Min;
//	public static String name = "g08";
	
	// benchmark function: g09
//	public static double range[][] = cop.g09Range();
//	public static int dim = range[0].length;
//	static int bm = 4;
//	public static double gBest[] = range[1];
//	public static double min = cop.g09Min;
//	public static String name = "g09";

	public static void main(String[] args) throws FileNotFoundException {
		double sols[][] = new double[iterations][20];
		PrintWriter out = new PrintWriter("NVPSO_"+ name + ".csv");
		PrintWriter rawData = new PrintWriter("rawData_" + name + ".csv");
		out.println("iteration, average, best, worst, std_dev");
		for (int run = 0; run < 20; run++) {				
			// Initialize swarm
			if (bm == 3 || bm == 10 || bm == 4) {
				for (int i = 0; i < dim; i++) {
					gBest[i] = (range[1][i] - range[0][i])/2;
				}
			} else if (bm == 0) {
				gBest = range[0];
			} else {
				gBest = range[1];
			}

			double defaultVel[] = new double[dim];
			Arrays.fill(defaultVel, 0.0);
			for (int i = 0; i < n; i++) {
				double pos[] = generatePosition();
				while (psi(pos) > 1) {
					pos = generatePosition();
				}
				swarm[i] = new Particle(pos, pos, pos, defaultVel, 0, 0);
			}
			
			// Initialize global best position
			double best = eval(gBest);
			for (int i = 0; i < n; i++) {
				if (eval(swarm[i].position) < best && checkDomain(swarm[i].position)) {
					gBest = swarm[i].position;
				}
			}
			
			// run simulator for x iterations
			for (int t = 0; t < iterations; t++) {
				// update particle positions
				for (int j = 0; j < n; j++) {
					updatePosition(swarm[j]);
					double beta = 0.99;
					while (psi(swarm[j].position) > 1 && beta >= 0) {
						beta = beta - 0.01;
						fixPosition(swarm[j], beta);
					}
					if (psi(swarm[j].position) > 1 || !checkDomain(swarm[j].position)) {
						swarm[j].kill();
					}
					// now update personal best
					if (((eval(swarm[j].position) - eval(swarm[j].pBest)) < 0.0) && !swarm[j].dead) {
						swarm[j].pBest = swarm[j].position;
					}
				}
				setGlobalBest();
				sols[t][run] = eval(gBest);
			}
		}
		
		// print data to file
		for (int i = 0; i < iterations; i++) {
			double av = 0.0;
			double min = sols[i][0];
			double max = sols[i][0];
			double std = 0;
			for (int j = 0; j < 20; j++) {
				av += sols[i][j];
				if (sols[i][j] < min) {
					min = sols[i][j];
				} else if (sols[i][j] > max) {
					max = sols[i][j];
				}
				rawData.print(sols[i][j] + ", ");
			}
			rawData.println();
			av = av/20;
			for (int j = 0; j < 20; j++) {
				std += (sols[i][j] - av)*(sols[i][j] - av);
			}
			std = Math.sqrt(std/20);
			out.println(i + ", " + av + ", " + min + ", " + max + ", " + std);
		}
		out.close();
		rawData.close();
	}
	
	// this function must be called PER PARTICLE
	// update velocity -> v(t+1)
	public static void updateVelocity(Particle p) {
		double newVel[] = new double[dim];
		for (int i = 0; i < dim; i++) {
			double r1 = Math.random();
			double r2 = Math.random();
			newVel[i] = w*p.velocity[i] + c1*r1*(p.pBest[i] - p.position[i]) + c2*r2*(gBest[i] - p.position[i]);
		}
		p.updateVelocity(newVel);
	}
	
	// this function must be called PER PARTICLE
	// update position -> x(t+1)
	public static void updatePosition(Particle p) {
		p.setOld();
		double newPos[] = new double[dim];
		updateVelocity(p);
		for (int i = 0; i < dim; i++) {
			newPos[i] = p.position[i] + p.velocity[i];
		}
		double a = shrinkageCoefficient(p, newPos);
		double pos[] = new double[dim];
		for (int i = 0; i < dim; i++) {
			pos[i] = p.position[i] + a*(newPos[i] - p.position[i]);
		}
		p.updatePosition(pos);
		p.alive();
	}
	
	// returns the shrinkage coefficient for dimension d
	public static double shrinkageCoefficient(Particle p, double x_dash[]) {
		double a[] = new double[dim];
		double a_dash = 0.0;
		double m;
		for (int d = 0; d < dim; d++) {
			if (p.position[d] > range[1][d]) {
				m = range[1][d];
			} else if (p.position[d] < range[0][d]){
				m = range[0][d];
			} else {
				m = x_dash[d];
			}
			a[d] = (m - p.position[d])/(x_dash[d] - p.position[d]);
			if (d == 0) {
				a_dash = a[d];
			} else if (a[d] < a_dash) {
				a_dash = a[d];
			}
		}
		return a_dash;
	}
	
	// generates an initial position ~U(x_min, x_max)
	public static double[] generatePosition() {
		double pos[] = new double[dim];
		for (int i = 0; i < dim; i++) {
			pos[i] = Math.random()*(range[1][i]-range[0][i]) + range[0][i];
		}
		return pos;
	}
	
	public static void fixPosition(Particle p, double beta) {
		double pos[] = new double[dim];
		for (int i = 0; i < dim; i++) {
			pos[i] = p.oldPosition[i] + beta*(p.position[i] - p.oldPosition[i]);
		}
		p.alive();
		p.updatePosition(pos);
	}
	
	// determine the global best position
	public static void setGlobalBest() {
		double best = eval(gBest);
		for (int i = 0; i < n; i++) {
			if (((eval(swarm[i].position) - best) < 0.0) && !swarm[i].dead && checkDomain(swarm[i].position)) {
				gBest = swarm[i].position;
				best = eval(swarm[i].position);
			}
		}
	}
	
	public static boolean checkDomain(double x[]) {
		boolean res = true;
		for (int i = 0; i < dim; i++) {
			if ((x[i] > range[1][i]) || (x[i] < range[0][i])) {
				res = false;
			}
		}
		return res;
	}
	
	// evaluate position according to benchmark function
	public static double eval(double x[]) {
		if (bm == 0) {
			return cop.g01(x);
		} else if (bm == 1) {
			return cop.g02(x);
		} else if (bm == 2) {
			return cop.g03(x);
		} else if (bm == 3) {
			return cop.g07(x);
		} else if (bm == 4) {
			return cop.g09(x);
		} else if (bm == 5) {
			return cop.g10(x);
		} else if (bm == 6) {
			return cop.g04(x);
		} else if (bm == 7) {
			return cop.g18(x);
		} else if (bm == 8) {
			return cop.g14(x);
		} else if (bm == 9) {
			return cop.g06(x);
		} else {
			return cop.g08(x);
		}
	}
	
	// evaluates the value of psi 
	public static double psi(double x[]) {
		if (bm == 0) {
			return cop.psiG01(x);
		} else if (bm == 1) {
			return cop.psiG02(x);
		} else if (bm == 2) {
			return cop.psiG03(x);
		} else if (bm == 3) {
			return cop.psiG07(x);
		} else if (bm == 4) {
			return cop.psiG09(x);
		} else if (bm == 5) {
			return cop.psiG10(x);
		} else if (bm == 6) {
			return cop.psiG04(x);
		} else if (bm == 7) {
			return cop.psiG18(x);
		} else if (bm == 8) {
			return cop.psiG14(x);
		} else if (bm == 9) {
			return cop.psiG06(x);
		} else {
			return cop.psiG08(x);
		}
	}
}