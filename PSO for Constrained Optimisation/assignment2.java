import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Arrays;

public class assignment2 {
	
	public static int n = 20; //swarm size -> are all swarms the same size?
	public static int iterations = 5000; // number of iterations used as stopping condition
	public static int tournament_size = 3;
	public static int archive_size = 300;
	public static int runs = 20;

	// control parameters
	public static double w; // inertia weight
	public static double c1; // acceleration coefficient
	public static double c2; // acceleration coefficient
	public static double c3; // acceleration coefficient
	public static double l; // lambda value: 0.2301950353739577
	
	public static int dim;
	public static int m;
	public static int bm;
	public static double[][] range;
	
	public static Particle swarm[][];
	public static double nBest[][];
	public static double archive[][];
	public static int aCnt = 0;
	
	public static void main(String[] args) throws FileNotFoundException {
		getControlParams();
		PrintWriter out, arch, var;
		for (int scenario = 0; scenario < 8; scenario++) { // total of 11 scenarios
			if (scenario == 0) {						
				dim = 13; //dimensions
				m = 2; // number of objectives -> there is a swarm for each objective
				bm = 0;
				range = copBenchmark.g01Range();
				out = new PrintWriter("MGPSO_g01.csv");
				arch = new PrintWriter("g01_archive.csv");
				var = new PrintWriter("g01_stddev.csv");
			} else if (scenario == 1) {				
				dim = 5; //dimensions
				m = 2; // number of objectives -> there is a swarm for each objective
				bm = 6;
				range = copBenchmark.g04Range();
				out = new PrintWriter("MGPSO_g04.csv");
				arch = new PrintWriter("g04_archive.csv");
				var = new PrintWriter("g04_stddev.csv");
			} else if (scenario == 2) {				
				dim = 10; //dimensions
				m = 2; // number of objectives -> there is a swarm for each objective
				bm = 3;
				range = copBenchmark.g07Range();
				out = new PrintWriter("MGPSO_g07.csv");
				arch = new PrintWriter("g07_archive.csv");
				var = new PrintWriter("g07_stddev.csv");
			} else if (scenario == 3) {				
				dim = 2; //dimensions
				m = 2; // number of objectives -> there is a swarm for each objective
				bm = 15;
				range = copBenchmark.g08Range();
				out = new PrintWriter("MGPSO_g08.csv");
				arch = new PrintWriter("g08_archive.csv");
				var = new PrintWriter("g08_stddev.csv");
			} else if (scenario == 4){	
				dim = 7; //dimensions
				m = 2; // number of objectives -> there is a swarm for each objective
				bm = 4;
				range = copBenchmark.g09Range();
				out = new PrintWriter("MGPSO_g09.csv");
				arch = new PrintWriter("g09_archive.csv");
				var = new PrintWriter("g09_stddev.csv");
			} else if (scenario == 5) {
				dim = 30; //dimensions
				m = 2; // number of objectives -> there is a swarm for each objective
				bm = 8;
				range = mopBenchmark.r();
				out = new PrintWriter("zdt1.csv");
				arch = new PrintWriter("zdt1_archive.csv");
				var = new PrintWriter("zdt1_stddev.csv");
			} else if (scenario == 6) {
				dim = 30; //dimensions
				m = 2; // number of objectives -> there is a swarm for each objective
				bm = 9;
				range = mopBenchmark.r();
				out = new PrintWriter("zdt2.csv");
				arch = new PrintWriter("zdt2_archive.csv");
				var = new PrintWriter("zdt2_stddev.csv");
			} else {
				dim = 30; //dimensions
				m = 2; // number of objectives -> there is a swarm for each objective
				bm = 10;
				range = mopBenchmark.r();
				out = new PrintWriter("zdt3.csv");
				arch = new PrintWriter("zdt3_archive.csv");
				var = new PrintWriter("zdt3_stddev.csv");
			}
			out.println("f1, f2");
			
			double sol[][] = new double[runs][2];
			double sol1[][] = new double[300][runs*2];
			double sol2[][] = new double[5000][runs];
			// run the simulation 20 times
			for (int r =  0; r < runs; r++) {
				getControlParams();
				swarm = new Particle[m][n]; // m swarms each containing n particles
				nBest = new double[m][dim]; // the best known position in each swarm
				archive = new double[archive_size][dim]; // archive contains n solutions
				System.out.println(bm + " | run: " + r);
				aCnt = 0;
				// initialize swarm
				double default_vel[] = new double[dim];
				for (int i = 0; i < m; i++) {
					for (int j = 0; j < n; j++) {
						double pos[] = generatePosition();
						Arrays.fill(default_vel, 0.0);
						swarm[i][j] = new Particle(pos, pos, pos, default_vel, l, i);
						if (j == 0) {
							nBest[i] = range[1];
						}
						updateNbest(i, pos);
					}
				}
				
				// each experiment runs for 5000 iterations
				double curr_best = Double.MAX_VALUE;
				for (int t = 0; t < iterations; t++) {
					// update personal best of particles
					for (int i = 0; i < m; i++) {
						for (int j = 0; j < n; j++) {
							if (!swarm[i][j].dead && eval(swarm[i][j].position, i) < eval(swarm[i][j].pBest, i)) {
								swarm[i][j].updateBest(swarm[i][j].position);
							}
							
							// update nBest
							updateNbest(i, swarm[i][j].pBest);
							
							// update archive with the solution
							if (nd(swarm[i][j].position) && aCnt < archive_size && !swarm[i][j].dead) {
								archive[aCnt] = swarm[i][j].position;
								aCnt++;
								removeDom(swarm[i][j].position);
								if (eval(swarm[i][j].position, 1) == 0) {
									curr_best = eval(swarm[i][j].position, 0);
								}
							} else if (nd(swarm[i][j].position) && !swarm[i][j].dead) {
								deleteArchive();
								archive[aCnt] = swarm[i][j].position;
								aCnt++;
								removeDom(swarm[i][j].position);
								if (eval(swarm[i][j].position, 1) == 0) {
									curr_best = eval(swarm[i][j].position, 0);
								}
							}

						}
					}
					
					for (int i = 0; i < m; i++) {
						for (int j = 0; j < n; j++) {
							// select a solution a from the archive using tournament selection
							double a[] = getGuide();
							// update velocity
							updateVelocity(swarm[i][j], a);
							// update position
							updatePosition(swarm[i][j]);
						}
					}
					sol2[t][r] = curr_best;
				}
				
				// populate arrays for experimental data storage
				double best[] = archive[0];
				double f2 = eval(archive[0], 1);
				for (int j = 0; j < aCnt; j++) {
					if (eval(archive[j], 1) < f2) {
						f2 = eval(archive[j], 1);
						best = archive[j];
					}
				}
				for (int j = 0; j < aCnt; j++) {
					if (eval(archive[j], 1) == f2 && eval(archive[j], 0) <= eval(best, 0)) {
						best = archive[j];
					}
				}
				sol[r][0] = eval(best, 0);
				sol[r][1] = eval(best, 1);
				for (int j = 0; j < aCnt; j++) {
					sol1[j][r*2] = eval(archive[j], 0);
					sol1[j][r*2 + 1] = eval(archive[j], 1);
				}
				
			}
			
			// print data to file
			int feasible = 0;
			double av = 0;
			double best = Double.MAX_VALUE;
			double worst = Double.MIN_VALUE;
			double std = 0;
			for (int i = 0; i < runs; i++) {
				out.println(sol[i][0] + ", " + sol[i][1]);
				if (!(sol[i][1] > 0)) {
					feasible++;
					av += sol[i][0];
					if (sol[i][0] < best || best == Double.MAX_VALUE) {
						best = sol[i][0];
					} else if (sol[i][0] > worst || worst == Double.MIN_VALUE) {
						worst = sol[i][0];
					}
				}
			}
			av = av/feasible;
			for (int i = 0; i < runs; i++) {
				if (!(sol[i][1] > 0)) {
					std += (sol[i][0] - av)*(sol[i][0] - av);
				}
			}
			std = Math.sqrt(std/feasible);
			out.println("average: " + av);
			out.println("best: " + best);
			out.println("worst: " + worst);
			out.println("std_dev: " + std);
			out.println("# feasible solutions: " + feasible);
			out.close();
			for (int i = 0; i < 300; i++) {
				for (int j = 0; j < runs*2; j++) {
					arch.print(sol1[i][j] + ", ");
				}
				arch.println();
			}
			arch.close();
			
			var.println("iteration, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17, r18, r19, r20, av, std_dev");
			for (int i = 0; i < 5000; i++) {
				double sum = 0;
				double dev = 0;
				var.print(i + ", ");
				for (int j = 0; j < runs; j++) {
					sum += sol2[i][j];
					var.print(sol2[i][j] + ", ");
				}
				for (int j = 0; j < runs; j++) {
					dev += (sol2[i][j] - sum/runs)*(sol2[i][j] - sum/runs);
				}
				var.println(sum/runs + ", " + Math.sqrt(dev/runs));
			}
			var.close();
		}
	}
	
	// generates an initial position ~U(x_min, x_max)
	public static double[] generatePosition() {
		double pos[] = new double[dim];
		for (int i = 0; i < dim; i++) {
			pos[i] = Math.random()*(range[1][i]-range[0][i]) + range[0][i];
		}
		return pos;
	}
	
	// update the neighbourhood best position
	public static void updateNbest(int neighbourhood, double pos[]) {
		if (eval(pos, neighbourhood) < eval(nBest[neighbourhood], neighbourhood) && checkDomain(pos)) {
			nBest[neighbourhood] = pos;
		}
	}
	
	// updates the velocity of a particle p
	public static void updateVelocity(Particle p, double a[]) {
		double newVel[] = new double[dim];
		for (int j = 0; j < dim; j++) {
			double r1 = Math.random();
			double r2 = Math.random();
			double r3 = Math.random();
			newVel[j] = w*p.velocity[j] + c1*r1*(p.pBest[j] - p.position[j]) + p.lambda*c2*r2*(nBest[p.neighbourhood][j] - p.position[j]) + (1-p.lambda)*c3*r3*(a[j] - p.position[j]);
		}
		p.updateVelocity(newVel);
	}
	
	// updates the position of a particle p
	public static void updatePosition(Particle p) {
		p.alive();
		double newPos[] = new double[dim];
		for (int i = 0; i < dim; i++) {
			newPos[i] = p.position[i] + p.velocity[i];
		}
		p.updatePosition(newPos);
		for (int i = 0; i < dim; i++) {
			if (newPos[i] > range[1][i] || newPos[i] < range[0][i]) {
				p.kill();
				break;
			}
		}
	}
	
	// evaluates position against a single objective
	public static double eval(double x[], int obj) {
		if (bm == 0) {
			return copBenchmark.g01(obj, x);
		} else if (bm == 1) {
			return copBenchmark.g02(obj, x);
		} else if (bm == 2) {
			return copBenchmark.g03(obj, x);
		} else if (bm == 3) {
			return copBenchmark.g07(obj, x);
		} else if (bm == 4) {
			return copBenchmark.g09(obj, x);
		} else if (bm == 5) {
			return copBenchmark.g10(obj, x);
		} else if (bm == 6){
			return copBenchmark.g04(obj, x);
		} else if (bm == 7){
			return copBenchmark.g18(obj, x);
		} else if (bm == 8) {
			return mopBenchmark.zdt1(obj, x);
		} else if (bm == 9) {
			return mopBenchmark.zdt2(obj, x);
		} else if (bm == 10) {
			return mopBenchmark.zdt3(obj, x);
		} else if (bm == 11) {
			return mopBenchmark.zdt4(obj,  x);
		} else if (bm == 12) {
			return mopBenchmark.zdt6(obj, x);
		} else if (bm == 13) {
			return copBenchmark.g14(obj, x);
		} else if (bm == 14) {
			return copBenchmark.g06(obj,  x);
		} else {
			return copBenchmark.g08(obj, x);
		}
	}
	
	// returns an archive guide.
	// randomly selects <tournament size> solutions from the archive and returns the solution with the largest crowding distance
	public static double[] getGuide() {
		int pool[] = new int[tournament_size];
		double cDist[] = cDist();
		String unique = "";
		for (int i = 0; i < tournament_size; i++) {
			int r = (int) (Math.random()*(aCnt));
			int counter = 0;
			while (unique.contains(Integer.toString(r) + ",") && counter < 10) {
				r = (int) (Math.random()*(aCnt));
				counter++;
			} 
			pool[i] = r;
			unique += Integer.toString(r) + ",";
		}
		double max_dist = cDist[pool[0]];
		int pos = 0;
		for (int i = 1; i < tournament_size; i++) {
			if (cDist[pool[i]] < max_dist) {
				max_dist = cDist[pool[0]];
				pos = i;
			}
		}
		return archive[pos];
	}
	
	// a solution is non-dominated if it improves on some objectives, if there is no entry in the archive that is better than the solution for all objective functions
	// i.e. a solution that is better than some neighbourhoood bests and worse than others
	// check if there is a solution in the archive that is better than solution for all objectives
	public static boolean nd(double solution[]) {
		for (int i = 0; i < aCnt; i++) {
			if (eval(archive[i], 0) == eval(solution, 0) && eval(archive[i], 1) == eval(solution, 1)) {
				return false;
			}
			int cnt = 0;
			for (int j = 0; j < m; j++) {
				if (eval(archive[i], j) <= eval(solution, j)) {
					cnt++;
				}
			}
			if (cnt == m) {
				return false;
			}
		}
		return true;
	}
	
	// removes any dominated solutions from the archive
	public static void removeDom(double[] solution) {
		for (int i = 0; i < aCnt; i++) {
			int cnt = 0;
			double tmp[] = archive[i];
			for (int j = 0; j < m; j++) {
				if (eval(solution, j) < eval(archive[i], j)) {
					cnt++;
				}
			}
			if (cnt == m) {
				for (int k = i; k < aCnt - 1; k++) {
					archive[k] = archive[k + 1];
				}
				aCnt--;
				i--;
				for (int g = 0; g < aCnt; g++) {
					if (archive[g] == tmp) {
						System.out.println("remove failed");
					}
				}
			}
		}
	}
	
	// euclidean distance between 2 points
	public static double dist(double pos1[], double pos2[]) {
		double sum = 0.0;
		for (int i = 0; i < pos1.length; i++) {
			sum += (pos1[i] - pos2[i])*(pos1[i] - pos2[i]);
		}
		return Math.sqrt(sum);
	}
	
	// removes the most crowded solution from the archive
	public static void deleteArchive() {
		// get crowding distances of solutions in the archive
		double cDist[] = cDist();
		// remove the most crowded solution from the archive
		double max = cDist[0];
		int pos = 0;
		for (int i = 1; i < aCnt; i++) {
			if (cDist[i] > max) {
				max = cDist[i];
				pos = i;
			}
		}
		for (int i = pos; i < aCnt - 1; i++) {
			archive[i] = archive[i+1];
		}
		aCnt--;
	}
	
	// returns the crowding distance of each point in archive
	// http://gpbib.cs.ucl.ac.uk/gecco2005/docs/p257.pdf
	public static double[] cDist() {
		double cDist[] = new double[aCnt];
		int order[] = new int[aCnt];
		for (int i = 0; i < m; i++) {
			double max_dist = Math.sqrt(range[1][i]*range[1][i] - range[0][i]*range[0][i]);
			double sorted[][] = sort(i, order);
			cDist[order[0]] += max_dist;
			cDist[order[aCnt-1]] += max_dist;
			for (int j = 1; j < aCnt-1; j++) {
				cDist[order[j]] += (dist(sorted[j-1], sorted[j]) + dist(sorted[j], sorted[j+1]))/2;
			}
		}
		return cDist;
	}
	
	// sorts solutions in the archive according to an objective function value
	public static double[][] sort(int obj, int order[]) {
		double sol[][] = new double[aCnt][dim];
		double eval[] = new double[aCnt];
		for (int i = 0; i < aCnt; i++) {
			eval[i] = eval(archive[i], obj);
		}
		Arrays.sort(eval);
		for (int i = 0; i < aCnt; i++) {
			for (int j = 0; j < aCnt; j++) {
				if (eval[i] == eval(archive[j], obj)) {
					sol[i] = archive[j];
					order[i] = j;
					break;
				}
			}
		}
		return sol;
	}
	
	public static void getControlParams() {
		while (!checkStability()) {
			c1 = Math.random() + 1;
			c2 = Math.random() + 1;
			c3 = Math.random() + 1;
			w = Math.random();
			l = Math.random();
		}
		System.out.println("c1: " + c1 + " | c2: " + c2 + " | c3: " + c3 + " | w: " + w + " | l: " + l);
	}
	
	public static boolean checkStability() {
		if (c1 + l*c2 + (1-l)*c3 > 0 && l > 0) {
			double x = (4*(1-w*w))/(1-w + ((c1*c1+l*l*c2*c2+(1-l)*(1-l)*c3*c3)*(1+w))/(3*(c1+l*c2+(1-l)*c3)*(c1+l*c2+(1-l)*c3)));
			if (c1 + l*c2 + (1-l)*c3 < x) {
				return true;
			}
		}
		return false;
	}
	
	
	// ensures that a point is within the domain of the search space
	public static boolean checkDomain(double x[]) {
		boolean res = true;
		for (int i = 0; i < dim; i++) {
			if ((x[i] > range[1][i]) || (x[i] < range[0][i])) {
				res = false;
			}
		}
		return res;
	}

}
