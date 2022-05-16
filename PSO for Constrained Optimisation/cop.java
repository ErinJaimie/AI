public class cop {
	public static double g01Min = -15;
	public static double g02Min = -0.80361910412559;
	public static double g03Min = -1.00050010001000;
	public static double g04Min = -3.066553867178332E004;
	public static double g06Min = -6961.81387558015;
	public static double g07Min = 24.30620906818;
	public static double g08Min = -0.0958250414180359;
	public static double g09Min = 680.630057374402;
	public static double g10Min = 7049.24802052867;
	public static double g14Min = -47.7648884594915;
	public static double g18Min = -0.866025403784439;
	
	public static double g01(double x[]) {
		double sum1 = 0;
		double sum2 = 0;
		double sum3 = 0;
		for (int i = 0; i < 4; i++) {
			sum1 += x[i];
			sum2 += x[i]*x[i];
		}
		for (int i = 4; i < 13; i++) {
			sum3 += x[i];
		}
		return 5*sum1 - 5*sum2 - sum3;
	}
	
	public static double psiG01(double x[]) {
		double psi = 1;
		double g[] = new double[9];
		g[0] = Math.max(0, 2*x[0] + 2*x[1] + x[9] +x[10] - 10);
		g[1] = Math.max(0, 2*x[0] + 2*x[2] + x[9] + x[11] - 10);
		g[2] = Math.max(0, 2*x[1] + 2*x[2] * x[10] + x[11] -10);
		g[3] = Math.max(0, -8*x[0] +x[9]);
		g[4] = Math.max(0, -8*x[1] + x[10]);
		g[5] = Math.max(0, -8*x[2] + x[11]);
		g[6] = Math.max(0, -2*x[3] - x[4] + x[9]);
		g[7] = Math.max(0, -2*x[5] - x[6] + x[10]);
		g[8] = Math.max(0, -2*x[7] - x[8] + x[11]);
		for (int i = 0; i < 9; i++) {
			psi *= Math.exp(g[i]);
		}
		return psi;
	}
	
	public static double[][] g01Range() {
		double range[][] = new double[2][13];
		for (int i = 0; i < 13; i++) {
			if (i < 9) {
				range[0][i] = 0;
				range[1][i] = 1;
			} else if (i < 12) {
				range[0][i] = 0;
				range[1][i] = 100;
			} else {
				range[0][i] = 0;
				range[1][i] = 1;
			}
		}
		return range;
	}
	
	public static double g02(double x[]) {
		double sum1 = 0.0;
		double prod1 = 1;
		double sum2 = 0;
		for (int i = 0; i < 20; i++) {
			sum1 += Math.pow(Math.cos(x[i]), 4);
			prod1 *= Math.pow(Math.cos(x[i]), 2);
			sum2 += i*x[i]*x[i];
		}
		return -1*Math.abs((sum1 - 2*prod1)/Math.sqrt(sum2));
	}
	
	public static double psiG02(double x[]) {
		double psi = 1;
		double g[] = new double[2];
		double prod = 1;
		double sum = 0;
		for (int i = 0; i < 20; i++) {
			prod *= x[i];
			sum += x[i];
		}
		g[0] = Math.max(0, 0.75 - prod);
		g[1] = Math.max(0, sum-7.5*20);
		for (int i = 0; i < 2; i++) {
			psi *= Math.exp(g[i]);
		}
		return psi;
	}
	
	public static double[][] g02Range() {
		double range[][] = new double[2][20];
		for (int i = 0; i < 20; i++) {
			range[0][i] = 0;
			range[1][i] = 10;
		}
		return range;
	}
	
	public static double g03(double x[]) {
		double prod1 = 1;
		for (int i = 0; i < 10; i++) {
			prod1 *= x[i];
		}
		return -1*Math.pow(Math.sqrt(10), 10)*prod1;
	}
	
	public static double psiG03(double x[]) {
		double psi = 1;
		double h[] = new double[2];
		double sum = 0;
		for (int i = 0; i < 10; i++) {
			sum += x[i]*x[i];
		}
		h[0] = Math.max(0, Math.abs(sum - 1 - 1.0E-100));
		h[1] = Math.max(0, Math.abs(-sum + 1 - 1.0E-100));
		for (int i = 0; i < 2; i++) {
			psi *= Math.exp(h[i]);
		}
		return psi;
	}
	
	public static double[][] g03Range() {
		double range[][] = new double[2][10];
		for (int i = 0; i < 10; i++) {
			range[0][i] = 0;
			range[1][i] = 1;
		}
		return range;
	}
	
	public static double g04(double x[]) {
		return 5.3578547*x[2]*x[2] + 0.8356891*x[0]*x[4] + 37.293239*x[0] - 40792.141;
	}
	
	public static double[][] g04Range() {
		double range[][] = new double[2][5];
		for (int i = 0; i < 5; i++) {
			if (i == 0) {
				range[0][i] = 78;
				range[1][i] = 102;
			} else if ( i == 1) {
				range[0][i] = 33;
				range[1][i] = 45;
			} else {
				range[0][i] = 27;
				range[1][i] = 45;
			}
		}
		return range;
	}
	
	public static double psiG04(double x[]) {
		double psi = 1;
		double g[] = new double[6];
		g[0] = Math.max(0, 85.334407 + 0.0056858*x[1]*x[4] + 0.0006262*x[0]*x[3] - 0.0022053*x[2]*x[4] - 92);
		g[1] = Math.max(0, -85.334407 - 0.0056858*x[1]*x[4] - 0.0006262*x[0]*x[3] + 0.0022053*x[2]*x[4]);
		g[2] = Math.max(0, 80.51249 + 0.0071317*x[1]*x[4] + 0.0029955*x[0]*x[1] + 0.0021813*x[2]*x[2] - 110);
		g[3] = Math.max(0, -80.51249 - 0.0071317*x[1]*x[4] - 0.0029955*x[0]*x[1] - 0.0021813*x[2]*x[2] + 90);
		g[4] = Math.max(0, 9.300961 + 0.0047026*x[2]*x[4] + 0.0012547*x[0]*x[2] + 0.0019085*x[2]*x[3] - 25);
		g[5] = Math.max(0, -9.300961 - 0.0047026*x[2]*x[4] - 0.0012547*x[0]*x[2] - 0.0019085*x[2]*x[3] + 20);
		for (int i = 0; i < 6; i++) {
			psi *= Math.exp(g[i]);
		}
		return psi;
	}
	
	public static double g06(double x[]) {
		return (x[0] - 10)*(x[0] - 10)*(x[0] - 10) + (x[1] - 20)*(x[1] - 20)*(x[1] - 20);
	}
	
	public static double[][] g06Range() {
		double range[][] = new double[2][2];
		range[0][0] = 13;
		range[1][0] = 100;
		range[0][1] = 0;
		range[1][1] = 100;
		return range;
	}
	
	public static double psiG06(double x[]) {
		double psi = 1;
		double g[] = new double[2];
		g[0] = Math.max(0,  -1*(x[0] - 5)*(x[0] - 5) - (x[1] - 5)*(x[1] - 5) + 100);
		g[1] = Math.max(0,  (x[0] - 6)*(x[0] - 6) + (x[1] - 5)*(x[1] - 5) - 82.81);
		for (int i = 0; i < 2; i++) {
			psi *= Math.exp(g[i]);
		}
		return psi;
	}
	
	public static double g07(double x[]) {
		return x[0]*x[0] + x[1]*x[1] + x[0]*x[1] - 14*x[0] - 16*x[1] + (x[2] - 10)*(x[2] - 10) + 4*(x[3] - 5)*(x[3] - 5) + 
				(x[4] - 3)*(x[4] - 3) + 2*(x[5] - 1)*(x[5] - 1) + 5*x[6]*x[6] + 7*(x[7] - 11)*(x[7] - 11) + 
				2*(x[8] - 10)*(x[8] - 10) + (x[9] - 7)*(x[9] - 7) + 45;
	}
	
	public static double[][] g07Range() {
		double range[][] = new double[2][10];
		for (int i = 0; i < 10; i++) {
			range[0][i] = -10;
			range[1][i] = 10;
		}
		return range;
	}
	
	public static double psiG07(double x[]) {
		double psi = 1;
		double g[] = new double[8];
		g[0] = Math.max(0,  -105 + 4*x[0] + 5*x[1] - 3*x[6] + 9*x[7]);
		g[1] = Math.max(0,  10*x[0] - 8*x[1] - 17*x[6] + 2*x[7]);
		g[2] = Math.max(0,  -8*x[0] + 2*x[1] + 5*x[8] - 2*x[9] - 12);
		g[3] = Math.max(0, 3*(x[0] - 2)*(x[0] - 2) + 4*(x[1] - 3)*(x[1] - 3) + 2*x[2]*x[2] - 7*x[3] - 120);
		g[4] = Math.max(0, 5*x[0]*x[0] + 8*x[1] + (x[2] - 6)*(x[2] - 6) - 2*x[3] -40);
		g[5] = Math.max(0, x[0]*x[0] + 2*(x[1] - 2)*(x[1] - 2) - 2*x[0]*x[1] + 14*x[4] - 6*x[5]);
		g[6] = Math.max(0,  0.5*(x[0] - 8)*(x[0] - 8) + 2*(x[1] - 4)*(x[1] - 4) + 3*x[4] - x[5] -30);
		g[7] = Math.max(0, -3*x[0] + 6*x[1] + 12*(x[8] - 8)*(x[8] - 8) - 7*x[9]);
		for (int i = 0; i < 8; i++) {
			psi *= Math.exp(g[i]);
		}
		return psi;
	}
	
	public static double g08(double x[]) {
		return -1*(((Math.pow(Math.sin(2*Math.PI*x[0]), 3))*(Math.sin(2*Math.PI*x[1])))/(Math.pow(x[0],  3)*(x[0] + x[1])));
	}
	
	public static double[][] g08Range() {
		double range[][] = new double[2][2];
		for (int i = 0; i < 2; i++) {
			range[0][i] = 0;
			range[1][i] = 10;
		}
		return range;
	}
	
	public static double psiG08(double x[]) {
		double psi = 1;
		double g[] = new double[2];
		g[0] = Math.max(0,  x[0]*x[0] - x[1] + 1);
		g[1] = Math.max(0,  1 - x[0] + (x[1] - 4)*(x[1] - 4));
		for (int i = 0; i < 2; i++) {
			psi *= Math.exp(g[i]);
		}
		return psi;
	}
	
	public static double g09(double x[]) {
		return (x[0] - 10)*(x[0] - 10) + 5*(x[1] - 12)*(x[1] - 12) + x[2]*x[2]*x[2]*x[2] + 3*(x[3] - 11)*(x[3] - 11)
				+ 10*Math.pow(x[4], 6) + 7*x[5]*x[5] + Math.pow(x[6], 4) - 4*x[5]*x[6] - 10*x[5] - 8*x[6];
	}
	
	public static double psiG09(double x[]) {
		double psi = 1;
		double g[] = new double[4];
		g[0] = Math.max(0,  -127 + 2*x[0]*x[0] + 3*Math.pow(x[1], 4) + x[2] + 4*x[3]*x[3] + 5*x[4]);
		g[1] = Math.max(0,  -282 + 7*x[0] + 3*x[1] + 10*x[2]*x[2] + x[3] - x[4]);
		g[2] = Math.max(0,  -196 +23*x[0] +x[1]*x[1] + 6*x[5]*x[5] - 8*x[6]);
		g[3] = Math.max(0, 4*x[0]*x[0] + x[1]*x[1] - 3*x[0]*x[1] + 2*x[2]*x[2] + 5*x[5] - 11*x[6]);
		for (int i = 0; i < 4; i++) {
			psi *= Math.exp(g[i]);
		}
		return psi;
	}
	
	public static double[][] g09Range() {
		double range[][] = new double[2][7];
		for (int i = 0; i < 7; i++) {
			range[0][i] = -1;
			range[1][i] = 10;

		}
		return range;
	}
	
	public static double g10(double x[]) {
		return x[0] + x[1] + x[2];
	}
	
	public static double psiG10(double x[]) {
		double psi = 1;
		double g[] = new double[6];
		g[0] = Math.max(0,  -1 + 0.0025*(x[3]+x[5]));
		g[1] = Math.max(0,  -1 + 0.0025*(x[4]+x[6]-x[3]));
		g[2] = Math.max(0,  -1 + 0.01*(x[7] - x[4]));
		g[3] = Math.max(0, -x[0]*x[5] + 833.33252*x[3] + 100*x[0] - 83333.333);
		g[4] = Math.max(0, -x[1]*x[6] + 1250*x[4] +x[1]*x[3] - 1250*x[3]);
		g[5] = Math.max(0, -1*x[2]*x[7] + 1250000 + x[2]*x[4] - 2500*x[4]);
		for (int i = 0; i < 6; i++) {
			psi *= Math.exp(g[i]);
		}
		return psi;
	}
	
	public static double[][] g10Range() {
		double range[][] = new double[2][8];
		range[0][0] = 100;
		range[1][0] = 10000;
		for (int i = 1; i < 8; i++) {
			if (i < 3) {
				range[0][i] = 1000;
				range[1][i] = 10000;
			} else {
				range[0][i] = 10;
				range[1][i] = 1000;
			}

		}
		return range;
	}
	
	public static double g14(double x[]) {
		double c[] = {-6.089, -17.164, -34.054, -5.914, -24.721, -14.986, -24.1, -10.708, -26.662, -22.179};
		double sum1 = 0;
		double sum2 = 0;
		for (int i = 0; i < 10; i++) {
			sum1 += x[i];
		}
		for (int i = 0; i < 10; i++) {
			sum2 += c[i] + Math.log(x[i]/sum1);
		}
		return sum2;
	}
	
	public static double psiG14(double x[]) {
		double psi = 1;
		double h[] = new double[3];
		h[0] = Math.abs(x[0] + 2*x[1] + 2*x[2] + x[5] + x[9] - 2);
		h[1] = Math.abs(x[3] + 2*x[4] + x[5] + x[6] - 1);
		h[2] = Math.abs(x[2] + x[6] + x[7] + 2*x[8] + x[9] - 1);
		for (int i = 0; i < 3; i++) {
			psi *= Math.exp(Math.max(0, h[i]));
		}
		return psi;
	}
	
	public static double[][] g14Range() {
		double range[][] = new double[2][10];
		for (int i = 0; i < 10; i++) {
			range[0][i] = 0;
			range[1][i] = 10;
		}
		return range;
	}
	
	public static double g18(double x[]) {
		return -0.5*(x[0]*x[3] - x[1]*x[2] + x[2]*x[8] - x[4]*x[8] + x[4]*x[7] - x[5]*x[6]);
	}
	
	public static double psiG18(double x[]) {
		double psi = 1;
		double g[] = new double[13];
		g[0] = Math.max(0,  x[2]*x[2] + x[3]*x[3] - 1);
		g[1] = Math.max(0,  x[8]*x[8] - 1);
		g[2] = Math.max(0,  x[4]*x[4] + x[5]*x[5] - 1);
		g[3] = Math.max(0, x[0]*x[0] + (x[1] - x[8])*(x[1] - x[8]) - 1);
		g[4] = Math.max(0, (x[0] - x[4])*(x[0] - x[4]) + (x[1] - x[5])*(x[1] - x[5]) - 1);
		g[5] = Math.max(0, (x[0] - x[6])*(x[0] - x[6]) + (x[1] - x[7])*(x[1] - x[7]) - 1);
		g[6] = Math.max(0, (x[2] - x[4])*(x[2] - x[4]) + (x[3] - x[5])*(x[3] - x[5]) - 1);
		g[7] = Math.max(0, (x[2] - x[6])*(x[2] - x[6]) + (x[3] - x[7])*(x[3] - x[7]) - 1);
		g[8] = Math.max(0, x[6]*x[6] + (x[7] - x[8])*(x[7] - x[8]) - 1);
		g[9] = Math.max(0, x[1]*x[2] - x[0]*x[3]);
		g[10] = Math.max(0, -1*x[2]*x[8]);
		g[11] = Math.max(0, x[4]*x[8]);
		g[12] = Math.max(0, x[5]*x[6] - x[4]*x[7]);
		for (int i = 0; i < 13; i++) {
			psi *= Math.exp(g[i]);
		}
		return psi;
	}
	
	public static double[][] g18Range() {
		double range[][] = new double[2][9];
		range[0][8] = 0;
		range[1][8] = 20;
		for (int i = 0; i < 8; i++) {
			range[0][i] = -10;
			range[1][i] = 10;
		}
		return range;
	}
}
