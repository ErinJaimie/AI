public class Particle {
	
	public double pBest[];
	public double position[];
	public double oldPosition[];
	public double velocity[];
	public double lambda;
	public int neighbourhood;
	public boolean dead; // kill a particle if it leaves the search space
	
	public Particle(double pos[], double old[], double best[], double vel[], double lam, int n) {
		pBest = best;
		position = pos;
		velocity = vel;
		lambda = lam;
		neighbourhood = n;
		dead = false;
		oldPosition = old;
	}
	
	public void updatePosition(double newPos[]) {
		position = newPos;
	}
	
	public void updateBest(double best[]) {
		pBest = best;
	}
	
	public void updateVelocity(double vel[]) {
		velocity = vel;
	}
	
	public void printPosition() {
		for (int i = 0; i < position.length; i++) {
			System.out.print(position[i] + ", ");
		}
		System.out.println();
	}
	
	public void kill() {
		dead = true;
	}
	
	public void alive() {
		dead = false;
	}
	
	public void setOld() {
		oldPosition = position;
	}
}
