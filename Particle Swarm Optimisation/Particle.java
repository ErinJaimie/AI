public class Particle {
	
	public double pBest[];
	public double position[];
	public double velocity[];
	public boolean dead; // kill a particle if it leaves the search space
	
	public Particle(double pos[], double best[], double vel[]) {
		pBest = best;
		position = pos;
		velocity = vel;
		dead = false;
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
}
