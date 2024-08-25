package Graph;

import java.io.Serializable;

/**
 * A 2D vector class for various mathematical operations.
 * Originally created for Data Structures, SP2 2017.
 * Reworked for Data Structures Advanced, SP5 2022
 * @author James Baumeister
 * @author Daniel Ablett
 * @version 1.1
 */
public class Position implements Serializable {

	public double x;
	public double y;

	public Position(double x, double y) {
		this.x = x;
		this.y = y;
	}

	/**
	 * Calculates the Euclidean distance between two 2D points.
	 * The formula is:
	 * d(p, q) = sqrt ( (qx - px)^2 + (qy - py)^2 )
	 * @param v2 The second point
	 * @return The distance between the points
	 */
	public double distance(Position v2) {
		// TODO
		
		return -1.0;
	}

	@Override
	public String toString() {
		return "Vector{" +
				"x = " + x +
				", y = " + y +
				'}';
	}
}
