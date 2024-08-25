package Graph;

import java.io.Serializable;
import java.util.ArrayList;

/**
 * A node represents one polygon on a terrain map. It has a number
 * of features: elevation, position and connections to other nodes.
 * DO NOT MODIFY THIS CLASS. A clean version will be used for marking.
 * Originally created for Data Structures, SP2 2017.
 * Reworked for Data Structures Advanced, SP5 2022
 * @author James Baumeister
 * @author Daniel Ablett
 * @version 5.0
 */
public class Node implements Serializable {

	private static final long serialVersionUID = -3085575079765859612L;

	private double elevation;
	private Position position;
	private ArrayList<Edge> edges;

	public Node(double elevation, Position position) {
		this.elevation = elevation;
		this.position = position;
		edges = new ArrayList<Edge>();
	}

	public double getElevation() {
		return elevation;
	}

	public void setElevation(double elevation) {
		this.elevation = elevation;
	}

	public Position getPosition() {
		return position;
	}

	public void setPosition(Position position) {
		this.position = position;
	}

	public ArrayList<Edge> getEdges() {
		return edges;
	}

	public void setEdges(ArrayList<Edge> edges) {
		this.edges = edges;
	}

	@Override
	public String toString() {
		return "Node{" +
				"elevation=" + elevation +
				", position=" + position +
				'}';
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) {
			return true;
		}
		if (o == null || getClass() != o.getClass()) {
			return false;
		}

		Node node = (Node) o;

		if (Double.compare(node.elevation, elevation) != 0) {
			return false;
		}
		if (position != null ? !position.equals(node.position) : node.position != null) {
			return false;
		}
		if (edges == null || ((Node) o).getEdges() == null)
			return false;

		for (Edge e : edges) {
			if (!((Node) o).getEdges().contains(e)) {
				return false;
			}
		}
		return true;
	}

	@Override
	public int hashCode() {
		int result;
		long temp;
		temp = Double.doubleToLongBits(elevation);
		result = (int) (temp ^ (temp >>> 32));
		result = 31 * result + (position != null ? position.hashCode() : 0);
		return result;
	}
}
