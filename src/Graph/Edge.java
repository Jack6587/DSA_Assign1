package Graph;

import java.io.Serializable;
import java.util.Map;

/**
 * Represents a path, or connection, between two nodes.
 * Originally created for Data Structures, SP2 2017.
 * Reworked for Data Structures Advanced, SP5 2022
 * @author James Baumeister
 * @author Daniel Ablett, Gun Lee
 * @version 5.1
 */
public class Edge implements Serializable {
	private Node fromNode;
	private Node toNode;
	private double weight;

	/**
	 * @param node1
	 * @param node2
	 */
	public Edge(Node fromNode, Node toNode) {
		this.fromNode = fromNode;
		this.toNode = toNode;
		this.weight = calculateWeight();
		
	}

	public Node getFromNode() {
		return fromNode;
	}

	public void setFromNode(Node node) {
		this.fromNode = node;
		this.weight = calculateWeight();
	}

	public Node getToNode() {
		return toNode;
	}

	public void setToNode(Node node) {
		this.toNode = node;
		this.weight = calculateWeight();
	}

	/**
	 * Calculates the weight of travel along this edge,
	 * taking into account each node's elevation.
	 * The weight should be calculated as:
	 * w(e) = d(p1, p2) * (0.01+|e2-e1|)
	 * where e1 is the elevation of the "from" node, and
	 * e2 is the elevation of the "to" node. d is a
	 * function that calculates the Euclidean distance
	 * between two 2D points, p1 and p2, that are the positions of each node.
	 */
	private double calculateWeight() {
		Position fromNodePosition = fromNode.getPosition(); // from node position - start point of the edge
		Position toNodePosition = toNode.getPosition(); // to node position - end point of the edge
		
		double fromNodeElevation = fromNode.getElevation(); // elevation of from node
		double toNodeElevation = toNode.getElevation(); // elevation of to node
		
		double distance = fromNodePosition.distance(toNodePosition); // Euclidean distance between from and to node
		
		double weight = distance * (0.01 + Math.abs(fromNodeElevation - toNodeElevation)); // calculates weight (distance multiplied by absolute difference between the two elevations 
		
		return weight; // returns calculated weight
	}

	public double getWeight() {
		return weight;
	}

	@Override
	public String toString() {
		return "Edge{" +
				"fromNode=" + fromNode +
				", toNode=" + toNode +
				", weight=" + weight +
				'}';
	}
}
