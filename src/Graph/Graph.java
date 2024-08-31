package Graph;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

/**
 * Represents a non-directional graph where each vertex
 * is a Node object. Connections between nodes are based
 * on the cartesian coordinate system.
 * @author James Baumeister on 30/4/17
 */
public class Graph {

	/**
	 * Connects all nodes, building their E, W, S, N, NE, SE, NW, SW edges,
	 * in that order. The nodes should form a 10x10 square grid, and the array
	 * is such that every i'th node % 10 = 9 is a right edge.
	 * See the assignment specification for more information.
	 * @param nodes An array of Node objects to be connected
	 * @return An array of connected Node objects
	 */
	public void connectNodes(Node[] nodes) {
		int size = 10; // size of the grid
		
		for (int i = 0; i < nodes.length; i++) {
			Node current = nodes[i];
			int row = i / size;
			int column = i % size; // if column = 0, no west, north-west or south-west
			
			if(column > 0) { // west nodes
				Node westNode = nodes[i - 1];
				current.getEdges().add(new Edge(current, westNode));
			}
			if(row > 0) { // north nodes
				Node northNode = nodes[i - size];
				current.getEdges().add(new Edge(current, northNode));
			}
			if(column < size - 1) { // east nodes
				Node eastNode = nodes[i + 1];
				current.getEdges().add(new Edge(current, eastNode));
			}
			if(row < size - 1) { // south nodes
				Node southNode = nodes[i + size];
				current.getEdges().add(new Edge(current, southNode));
			}
			if(row > 0 && column > 0) { // north-west nodes
				Node northWestNode = nodes[i - size - 1];
				current.getEdges().add(new Edge(current, northWestNode));
			}
			if(row > 0 && column < size - 1) { // north-east nodes
				Node northEastNode = nodes[i - size + 1];
				current.getEdges().add(new Edge(current, northEastNode));
			}
			if(row < size - 1 && column > 0) { // south-west nodes
				Node southWestNode = nodes[i + size - 1];
				current.getEdges().add(new Edge(current, southWestNode));
			}
			if(row < size - 1 && column < size - 1) { // south-east nodes
				Node southEastNode = nodes[i + size + 1];
				current.getEdges().add(new Edge(current, southEastNode));
			}

		}

	}


	/**
	 * Searches for an edge from the source node to the destination.
	 * @param source The source, or first, node
	 * @param destination The destination, or second, node
	 * @return The edge between the nodes, or null if not found
	 */
	public Edge getEdge(Node source, Node destination) {
		ArrayList<Edge> edges = source.getEdges();
		
		for(Edge edge : edges) {
			if(edge.getToNode().equals(destination)) {
				return edge;
			}
		}

		return null;
	}

	/**
	 * From an array of Node objects, this calculates the total cost of
	 * travelling (i.e. sum of weights) from the first to the last nodes.
	 * @param vertices An array of Node objects representing a path
	 * @return The total cost of travel, or if no edge Double.POSITIVE_INFINITY. 
	 */
	public double calculateTotalWeight(Node[] vertices) {
		
		if(vertices.length == 0) {
			return 0;
		}
		
		double totalWeight = 0;
		
		for(int i = 0; i < vertices.length - 1; i++) {
			Node fromNode = vertices[i];
			Node toNode = vertices[i + 1];
			
			Edge connectingEdge = null;
			
			for(Edge edge : fromNode.getEdges()) {
				if(edge.getToNode().equals(toNode)) {
					connectingEdge = edge;
				}
			}
			
			if(connectingEdge == null) {
				return Double.POSITIVE_INFINITY;
			}
			
			totalWeight += connectingEdge.getWeight();
		}
		
		return totalWeight;
	}


	/**
	 * Performs a breadth-first search of the graph and returns the shortest
	 * path from one node to another.
	 * @param start The node from which to start searching
	 * @param target The target node to which a path is built
	 * @return An array of Node objects representing the path from start to target, in that order
	 */
	public Node[] breadthFirstSearch(Node start, Node target) {
		Queue<Node> queue = new LinkedList<>();
		Map<Node, Node> path = new HashMap<>();
		Set<Node> visited = new HashSet<>();
		
		queue.offer(start);
		visited.add(start);
		path.put(start, null);

		while(!queue.isEmpty()) {
			Node current = queue.poll();
			
			if(current.equals(target)) {
				List<Node> result = new ArrayList<>();
				
				Node step = target;
				while(step != null) {
					result.add(step);
					step = path.get(step);
				}
				
				return result;
			}
			
			for(Edge edge : current.getEdges()) {
				Node neighbour = edge.getToNode();
				if(!visited.contains(neighbour)) {
					queue.offer(neighbour);
					visited.add(neighbour);
					path.put(neighbour, current);
				}
			}
		}
		
		return new Node[0];
	}


	/**
	 * Performs a depth-first search of the graph and returns the first-found
	 * path from one node to another.
	 * @param start The node from which to start searching
	 * @param target The target node to which a path is built
	 * @return An array of Node objects representing the path from start to target, in that order
	 */
	public Node[] depthFirstSearch(Node start, Node target) {
		// TODO

		return null;
	}


	/**
	 * Performs a search of the graph using Dijkstra's algorithm, which takes into
	 * account the edge weight. It should return the least-costly path from one node
	 * to another.
	 * @param start The node from which to start searching
	 * @param target The target node to which a path is built
	 * @return An array of Node objects representing the path from start to target, in that order
	 */
	public Node[] dijkstrasSearch(Node start, Node target) {
		// TODO

		return null;
	}

}
