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

import GUI.GraphLoader;

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
			int row = i / size; // calculates the current row index for the given node
			int column = i % size; // if column = 0, no west, north-west or south-west. Calculates column index for the given node
			
			if(column < size - 1) { // east nodes - basically saying that a node is not in the last column, so would have a node to the right (east)
				Node eastNode = nodes[i + 1]; // this should get the node to the right
				current.getEdges().add(new Edge(current, eastNode)); // edge now exists between current node and east node
			}
			if(column > 0) { // west nodes - opposite of east - as long as the node is not in the first column, it would have a node to the left
				Node westNode = nodes[i - 1]; // this should get the node to the left
				current.getEdges().add(new Edge(current, westNode)); // every time this is called, an edge is created and added to the current node's list of edges
			}
			if(row < size - 1) { // south nodes - if not in the last row (bottom)
				Node southNode = nodes[i + size]; // this should get the node directly below
				current.getEdges().add(new Edge(current, southNode));
			}
			if(row > 0) { // north nodes - if not in the first row (top)
				Node northNode = nodes[i - size]; // this should get the node directly above
				current.getEdges().add(new Edge(current, northNode));
			}
			if(row > 0 && column < size - 1) { // north-east nodes - if not in the first row (same as North node case) or last column (same as East node case)
				Node northEastNode = nodes[i - size + 1]; // gets node directly above and one to the right
				current.getEdges().add(new Edge(current, northEastNode));
			}
			if(row < size - 1 && column < size - 1) { // south-east nodes - if not in the last row and the last column
				Node southEastNode = nodes[i + size + 1]; // gets node directly below and one to the right
				current.getEdges().add(new Edge(current, southEastNode));
			}
			if(row > 0 && column > 0) { // north-west nodes - if not in the first row and first column
				Node northWestNode = nodes[i - size - 1]; // gets node directly above and one to the left
				current.getEdges().add(new Edge(current, northWestNode));
			}
			if(row < size - 1 && column > 0) { // south-west nodes - if not in the last row and first column
				Node southWestNode = nodes[i + size - 1]; // gets node directly below and one to the left
				current.getEdges().add(new Edge(current, southWestNode));
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
			if(edge.getFromNode().equals(source) && edge.getToNode().equals(destination)) {
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
				
				Collections.reverse(result);
				
				return result.toArray(new Node[0]);
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
		if(start == null || target == null) {
			return null;
		}
		
		Set<Node> visited = new HashSet<>();
		List<Node> path = new ArrayList<>();

		if(dfsHelper(start, target, visited, path)) {
			return path.toArray(new Node[0]);
		} else {
			return new Node[0];
		}

	}
	
	public boolean dfsHelper(Node current, Node target, Set<Node> visited, List<Node> path) {
		if(visited.contains(current)) {
			return false;
		}
		
		visited.add(current);
		path.add(current);
		
		if(current.equals(target)) {
			return true;
		}
		
		for(Edge edge : current.getEdges()) {
			Node neighbour = edge.getToNode();
			
			if(!visited.contains(neighbour)) {
				if(dfsHelper(neighbour, target, visited, path)) {
					return true;
				}
			}
		}
		
		path.remove(path.size() - 1);
		return false;
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
		Node[] nodes = GraphLoader.getNodes();
		
		int numV = nodes.length;
		Map<Node, Integer> nodeIndex = new HashMap<>();
		for(int i = 0; i < numV; i++) {
			nodeIndex.put(nodes[i], i);
		}
		
		int startNodeValue = nodeIndex.get(start);
		int targetNodeValue = nodeIndex.get(target);
		Set<Integer> vMinusS = new HashSet<>();
		
		int[] pred = new int[numV];
		double[] dist = new double[numV];
		
		for(int i = 0; i < numV; i++) {
			if(i == startNodeValue) {
				dist[i] = 0;
			}
			else {
				dist[i] = Double.POSITIVE_INFINITY;
			}
			pred[i] = -1; // no predecessors yet
			vMinusS.add(i);
		}
		
		while(!vMinusS.isEmpty()) {
			double minDist = Double.POSITIVE_INFINITY;
			int currentIndex = -1;
			for(int v : vMinusS) {
				if(dist[v] < minDist) {
					minDist = dist[v];
					currentIndex = v;
				}
			}
			
			vMinusS.remove(currentIndex); // Node u is removed which marks it as visited
			Node currentNode = nodes[currentIndex];
			
			for(Edge edge : currentNode.getEdges()) {
				Node neighbour = edge.getToNode();
				int neighbourIndex = nodeIndex.get(neighbour);
				
				if(vMinusS.contains(neighbourIndex)) {
					double weight = edge.getWeight();
					if(dist[currentIndex] + weight < dist[neighbourIndex]) {
						dist[neighbourIndex] = dist[currentIndex] + weight;
						pred[neighbourIndex] = currentIndex;
					}
				}
				
			}
		}
		
		List<Node> path = new ArrayList<>();
		int currentNodeIndex = targetNodeValue;
		while(currentNodeIndex != -1) {
			path.add(nodes[currentNodeIndex]);
			currentNodeIndex = pred[currentNodeIndex];
		}
		Collections.reverse(path);
		
		return path.toArray(new Node[0]);
	}

}
