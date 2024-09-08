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
		ArrayList<Edge> edges = source.getEdges(); // obtains the list of edges from the given source node
		
		for(Edge edge : edges) { // iterates over each edge in the obtained list
			if(edge.getFromNode().equals(source) && edge.getToNode().equals(destination)) { 
				// if both the starting point (fromNode) and destination point (toNode) of the edge match the given source and destination nodes, return the edge
				return edge;
			}
		}

		return null; // if no edge is found
	}

	/**
	 * From an array of Node objects, this calculates the total cost of
	 * travelling (i.e. sum of weights) from the first to the last nodes.
	 * @param vertices An array of Node objects representing a path
	 * @return The total cost of travel, or if no edge Double.POSITIVE_INFINITY. 
	 */
	public double calculateTotalWeight(Node[] vertices) {
		
		if(vertices.length == 0) { // if there are no nodes in the given array, the weight is 0
			return 0;
		}
		
		double totalWeight = 0; // initialises a totalWeight variable that will be added to and returned
		
		for(int i = 0; i < vertices.length - 1; i++) { 
			// iterates over the nodes array -1 to find edges. The reason I don't loop over the last one is that I am
			// looking for an edge between the current node and the next node - making it kind of pointless to check the final node if there is nothing after it
			Node fromNode = vertices[i]; // gets current node - so the from node of the edge
			Node toNode = vertices[i + 1]; // gets the next node - so the destination node of the edge
			
			Edge connectingEdge = null; // currently holds the edge between from and to node
			
			for(Edge edge : fromNode.getEdges()) { // iterates over all edges of the from node to find a connection to the to node
				if(edge.getToNode().equals(toNode)) { // if the edge points to the next node
					connectingEdge = edge; // store this edge
				}
			}
			
			if(connectingEdge == null) {
				// no valid path found
				return Double.POSITIVE_INFINITY;
			}
			
			totalWeight += connectingEdge.getWeight(); // add weight of the edge to the total weight
		}
		
		return totalWeight; // return the total calculated weight
	}


	/**
	 * Performs a breadth-first search of the graph and returns the shortest
	 * path from one node to another.
	 * @param start The node from which to start searching
	 * @param target The target node to which a path is built
	 * @return An array of Node objects representing the path from start to target, in that order
	 */
	public Node[] breadthFirstSearch(Node start, Node target) {
		Queue<Node> queue = new LinkedList<>(); // creates a queue to manage nodes needing to be explored
		Map<Node, Node> path = new HashMap<>(); // a map that monitors each node's predecessor from the start to determine the path
		Set<Node> visited = new HashSet<>(); // a set that monitors all nodes that have already been visited - this is so we don't go back and visit nodes that have already been checked
		
		// these add the start node to the queue, and mark it as visited
		queue.offer(start);
		visited.add(start);
		
		path.put(start, null); // starts the path map - the first node, start, has no predecessors, so that value is null

		while(!queue.isEmpty()) { // queue continues until no more nodes to check
			Node current = queue.poll(); // remove the node at the front of the queue to check
			
			if(current.equals(target)) {
				List<Node> result = new ArrayList<>(); // path is reconstructed (from start to target) if target is found
				
				Node step = target; // we start from the target node and work backwards
				while(step != null) {
					result.add(step); // adds step to the result list to be returned
					step = path.get(step); // gets the predecessor node (the value component of the map structure)
				}
				
				Collections.reverse(result); // reverses the list so that the path is ordered from start to target (instead of the other way around, the way it was constructed)
				
				return result.toArray(new Node[0]); // converts result list to an array to return
			}
			
			for(Edge edge : current.getEdges()) { // iterate over all nodes connected to the current node
				Node neighbour = edge.getToNode(); // gets a neighbouring node
				if(!visited.contains(neighbour)) { // if this neighbour has not been visited yet
					queue.offer(neighbour); // add it to the queue 
					visited.add(neighbour); // mark as visited 
					path.put(neighbour, current); // record our current node as the predecessor of this neighbour
				}
			}
		}
		
		return new Node[0]; // returns an empty array if no path is determined
	}


	/**
	 * Performs a depth-first search of the graph and returns the first-found
	 * path from one node to another.
	 * @param start The node from which to start searching
	 * @param target The target node to which a path is built
	 * @return An array of Node objects representing the path from start to target, in that order
	 */
	public Node[] depthFirstSearch(Node start, Node target) {
		if(start == null || target == null) { // if start or target is null, do nothing
			return null;
		}
		
		Set<Node> visited = new HashSet<>(); // set that manages visited nodes to ensure they are not revisited. A set also ensures duplicate entries do not exist
		List<Node> path = new ArrayList<>(); // list that stores the path from start to target

		if(dfsHelper(start, target, visited, path)) { // this calls the helper function to perform DFS
			return path.toArray(new Node[0]); // if path is found, make the list an array and return
		} else {
			return new Node[0]; // otherwise, return an empty array
		}

	}
	
	public boolean dfsHelper(Node current, Node target, Set<Node> visited, List<Node> path) {
		if(visited.contains(current)) { // if the current node has already been visited, exit
			return false;
		}
		
		// add the current node to visited, and add it to the path
		visited.add(current);
		path.add(current);
		
		if(current.equals(target)) { // this indicates a path has been found if current = target
			return true;
		}
		
		for(Edge edge : current.getEdges()) { // iterate over edges of current node
			Node neighbour = edge.getToNode(); // get a neighbouring node
			
			if(!visited.contains(neighbour)) { // if it has not been visited, perform DFS on the node
				if(dfsHelper(neighbour, target, visited, path)) { // recursion - we jump into dfsHelper again and operate on the neighbour, until we reach the target
					return true;
				}
			}
		}
		
		path.remove(path.size() - 1); // this is used for backtracking
		// we remove the most recently added node if we discover the current path does not lead to the target node.
		// We then jump out of this method and back to the previous for loop iterating over edges
		return false; // target node was not found in this path
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
		Node[] nodes = GraphLoader.getNodes(); // loads all the nodes from the graph using GraphLoader class
		
		int numV = nodes.length; // number of nodes in the graph
		Map<Node, Integer> nodeIndex = new HashMap<>(); // used to keep track of each node and its corresponding value, so that we can search and iterate in a simpler way
		for(int i = 0; i < numV; i++) {
			nodeIndex.put(nodes[i], i); // populates the map with the node (key) and its index (value)
		}
		
		int startNodeValue = nodeIndex.get(start); // gets the value of the start node
		int targetNodeValue = nodeIndex.get(target); // gets the value of the target node - useful for future referral
		Set<Integer> vMinusS = new HashSet<>(); // set of all nodes that have not been visited
		
		
		int[] pred = new int[numV]; // stores predecessor of each node
		double[] dist = new double[numV]; // array that keeps track of the shortest distance established from the start node to EACH node. 
		// E.g. Node 20 would have a distance value that corresponds to how far it is from the start node
 		
		for(int i = 0; i < numV; i++) { // for loop to initialise the dist and pred arrays
			if(i == startNodeValue) {
				dist[i] = 0; // distance to the start node is obviously 0
			}
			else {
				dist[i] = Double.POSITIVE_INFINITY; // all unvisited nodes' distances are initially set to infinity
			}
			pred[i] = -1; // no predecessors yet (for every node)
			vMinusS.add(i); // each node is added to unvisited
		}
		
		while(!vMinusS.isEmpty()) { // while there are nodes not yet visited, find the node with the smallest distance
			double minDist = Double.POSITIVE_INFINITY;
			int currentIndex = -1;
			for(int v : vMinusS) {
				if(dist[v] < minDist) {
					minDist = dist[v];
					currentIndex = v;
				}
			}
			
			vMinusS.remove(currentIndex); // Current node is removed which marks it as visited
			Node currentNode = nodes[currentIndex]; // obtain the node of the corresponding index
			
			for(Edge edge : currentNode.getEdges()) { // explore all neighbours
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
