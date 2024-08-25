package Testing;

import Graph.Graph;
import Graph.Node;
import Graph.Edge;
import GUI.GraphLoader;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static junit.framework.TestCase.fail;


/**
 * @author James Baumeister 4/5/17
 * @author Daniel Ablett 2022
 * @author Gun Lee 2023
 * Most tests assume that connectNodes is working.
 * Complete that method first.
 */
public class GraphTest extends DSUnitTesting {
	
	Graph g;
	Node[] nodes;
	
	@Before
	public void initialise() {
		GraphLoader gl = new GraphLoader();
		nodes = gl.getNodes();
		g = new Graph();
	}
	
	@Test
	public void connectNodes1() {
		AssignmentMarker.marks.put("Graph:connectNodes1", 5.0f);
		
		for (Node n : nodes) {
			if (n.getEdges().size() != 0) {
				fail("Nodes should not be connected before connectNodes() is called.");
			}
		}
		
		g.connectNodes(nodes);
		
		for (Node n : nodes) {
			if (n.getEdges().size() <= 0) {
				fail("Nodes are not connected.");
			}
		}
	}
	
	@Test
	public void connectNodes2() {
		AssignmentMarker.marks.put("Graph:connectNodes2", 5.0f);
				
		g.connectNodes(nodes);
		
		List<Node> nl = Arrays.asList(nodes);
				
		// Test a couple of cases
		Assert.assertEquals("Node 0 should have three edges", 3, nodes[0].getEdges().size());
		Assert.assertEquals("Node 0's east edge", 1, nl.indexOf(nodes[0].getEdges().get(0).getToNode()) );
		Assert.assertEquals("Node 0's south edge", 10, nl.indexOf(nodes[0].getEdges().get(1).getToNode()) );
		Assert.assertEquals("Node 0's south-east edge", 11, nl.indexOf(nodes[0].getEdges().get(2).getToNode()) );
		
		Assert.assertEquals("Node 99 should have three edges", 3, nodes[99].getEdges().size());
		Assert.assertEquals("Node 99's west edge", 98, nl.indexOf(nodes[99].getEdges().get(0).getToNode()) );
		Assert.assertEquals("Node 99's north edge", 89, nl.indexOf(nodes[99].getEdges().get(1).getToNode()) );
		Assert.assertEquals("Node 99's north-west edge", 88, nl.indexOf(nodes[99].getEdges().get(2).getToNode()) );
		
		Assert.assertEquals("Node 45 should have 8 edges", 8, nodes[45].getEdges().size());
		Assert.assertEquals("Node 45's east edge", 46, nl.indexOf(nodes[45].getEdges().get(0).getToNode()) );
		Assert.assertEquals("Node 45's west edge", 44, nl.indexOf(nodes[45].getEdges().get(1).getToNode()) );
		Assert.assertEquals("Node 45's south edge", 55, nl.indexOf(nodes[45].getEdges().get(2).getToNode()) );
		Assert.assertEquals("Node 45's north edge", 35, nl.indexOf(nodes[45].getEdges().get(3).getToNode()) );
		Assert.assertEquals("Node 45's north-east edge", 36, nl.indexOf(nodes[45].getEdges().get(4).getToNode()) );
		Assert.assertEquals("Node 45's south-east edge", 56, nl.indexOf(nodes[45].getEdges().get(5).getToNode()) );
		Assert.assertEquals("Node 45's north-west edge", 34, nl.indexOf(nodes[45].getEdges().get(6).getToNode()) );
		Assert.assertEquals("Node 45's south-west edge", 54, nl.indexOf(nodes[45].getEdges().get(7).getToNode()) );
		
		Assert.assertEquals("Node 67 should have 8 edges", 8, nodes[67].getEdges().size());
		Assert.assertEquals("Node 67's east edge", 68, nl.indexOf(nodes[67].getEdges().get(0).getToNode()) );
		Assert.assertEquals("Node 67's west edge", 66, nl.indexOf(nodes[67].getEdges().get(1).getToNode()) );
		Assert.assertEquals("Node 67's south edge", 77, nl.indexOf(nodes[67].getEdges().get(2).getToNode()) );
		Assert.assertEquals("Node 67's north edge", 57, nl.indexOf(nodes[67].getEdges().get(3).getToNode()) );
		Assert.assertEquals("Node 67's north-east edge", 58, nl.indexOf(nodes[67].getEdges().get(4).getToNode()) );
		Assert.assertEquals("Node 67's south-east edge", 78, nl.indexOf(nodes[67].getEdges().get(5).getToNode()) );
		Assert.assertEquals("Node 67's north-west edge", 56, nl.indexOf(nodes[67].getEdges().get(6).getToNode()) );
		Assert.assertEquals("Node 67's south-west edge", 76, nl.indexOf(nodes[67].getEdges().get(7).getToNode()) );
	}
	
	@Test
	public void getEdge1() {
		AssignmentMarker.marks.put("Graph:getEdge1", 3.0f);
		Node n1 = nodes[0];
		Node n2 = nodes[1];
		
		// Some fake edges
		Edge e1 = new Edge(nodes[4], nodes[6]);
		Edge e2 = new Edge(nodes[10], nodes[11]);
		Edge e3 = new Edge(nodes[3], nodes[2]);
		Edge e4 = new Edge(nodes[8], n2);
		
		// The actual edge
		Edge e = new Edge(n1, n2);
		
		ArrayList<Edge> edges1 = new ArrayList<Edge>();
		edges1.add(e1);
		edges1.add(e2);
		edges1.add(e3);
		edges1.add(e4);
		edges1.add(e);
		n1.setEdges(edges1);		
		
		Assert.assertEquals("Correct edge", e, g.getEdge(n1, n2));
	}
	
	@Test
	public void getEdge2() {
		AssignmentMarker.marks.put("Graph:getEdge2", 2.0f);
		Node n1 = nodes[0];
		Node n2 = nodes[1];
				
		// The actual edge
		Edge e = new Edge(n1, n2);
		
		ArrayList<Edge> edges1 = new ArrayList<Edge>();
		edges1.add(e);
		n1.setEdges(edges1);
				
		Assert.assertEquals("Correct edge", e, g.getEdge(n1, n2));
		Assert.assertEquals("Non-existing edge", null, g.getEdge(n2, n1));
	}
	
	@Test
	public void calculateTotalWeight1() {
		AssignmentMarker.marks.put("Graph:calculateTotalWeight1", 5.0f);
		g.connectNodes(nodes);
		
		Node[] path = {nodes[23], nodes[24], nodes[35], nodes[46], nodes[57]};
		
		try {
			Assert.assertEquals("Cost of a path.", 71.028, g.calculateTotalWeight(path),
					0.01);
		}
		catch (NullPointerException e) {
			fail("Null pointer: connectNodes may not be functioning correctly.");
		}
	}
	
	@Test
	public void calculateTotalWeight2() {
		AssignmentMarker.marks.put("Graph:calculateTotalWeight2", 5.0f);
		g.connectNodes(nodes);
		
		// This is invalid path with no edge connecting, hence should not return the weight of non-existing edge 
		Node[] invalid_path = {nodes[23], nodes[46]};
		
		try {
			Assert.assertEquals("Cost of a path with non-existing edge expected not to be calculated (i.e., return infinity)", Double.POSITIVE_INFINITY, g.calculateTotalWeight(invalid_path), 0.001);
		}
		catch (NullPointerException e) {  
			fail("Null pointer: not handling non-existing edge.");
		}
	}
	
	
	
	// BFS 
	
	@Test
	public void breadthFirstSearch1() {
		AssignmentMarker.marks.put("Graph:breadthFirstSearch1", 5.0f);
		g.connectNodes(nodes);
		
		Node[] exp_path = {nodes[0], nodes[10], nodes[20], nodes[30], nodes[40], nodes[50], nodes[60]};
		Node[] actual_path = g.breadthFirstSearch(nodes[0], nodes[60]);
		
		if (!paths_equal(exp_path, actual_path))
			fail("BFS1: incorrect path.");
	}
	
	@Test
	public void breadthFirstSearch2() {
		AssignmentMarker.marks.put("Graph:breadthFirstSearch2", 5.0f);
		g.connectNodes(nodes);
		
		Node[] exp_path = new Node[]{nodes[71], nodes[61], nodes[51], nodes[42], nodes[33], nodes[24], nodes[15]};
		Node[] actual_path = g.breadthFirstSearch(nodes[71], nodes[15]);
		
		if (!paths_equal(exp_path, actual_path))
			fail("BFS2: incorrect path.");
	}
	
	@Test
	public void breadthFirstSearch3() {
		AssignmentMarker.marks.put("Graph:breadthFirstSearch3", 5.0f);
		g.connectNodes(nodes);
		
		Node[] exp_path = new Node[]{nodes[3], nodes[4], nodes[5], nodes[6], nodes[17]};
		Node[] actual_path = g.breadthFirstSearch(nodes[3], nodes[17]);
		
		if (!paths_equal(exp_path, actual_path))
			fail("BFS3: incorrect path.");
	}
	
	
	
	
	private boolean paths_equal(Node[] exp, Node[] act) {
		if (exp.length != act.length) {
			print_path("Expected: " , exp);
			print_path("Actual  : " , act);
			return false;//fail(test_num + ": Paths are of different lengths.");
		}
		
		for (int i = 0; i < exp.length; i++) {
			if (exp[i] != act[i]) {
				print_path("Expected: " , exp);
				print_path("Actual  : " , act);
				return false;//fail(test_num + ": Actual path is different than expected.");
			}
		}
		
		return true;
	}
	
	private void print_path(String header, Node[] path) {
		System.out.print(header);
		List<Node> ns = Arrays.asList(nodes);
		for (Node n : path)
			System.out.print(ns.indexOf(n) + " ");
		System.out.println();
	}
	
	
	
	
	// DFS
	
	@Test
	public void depthFirstSearch1() {
		AssignmentMarker.marks.put("Graph:depthFirstSearch1", 5.0f);
		g.connectNodes(nodes);
		
		Node[] exp_path = {nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7],
				nodes[8], nodes[9], nodes[19], nodes[18]};
		Node[] actual_path = g.depthFirstSearch(nodes[0], nodes[18]);
		
		if (!paths_equal(exp_path, actual_path))
			fail("DFS1: incorrect path.");
	}
	
	@Test
	public void depthFirstSearch2() {
		AssignmentMarker.marks.put("Graph:depthFirstSearch2", 5.0f);
		g.connectNodes(nodes);
		
		Node[] exp_path = new Node[]{nodes[63], nodes[64], nodes[65], nodes[66], nodes[67], nodes[68], nodes[69], nodes[79],
				nodes[78], nodes[77], nodes[76], nodes[75], nodes[74], nodes[73], nodes[72], nodes[71], nodes[70],
				nodes[60], nodes[61], nodes[62], nodes[52]};
		Node[] actual_path = g.depthFirstSearch(nodes[63], nodes[52]);
		
		if (!paths_equal(exp_path, actual_path))
			fail("DFS2: incorrect path.");
	}
	
	@Test
	public void depthFirstSearch3() {
		AssignmentMarker.marks.put("Graph:depthFirstSearch3", 5.0f);
		g.connectNodes(nodes);
		
		Node[] exp_path = new Node[]{nodes[40], nodes[41], nodes[42], nodes[43], nodes[44], nodes[45], nodes[46], nodes[47],
				nodes[48], nodes[49], nodes[59], nodes[58]};
		Node[] actual_path = g.depthFirstSearch(nodes[40], nodes[58]);
		
		if (!paths_equal(exp_path, actual_path))
			fail("DFS3: incorrect path.");
	}
	
	
	
	// Dijkstra	
	
	@Test
	public void dijkstrasSearch1() {
		AssignmentMarker.marks.put("Graph:dijkstrasSearch1", 5.0f);
		g.connectNodes(nodes);

		Node[] exp_path = {nodes[40], nodes[50], nodes[60], nodes[71], nodes[72], nodes[73], nodes[74],
				nodes[75], nodes[65], nodes[55], nodes[46], nodes[47], nodes[58]};
		Node[] actual_path = g.dijkstrasSearch(nodes[40], nodes[58]);

		if (!paths_equal(exp_path, actual_path))
			fail("Dijkstra1: incorrect path.");
	}
	
	
	@Test
	public void dijkstrasSearch2() {
		AssignmentMarker.marks.put("Graph:dijkstrasSearch2", 5.0f);
		g.connectNodes(nodes);

		Node[] exp_path = new Node[]{nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[14], nodes[15],
				nodes[16], nodes[26], nodes[37], nodes[28], nodes[18]};
		Node[] actual_path = g.dijkstrasSearch(nodes[0], nodes[18]);

		if (!paths_equal(exp_path, actual_path))
			fail("Dijkstra2: incorrect path.");
	}
	
	@Test
	public void dijkstrasSearch3() {
		AssignmentMarker.marks.put("Graph:dijkstrasSearch3", 5.0f); // lower marks since alternative choices are available if not strict
		g.connectNodes(nodes);

		Node[] exp_path = new Node[]{nodes[96], nodes[95], nodes[84], nodes[83], nodes[82], nodes[72], nodes[71],
				nodes[60], nodes[50], nodes[40], nodes[30]};
		Node[] alt_path = new Node[]{nodes[96], nodes[95], nodes[84], nodes[83], nodes[82], nodes[81], nodes[71],
				nodes[60], nodes[50], nodes[40], nodes[30]};
		Node[] actual_path = g.dijkstrasSearch(nodes[96], nodes[30]);

		if (!paths_equal(exp_path, actual_path) && !paths_equal(alt_path, actual_path))
			fail("Dijkstra3: incorrect path.");
	}
	
	@Test
	public void dijkstrasSearch4() {
		AssignmentMarker.marks.put("Graph:dijkstrasSearch4", 5.0f);
		g.connectNodes(nodes);

		Node[] exp_path = new Node[]{nodes[71], nodes[60], nodes[50], nodes[40], nodes[30], nodes[20], nodes[11], nodes[2],
				nodes[3], nodes[4], nodes[14], nodes[15]};
		Node[] actual_path = g.dijkstrasSearch(nodes[71], nodes[15]);

		if (!paths_equal(exp_path, actual_path))
			fail("Dijkstra4: incorrect path.");
		
	}
	
	
}
