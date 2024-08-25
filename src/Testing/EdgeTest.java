package Testing;

import Graph.Edge;
import Graph.Node;
import Graph.Position;

import static junit.framework.TestCase.fail;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * @author James Baumeister 3/5/17
 * @author Daniel Ablett 2022
 * @author Gun Lee 2023
 */
public class EdgeTest extends DSUnitTesting {
	
	private Node node1;
	private Node node2;
	private Edge e;
	
	@Before
	public void initialise() {
		node1 = new Node(10, new Position(5, 5));
		node2 = new Node(50, new Position(23, 10));
	}
	
	@Test
	public void constructor() {
		AssignmentMarker.marks.put("Edge:constructor", 5.0f);
		e = new Edge(node1, node2);
		
		try {
			Assert.assertEquals("Weight to be set in the constructor",
					747.448, e.getWeight(), 0.01);
		}		
		catch (NullPointerException e) {
			fail("Null pointer: weight not set in the constructor");
		}
	}
	
	@Test
	public void calculateWeight() {
		AssignmentMarker.marks.put("Edge:calculateWeight", 5.0f);
		e = new Edge(node1, node2);
		
		Assert.assertEquals("Weight of the edge",
				747.448, e.getWeight(), 0.01);
		
		node2 = new Node(15, new Position(2, 2));
		e = new Edge(node1, node2);
		Assert.assertEquals("Weight of the edge",
				21.256, e.getWeight(), 0.01);
		
		e = new Edge(node2, node1);
		Assert.assertEquals("Weight of the edge",
				21.256, e.getWeight(), 0.01);
	}
}