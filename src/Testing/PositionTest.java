package Testing;

import Graph.Position;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author James Baumeister on 4/5/17
 */
public class PositionTest extends DSUnitTesting {
	
	@Test
	public void distance() {
		AssignmentMarker.marks.put("Position:distance", 5.0f);
		
		Position v1 = new Position(1, 1);
		Position v2 = new Position(1, 4);
		
		Assert.assertEquals("Euclidean distance calculation", 3.0, v1.distance(v2), 0.01);
		Assert.assertEquals("Euclidean distance calculation", 3.0, v2.distance(v1), 0.01);
		Assert.assertEquals("Euclidean distance calculation", 0.0, v1.distance(v1), 0.01);
		
		v2.x = 4;
		
		Assert.assertEquals("Euclidean distance calculation", 4.24, v1.distance(v2), 0.01);
		
		v2.x = -4;
		
		Assert.assertEquals("Euclidean distance calculation", 5.83, v1.distance(v2), 0.01);
	}
}
