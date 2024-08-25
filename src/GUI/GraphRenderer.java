package GUI;

import Graph.*;
import javafx.animation.AnimationTimer;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.paint.Color;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import static GUI.GraphGUI.canvas;
import static GUI.GraphGUI.player_canvas;
import static GUI.Player.position;

/**
 * Renders the simple graph in a simple GUI. The player marker (blue circle)
 * will move along the path returned by one of the search methods in Graph.
 */
class GraphRenderer {

	// For internal testing purposes
	// Change to test different paths
	private final int START_NODE = 0;
	private final int TARGET_NODE = 99;

	// Change this value to increase or decrease animation speed
	private final int ANIMATION_DELAY = 500;

	// Change this string to match the EXACT graph search method you want to animate
	// "breadthFirstSearch", "depthFirstSearch", "dijkstrasSearch"
	private final String METHOD = "breadthFirstSearch"; 

	private GraphicsContext gc;
	private GraphicsContext pgc;
	private Node[] nodes = new GraphLoader().getNodes();

	private Graph graph = new Graph();

	GraphRenderer() {
		gc = canvas.getGraphicsContext2D();
		pgc = player_canvas.getGraphicsContext2D();
		drawNodes();
		drawPlayer(GraphLoader.getNodes()[START_NODE]);
		
		startAnimation(METHOD);
	}
	
	public void startAnimation(String method) {
		new Thread(() -> {
			Method theMethod;
			try {
				Thread.sleep(ANIMATION_DELAY);

				try {
					theMethod = graph.getClass().getMethod(method,
							GraphLoader.getNodes()[START_NODE].getClass(),
							GraphLoader.getNodes()[TARGET_NODE].getClass());
					nodes = (Node[]) theMethod.invoke(graph, GraphLoader.getNodes()[START_NODE], GraphLoader.getNodes()
							[TARGET_NODE]);
				}
				catch (NoSuchMethodException e) {
					System.out.println("Invalid search method");
					System.exit(1);
				}
				catch (InvocationTargetException | SecurityException | IllegalAccessException e) {
					e.printStackTrace();
				}

				if (nodes == null) {
					return;
				}

				for (int i = 0; i < nodes.length; i++) {
					Thread.sleep(ANIMATION_DELAY);
					Player.position = nodes[i];
				}
			}
			catch (InterruptedException e) {
				e.printStackTrace();
			}
		}).start();
	}

	private void drawNodes() {
		gc.setFill(Color.RED);
		gc.setStroke(Color.BLACK);
		gc.setLineWidth(1);

		graph.connectNodes(nodes);

		for (Node n : nodes) {
			for (Edge e : n.getEdges()) {
				gc.strokeLine(n.getPosition().x, n.getPosition().y, e.getToNode().getPosition().x,
						e.getToNode().getPosition().y);

			}
		}

		for (Node n : nodes) {
			gc.fillOval(n.getPosition().x-2, n.getPosition().y-2, 4, 4);
		}
	}

	private void drawPlayer(Node player) {
		new Player(player);
		animatePlayer();
	}

	private void animatePlayer() {
		new AnimationTimer() {

			Position current = new Position(-1, -1);

			public void handle(long currentNanoTime) {
				if (updated()) {
					pgc.clearRect(0, 0, pgc.getCanvas().getWidth(), pgc.getCanvas().getHeight());
					pgc.setFill(Color.AQUAMARINE);
					pgc.fillOval(position.getPosition().x - 8,
							position.getPosition().y - 8,
							16, 16);
					current = position.getPosition();
				}
			}

			private boolean updated() {
				return !position.getPosition().equals(current);
			}
		}.start();
	}
}
