package GUI;

import javafx.fxml.FXML;
import javafx.event.ActionEvent;
import javafx.scene.control.Button;
import javafx.scene.control.ToggleButton;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;

import static GUI.GraphGUI.canvas;

/**
 * This class serves as the main controller for the assignment GUI.
 * All the logic associated with the GUI elements is specified here.
 * @author James Baumeister
 * @version 1.0
 */
public class GraphGUIController {

	@FXML
	private ImageView mapImageView;

	@FXML
	private Button elevation, polygons, bfs, dfs, dijkstra;

	@FXML
	private ToggleButton nodes;
	
	/*
	GraphGUIController() {
		new GraphRenderer();
	}*/

	/**
	 * Handles the pressing of all the image-switching buttons in the
	 * tool bar.
	 * <p>
	 * The event is queried to find the pressed button and the image
	 * is updated accordingly.
	 * @param event the originating event
	 */
	@FXML
	protected void handleMapTypeButtonAction(ActionEvent event) {
		Button b = (Button) event.getSource();

		if (b == polygons) {
			mapImageView.setImage(new Image("resources/simple_map_polygons.png"));
		}
		else if (b == elevation) {
			mapImageView.setImage(new Image("resources/simple_map_elevations.png"));
		}
	}

	/**
	 * Handles the pressing of the node visibility toggle button in
	 * the tool bar.
	 * <p>
	 * If the toggle is selected on press, visibility is turned off,
	 * otherwise visibility is turned on.
	 * @param event the originating event
	 */
	@FXML
	protected void handleToggleNodesButtonAction(ActionEvent event) {
		if (nodes.isSelected()) {
			canvas.setVisible(false);
		}
		else {
			canvas.setVisible(true);
		}

	}
	
	/**
	 * Handles pressing BFS button in the tool bar.
	 * @param event
	 */
	@FXML
	protected void handleBFS(ActionEvent event) {
	}
	
	/**
	 * Handles pressing DFS button in the tool bar.
	 * @param event
	 */
	@FXML
	protected void handleDFS(ActionEvent event) {
	}
	
	/**
	 * Handles pressing BFS button in the tool bar.
	 * @param event
	 */
	@FXML
	protected void handleDijkstra(ActionEvent event) {
	}
}
