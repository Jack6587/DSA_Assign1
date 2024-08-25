module dsaAssign1 {
	requires javafx.graphics;
	requires javafx.fxml;
	requires javafx.controls;
	requires junit;
	opens GUI to javafx.graphics, javafx.fxml;
	opens Testing to junit, hamcrest;
}