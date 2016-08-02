import java.applet.*;
import java.awt.*;

// Applet code for the "Hello, world!" example.
// This should be saved in a file named as "HelloWorld.java".
public class HelloWorld extends Applet {

    // Print a message on the screen (x=20, y=10).
    public void paint(Graphics g) {
        g.drawString("Hello, world!", 20, 10);

        // Draws a circle on the screen (x=40, y=30).
        g.drawArc(40, 30, 20, 20, 0, 360);

      // Draws a rectangle on the screen (x1=100, y1=100, x2=300,y2=300).
        g.drawRect(100, 100, 300, 300);

      // Draws a square on the screen (x1=100, y1=100, x2=200,y2=200).
        g.drawRect(100, 100, 200, 200);



    }
}


