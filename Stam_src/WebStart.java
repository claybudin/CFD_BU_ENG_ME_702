/*
 * WebStart.java
 * Alexander McKenzie
 * 12 March, 2004
 *
 */

/*
    from: http://www.multires.caltech.edu/teaching/demos/java/stablefluids.htm

    NOTES:
        Designed to run in broswer - Chrome won't work, Firefox works well
        Load FluidSolver.html into browser with drag-n-drop
        Don't need to run local web server from dir:
            python -m SimpleHTTPServer 8000
            but it will work as: http://localhost:8000/FluidSolver.html
        Need to add file:/ and/or http://localhost:8000 to Java safe list (Java in Sys Prefs)
        Donwloaded JDK 8 and installed for Mac
            Compile: javac *.java
            Create jar: jar cvf FluidSolver.jar [FW]*.class
        As applet, doesn't need a main() routine
        Seems like java applet gets cached, clearing in Java Control Panel doesn't help, have
          to restart Firefox
*/

import java.applet.Applet;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;

import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

/**
 * Applet display interface for fluid solver.
 *
 * @author Alexander McKenzie
 * @version 1.0
 **/

public class WebStart extends Applet
    implements MouseListener, MouseMotionListener, KeyListener, Runnable
{

    // frame dimensions (dxd pixels)
    int d = 300;

    // solver variables
    int n = 100; //60;
    float dt = 0.2f;
    FluidSolver fs = new FluidSolver();

    // flag to display velocity field
    boolean vkey = false;

    // drawing thread
    Thread artist = null;

    // mouse position
    int x, xOld;
    int y, yOld;

    // cell index
    int i, j;

    // cell dimensions
    int dg, dg_2;

    // cell position
    int dx, dy;

    // fluid velocity
    int u, v;

    int c;
    int button;

    BufferedImage bi;
    Graphics2D big;

    public void reset(){
        // calculate cell deimensions
        dg   = d  / n;
        dg_2 = dg / 2;

        fs.setup(n, dt);
    }

    public void init()
    {
        addMouseMotionListener(this);
        addMouseListener(this);
        addKeyListener(this);
        setFocusable(true);
        bi = (BufferedImage) createImage(d, d);
        big = bi.createGraphics();
    }

    public void start()
    {
        if (artist == null)
        {
            artist = new Thread(this);
            artist.start();
            reset();
        }
    }

    public void stop(){ artist = null; }

    public void run()
    {
        while (artist != null)
        {
            try
            {
                Thread.sleep(20);
            }
            catch (InterruptedException e)
            {
            }
            repaint();
        }
        artist = null;
    }

    public void update(Graphics g){ paint(g); }

    public void paint(Graphics g)
    {
        Graphics2D g2 = (Graphics2D) g;

        // clear screen
        big.setColor(Color.white);
        big.fillRect(0, 0, d, d);

        // solve fluid
        fs.velocitySolver();
        fs.densitySolver();

        for (int i = 1; i <= n; i++)
        {
            // x position of current cell
            dx = (int)( (i - 0.5f) * dg );
            for (int j = 1; j <= n; j++)
            {
                // y position of current cell
                dy = (int)( (j - 0.5f) * dg );

                // draw density
                if (fs.d[I(i, j)] > 0)
                {
                    c = (int) ( (1.0 - fs.d[I(i, j)]) * 255);
                    if (c < 0) c = 0;
                    big.setColor(new Color(c, c, c));
                    big.fillRect(dx-dg_2, dy-dg_2, dg, dg);
                }

                // draw velocity field
                if (vkey && i % 5 == 1 && j % 5 == 1)
                {
                    u = (int)( 50 * fs.u[I(i,j)] );
                    v = (int)( 50 * fs.v[I(i,j)] );
                    big.setColor(Color.red);
                    big.drawLine(dx, dy, dx+u, dy+v);
                }
            }
        }

        // draw status
        big.setColor(new Color(255,153,51));
        big.drawString("Grid: "+n+"x"+n, 5, 15);
        big.drawString("Timestep: "+ dt, 5, 30);
        //big.drawString("TEST", 5, 45);

        g2.drawImage(bi, null, 0, 0);
    }

    public void keyPressed(KeyEvent e)
    {
        // set flag for drawing velocity field
        if (e.getKeyChar() == 'v')
        {
            vkey = !vkey;
        }

        // reset solver
        if (e.getKeyChar() == 'r')
        {
            fs.reset();
        }

        // increase fluid grid size and reset applet
        if (e.getKeyChar() == ']')
        {
            if(n == d) return;

            // calculate next ideal grid size
            int i = n+1;
            while(d%i != 0){
                i++;
            }
            n = i;

            reset();
        }

        // reduce fluid grid size and reset applet
        if (e.getKeyChar() == '[')
        {
            if(n < 10) return;

            // calculate previous ideal grid size
            int i = n-1;
            while(d%i != 0){
                i--;
            }
            n = i;

            reset();
        }

        // increase timestep
        if (e.getKeyChar() == '.')
        {
            if(dt > 1) return;

            dt += 0.05f;

            // kill fp errors
            dt = (float) Math.round(dt * 100);
            dt /= 100;

            fs.dt = dt;
        }

        // reduce timestep
        if (e.getKeyChar() == ',')
        {
            if(dt < 0.1f) return;

            dt -= 0.05f;

            // kill fp errors
            dt = (float) Math.round(dt * 100);
            dt /= 100;

            fs.dt = dt;
        }
    }

    public void mousePressed(MouseEvent e)
    {
        // save button event
        button = e.getButton();

        // update mouse position
        xOld = x;
        yOld = y;
        x = e.getX();
        y = e.getY();

        updateLocation(e);
    }

    public void mouseDragged(MouseEvent e)
    {
        // update mouse position
        xOld = x;
        yOld = y;
        x = e.getX();
        y = e.getY();

        updateLocation(e);
    }

    public void updateLocation(MouseEvent e)
    {
        // get index for fluid cell under mouse position
        i = (int) ((x / (float) d) * n + 1);
        j = (int) ((y / (float) d) * n + 1);

        // set boundries
        if (i > n) i = n;
        if (i < 1) i = 1;
        if (j > n) j = n;
        if (j < 1) j = 1;

        // add density or velocity
        if (button == 1) fs.dOld[I(i, j)] = 100;
        if (button == 3 && e.getID() == MouseEvent.MOUSE_DRAGGED)
        {
            fs.uOld[I(i, j)] = (x - xOld) * 5;
            fs.vOld[I(i, j)] = (y - yOld) * 5;
        }
    }

    // util function for indexing
    private int I(int i, int j){ return i + (n + 2) * j; }

    // fulfill mouse interface requirements
    public void mouseReleased(MouseEvent e){}
    public void mouseMoved(MouseEvent e){}
    public void mouseClicked(MouseEvent e){}
    public void mouseExited(MouseEvent e){}
    public void mouseEntered(MouseEvent e){}

    // fulfill key interface requirements
    public void keyTyped(KeyEvent e){}
    public void keyReleased(KeyEvent e){}
}
