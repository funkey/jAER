/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.unizh.ini.jaer.projects.labyrinth;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;
import java.util.Observable;
import java.util.Observer;
import javax.media.opengl.GL;
import javax.media.opengl.GLAutoDrawable;
import javax.media.opengl.GLCanvas;
import javax.media.opengl.glu.*;
import net.sf.jaer.chip.AEChip;
import net.sf.jaer.event.BasicEvent;
import net.sf.jaer.event.EventPacket;
import net.sf.jaer.eventprocessing.*;
import net.sf.jaer.eventprocessing.filter.*;
import net.sf.jaer.eventprocessing.tracking.RectangularClusterTracker;
import net.sf.jaer.eventprocessing.tracking.RectangularClusterTracker.Cluster;
import net.sf.jaer.graphics.ChipCanvas;
import net.sf.jaer.graphics.FrameAnnotater;

/**
 * Specialized tracker for ball location.
 *
 * @author tobi
 */
public class LabyrinthBallTracker extends EventFilter2D implements FrameAnnotater, Observer {

    public static String getDescription() {
        return "Ball tracker for labyrinth game";
    }
    // filters and filter chain
    FilterChain filterChain;
    private RectangularClusterTracker.Cluster ball = null;
    RectangularClusterTracker tracker;
    LabyrinthMap map;
    
    // starting ball location on reset
    private Point2D.Float startingLocation = new Point2D.Float(getFloat("startingX", 50), getFloat("startingY", 100));
    // private fields, not properties
    BasicEvent startingEvent = new BasicEvent();
    int lastTimestamp = 0;
    GLCanvas glCanvas = null;
    private ChipCanvas canvas;

    public LabyrinthBallTracker(AEChip chip) {
        super(chip);
        filterChain = new FilterChain(chip);
        map=new LabyrinthMap(chip);
        filterChain.add(map);
        filterChain.add(new BackgroundActivityFilter(chip));
//        filterChain.add(new CircularConvolutionFilter(chip));
//        filterChain.add(new SubSamplingBandpassFilter(chip)); // TODO preferences should save enabled state of filters
        filterChain.add((tracker = new RectangularClusterTracker(chip)));
        tracker.addObserver(this);
        setEnclosedFilterChain(filterChain);
        String s = " Labyrinth Tracker";
        setPropertyTooltip(s, "startingLocation", "pixel location of starting location for tracker cluster");
        if (chip.getCanvas() != null && chip.getCanvas().getCanvas() != null) {
            glCanvas = (GLCanvas) chip.getCanvas().getCanvas();
        }
    }

    @Override
    public EventPacket<?> filterPacket(EventPacket<?> in) {
        out = getEnclosedFilterChain().filterPacket(in);
        if (tracker.getNumClusters() > 0) {
            // find most likely ball cluster from all the clusters. This is the one with most mass.
            float max = Float.MIN_VALUE;
            for (Cluster c : tracker.getClusters()) {
                if (!c.isVisible()) {
                    continue;
                }
                Point2D.Float l = c.getLocation();
                final int b = 5;
                if (l.x < -b || l.x > chip.getSizeX() + b || l.y < -b || l.y > chip.getSizeY() + b) {
                    continue;
                }
                float mass = c.getMass();
                if (mass > max) {
                    max = mass;
                    ball = c;
                }
            }
        } else {
            ball = null;
        }
        if (!in.isEmpty()) {
            lastTimestamp = in.getLastTimestamp();
        }
        return out;
    }

    @Override
    public void resetFilter() {
        createBall(startingLocation);
    }

    protected void createBall(Point2D.Float location) {
        getEnclosedFilterChain().reset();
        // TODO somehow need to spawn an initial cluster at the starting location
        Cluster b = tracker.createCluster(new BasicEvent(lastTimestamp, (short) location.x, (short) location.y, (byte) 0));
        b.setMass(10000); // some big number
        tracker.getClusters().add(b);
    }

    @Override
    public void initFilter() {
    }
    private GLU glu = new GLU();
    private GLUquadric quad = null;

    @Override
    public void annotate(GLAutoDrawable drawable) {
        GL gl = drawable.getGL();

        // annotate starting ball location
        gl.glColor3f(0, 0, 1);
        gl.glLineWidth(3f);
        gl.glBegin(GL.GL_LINES);
        final int CROSS_SIZE = 3;
        gl.glVertex2f(startingLocation.x - CROSS_SIZE, startingLocation.y - CROSS_SIZE);
        gl.glVertex2f(startingLocation.x + CROSS_SIZE, startingLocation.y + CROSS_SIZE);
        gl.glVertex2f(startingLocation.x - CROSS_SIZE, startingLocation.y + CROSS_SIZE);
        gl.glVertex2f(startingLocation.x + CROSS_SIZE, startingLocation.y - CROSS_SIZE);
        gl.glEnd();

        if (ball != null) {
            gl.glColor4f(.25f, .25f, .25f, .3f);
            gl.glPushMatrix();
            gl.glTranslatef(ball.location.x, ball.location.y, 0);
            if (glu == null) {
                glu = new GLU();
            }
            if (quad == null) {
                quad = glu.gluNewQuadric();
            }
            glu.gluQuadricDrawStyle(quad, GLU.GLU_LINE);
            glu.gluDisk(quad, 0, ball.getRadius(), 16, 1);
            gl.glPopMatrix();
        }

    }

    /**
     * @return the ball
     */
    public RectangularClusterTracker.Cluster getBall() {
        return ball;
    }

    /**
     * @return the startingLocation
     */
    public Point2D.Float getStartingLocation() {
        return startingLocation;
    }

    /**
     * @param startingLocation the startingLocation to set
     */
    public void setStartingLocation(Point2D.Float startingLocation) {
        this.startingLocation = startingLocation;
        putFloat("startingX", startingLocation.x);
        putFloat("startingY", startingLocation.y);
    }

    @Override
    public void update(Observable o, Object arg) {
        if (arg instanceof UpdateMessage) {
            UpdateMessage m = (UpdateMessage) arg;
            callUpdateObservers(m.packet, m.timestamp); // pass on updates from tracker
        }
    }

    /** Resets tracker and creates a ball at the specified location with large initial mass.
     * 
     * @param pf the starting location. 
     */
    public void setBallLocation(Point2D.Float pf) {
        resetFilter();
        createBall(pf);
    }
}