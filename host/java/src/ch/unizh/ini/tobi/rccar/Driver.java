/*
 * Driver.java
 *
 * Created on February 27, 2007, 9:51 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 *
 *
 *Copyright February 27, 2007 Tobi Delbruck, Inst. of Neuroinformatics, UNI-ETH Zurich
 */

package ch.unizh.ini.tobi.rccar;

import ch.unizh.ini.caviar.chip.*;
import ch.unizh.ini.caviar.event.*;
import ch.unizh.ini.caviar.event.EventPacket;
import ch.unizh.ini.caviar.eventprocessing.*;
import ch.unizh.ini.caviar.eventprocessing.label.*;
import ch.unizh.ini.caviar.eventprocessing.label.SimpleOrientationFilter;
import ch.unizh.ini.caviar.eventprocessing.tracking.HoughLineTracker;
import ch.unizh.ini.caviar.graphics.FrameAnnotater;
import ch.unizh.ini.caviar.hardwareinterface.*;
import ch.unizh.ini.caviar.util.filter.LowpassFilter;
import java.awt.Graphics2D;
import java.util.logging.*;
import java.util.prefs.*;
import javax.media.opengl.*;
import javax.media.opengl.GLAutoDrawable;
import javax.media.opengl.glu.*;

/**
 * Drives the RC car by either centering event activity or non activity, depending on flipSteering switch.
 
 * @author tobi
 */
public class Driver extends EventFilter2D implements FrameAnnotater{
    
    static Logger log=Logger.getLogger("Driver");
    static Preferences prefs=Preferences.userNodeForPackage(Driver.class);
    private SiLabsC8051F320_USBIO_CarServoController servo;
    private LowpassFilter filter=new LowpassFilter();
    private float gain=prefs.getFloat("Driver.gain",1);
    private float lpCornerFreqHz=prefs.getFloat("Driver.lpCornerFreqHz",1);
    private boolean flipSteering=prefs.getBoolean("Driver.flipSteering",false);
    private HoughLineTracker houghLineTracker;
    private float steerInstantaneous=0.5f; // instantaneous value, before filtering
    private float steerCommand=0.5f; // actual command, as modified by filtering
    private float speed;
    private int sizex;
    private final int STEERING_SERVO=0, SPEED_SERVO=1;
    private float radioSteer=0.5f, radioSpeed=0.5f;
    
    /** Creates a new instance of Driver */
    public Driver(AEChip chip) {
        super(chip);
        chip.getCanvas().addAnnotator(this);
        initFilter();
    }
    
    EventPacket oriPacket=null, outOri=new EventPacket<OrientationEvent>(OrientationEvent.class);
    
    public EventPacket<?> filterPacket(EventPacket<?> in) {
        if(!isFilterEnabled()) return in;
        checkServo();
        houghLineTracker.filterPacket(in);
        
        int n=in.getSize();
        if(n==0) return in;
        
        // get values from radio receiver (user sets speed or steers)
        radioSteer=servo.getRadioSteer();
        radioSpeed=servo.getRadioSpeed();

        sizex=getChip().getSizeX();// must do this here in case chip has changed
        // compute instantaneous position of line according to hough line tracker (which has its own lowpass filter)
        double rhoPixels=(float)houghLineTracker.getRhoPixelsFiltered();  // distance of line from center of image
        double thetaRad=(float)houghLineTracker.getThetaDegFiltered()/180*Math.PI; // angle of line, pi/2 is horizontal
        double hDistance=rhoPixels*Math.cos(thetaRad); // horizontal distance of line from center in pixels
        steerInstantaneous=(float)(hDistance/sizex); // as fraction of image
        
        // apply proportional gain setting, reduce by speed of car, center at 0.5f
        steerInstantaneous=(steerInstantaneous*(1-radioSpeed))*gain+0.5f; 
        steerCommand=filter.filter(steerInstantaneous,in.getLastTimestamp()); // lowpass filter

        if(servo.isOpen()){
            servo.setSteering(steerCommand); // 1 steer right, 0 steer left
        }
        return in;
    }
    
    private void checkServo(){
        if(servo==null){
            servo=new SiLabsC8051F320_USBIO_CarServoController();
        }
    }
    
    public Object getFilterState() {
        return null;
    }
    
    public void resetFilter() {
    }
    
    public void initFilter() {
        filter.set3dBFreqHz(lpCornerFreqHz);
        houghLineTracker=(HoughLineTracker)(chip.getFilterChain().findFilter(HoughLineTracker.class));
        
        setEnclosedFilter(houghLineTracker);
    }
    
    public void annotate(float[][][] frame) {
    }
    
    public void annotate(Graphics2D g) {
    }
    
    GLU glu=null;
    GLUquadric wheelQuad;
    
    
    public void annotate(GLAutoDrawable drawable) {
        if(!isFilterEnabled()) return;
        
        houghLineTracker.annotate(drawable);
        
        GL gl=drawable.getGL();
        if(gl==null) return;
        final int radius=30;
        
        // draw steering wheel
        if(glu==null) glu=new GLU();
        if(wheelQuad==null) wheelQuad = glu.gluNewQuadric();
        gl.glPushMatrix();
        {
            gl.glTranslatef(chip.getSizeX()/2, (chip.getSizeY())/2,0);
            gl.glLineWidth(6f);
            glu.gluQuadricDrawStyle(wheelQuad,GLU.GLU_FILL);
            glu.gluDisk(wheelQuad,radius,radius+1,16,1);
        }
        gl.glPopMatrix();
        
        // draw steering vector, including external radio input value
        
        gl.glPushMatrix();
        {
            gl.glColor3f(1,1,1);
            gl.glTranslatef(chip.getSizeX()/2,chip.getSizeY()/2,0);
            gl.glLineWidth(6f);
            gl.glBegin(GL.GL_LINES);
            {
                gl.glVertex2f(0,0);
                double a=2*(filter.getValue()-0.5f); // -1 to 1
                a=Math.atan(a);
                float x=radius*(float)Math.sin(a);
                float y=radius*(float)Math.cos(a);
                gl.glVertex2f(x,y);
                if(servo!=null && servo.isOpen()){
                    gl.glColor3f(1,0,0);
                    gl.glVertex2f(0,0);
                    a=2*(radioSteer-0.5f); // -1 to 1
                    a=Math.atan(a);
                    x=radius*(float)Math.sin(a);
                    y=radius*(float)Math.cos(a);
                    gl.glVertex2f(x,y);
                }
            }
            gl.glEnd();
        }
        gl.glPopMatrix();
        
        // draw external speed value
        if(servo!=null && servo.isOpen()){
            gl.glPushMatrix();
            {
                gl.glColor3f(1,1,1);
                gl.glTranslatef(1,chip.getSizeY()/2,0);
                gl.glLineWidth(15f);
                gl.glBegin(GL.GL_LINES);
                {
                    gl.glVertex2f(0,0);
                    gl.glVertex2f(0,chip.getSizeY()*(radioSpeed-0.5f));
                }
                gl.glEnd();
            }
            gl.glPopMatrix();
        }
        
        
    }
    
    public float getGain() {
        return gain;
    }
    
    /** Sets steering gain */
    public void setGain(float gain) {
        this.gain = gain;
        prefs.putFloat("Driver.gain",gain);
    }
    
    public float getLpCornerFreqHz() {
        return lpCornerFreqHz;
    }
    
    public void setLpCornerFreqHz(float lpCornerFreqHz) {
        this.lpCornerFreqHz = lpCornerFreqHz;
        prefs.putFloat("Driver.lpCornerFreqHz",lpCornerFreqHz);
        filter.set3dBFreqHz(lpCornerFreqHz);
    }
    
    public boolean isFlipSteering() {
        return flipSteering;
    }
    
    /** If set true, then drive towards events (road is textured), if false, drive away from events (side is textured). */
    public void setFlipSteering(boolean flipSteering) {
        this.flipSteering = flipSteering;
        prefs.putBoolean("Driver.flipSteering",flipSteering);
    }
    
}
