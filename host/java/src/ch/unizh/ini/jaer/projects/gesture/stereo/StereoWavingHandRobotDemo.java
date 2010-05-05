/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.unizh.ini.jaer.projects.gesture.stereo;
import ch.unizh.ini.jaer.hardware.pantilt.PanTilt;
import ch.unizh.ini.jaer.projects.gesture.virtualdrummer.BlurringFilter2DTracker.Cluster;
import java.awt.geom.*;
import java.util.List;
import java.util.Observable;
import java.util.Observer;
import javax.media.opengl.GL;
import javax.media.opengl.GLAutoDrawable;
import net.sf.jaer.chip.AEChip;
import net.sf.jaer.event.EventPacket;
import net.sf.jaer.eventprocessing.EventFilter2D;
import net.sf.jaer.eventprocessing.FilterChain;
import net.sf.jaer.graphics.FrameAnnotater;
import net.sf.jaer.hardwareinterface.*;
import net.sf.jaer.util.filter.LowpassFilter2d;
/**
 *  Simple demonstration of stereo disparity using pantilt servo to wave hand to show output of stereo blob tracking.
 * @author tobi delbruck, junheang lee
 *
 * This is part of jAER
<a href="http://jaer.wiki.sourceforge.net">jaer.wiki.sourceforge.net</a>,
licensed under the LGPL (<a href="http://en.wikipedia.org/wiki/GNU_Lesser_General_Public_License">http://en.wikipedia.org/wiki/GNU_Lesser_General_Public_License</a>.
 */
public class StereoWavingHandRobotDemo extends EventFilter2D implements FrameAnnotater,Observer{
    private float servoLimit = getPrefs().getFloat("StereoWavingHandRobotDemo.servoLimit",0.25f);
    private float disparityScaling = getPrefs().getFloat("StereoWavingHandRobotDemo.disparityScaling",0.05f);
    private float tauArmMs=getPrefs().getFloat("StereoWavingHandRobotDemo.tauArmMs",20);
    private FilterChain chain = null;
    BlurringFilterStereoTracker tracker = null;
    PanTilt panTilt = null;
    LowpassFilter2d filter=new LowpassFilter2d(tauArmMs);
    private final int SERVO_RETRY_INTERVAL=300;

    public StereoWavingHandRobotDemo (AEChip chip){
        super(chip);
        setPropertyTooltip("servoLimit","limits servos to this value around 0.5f");
        setPropertyTooltip("disparityScaling","disparity values are scaled by this factor to arrive at servo position");
        tracker = new BlurringFilterStereoTracker(chip);
        tracker.addObserver(this);
        chain=new FilterChain(chip);
        chain.add(tracker);
        setEnclosedFilterChain(chain);
    }
    @Override
    public EventPacket<?> filterPacket (EventPacket<?> in){
        in = tracker.filterPacket(in);
//        moveArm();  // rely on updates from tracker at updateIntervalMs to move arm
        return in;
    }

    @Override
    public void resetFilter (){
        setPanTilt(.5f,.5f);
        filter.setInternalValue2d(.5f,.5f);
    }

    private float clipServo(float f){
        if(f<0.5f-servoLimit) f=.5f-servoLimit; else if(f>.5f+servoLimit) f=.5f+servoLimit;
        return f;
    }

    private int servoRetryCounter=0;

    private void setPanTilt(float pan, float tilt){
        pan=clipServo(pan);
        tilt=clipServo(tilt);
        if(panTilt!=null && servoRetryCounter<0){ // don't keep retrying every time
            try{
                panTilt.setPanTiltValues(pan,tilt);
            }catch(HardwareInterfaceException e){
                log.warning(e.toString());
                servoRetryCounter=SERVO_RETRY_INTERVAL;
            }
        }
        servoRetryCounter--;
    }

    @Override
    public void initFilter (){
        panTilt=new PanTilt();
        resetFilter();
    }

    public void annotate (GLAutoDrawable drawable){
        GL gl=drawable.getGL();
        gl.glPushMatrix();
        gl.glColor3f(0,0,1);
        gl.glLineWidth(5);
        Point2D.Float p=filter.getValue2d();
        final int W=10;
        float x=p.x*chip.getSizeX();
        float y=p.y*chip.getSizeY();
        gl.glRectf(x-W,y-W,x+W,y+W);

        gl.glPopMatrix();
    }

    @Override
    public synchronized void setFilterEnabled (boolean yes){
        super.setFilterEnabled(yes);
        if ( yes ){
            if ( panTilt == null ){
                panTilt = new PanTilt();
            }
        } else{
            if ( panTilt != null ){
                panTilt.close();
            }
        }
    }

    private Cluster findCluster (){
        List<Cluster> clusters = tracker.getClusters();
        Cluster sc = null;
        for ( Cluster c:clusters ){
          
            if ( c.isVisible() ){
                sc = c;
                break;
            }
        }
        return sc;
    }

    private void moveArm (int timestamp){
        Cluster sc = findCluster();
        if ( sc != null ){
//            float disp = sc.getDisparity();
            Point2D.Float p = sc.getLocation();
//            Point3D p3 = sc.getLocation3dm();
            float pan=(float)p.getX()/chip.getSizeX();
            float tilt=(float)p.getY()/chip.getSizeY();
            Point2D.Float pt=filter.filter2d(pan,tilt,timestamp);
//            float tilt=sc.getDisparity()*disparityScaling;
            setPanTilt(pt.x,pt.y);
        }
    }

    synchronized public void update (Observable o,Object arg){
        if ( arg instanceof UpdateMessage ){
            UpdateMessage msg = (UpdateMessage)arg;
            moveArm(msg.timestamp);
//            log.info("update from "+o+" with message "+arg.toString());
        }
    }

    /**
     * @return the servoLimit
     */
    public float getServoLimit (){
        return servoLimit;
    }

    /**
     * @param servoLimit the servoLimit to set
     */
    public void setServoLimit (float servoLimit){
        servoLimit=servoLimit<0?0:servoLimit;
        servoLimit=servoLimit>1?1:servoLimit;
        this.servoLimit = servoLimit;
        getPrefs().putFloat("StereoWavingHandRobotDemo.servoLimit",servoLimit);
    }

    /**
     * @return the disparityScaling
     */
    public float getDisparityScaling (){
        return disparityScaling;
    }

    /**
     * @param disparityScaling the disparityScaling to set
     */
    public void setDisparityScaling (float disparityScaling){
        this.disparityScaling = disparityScaling;
        getPrefs().putFloat("StereoWavingHandRobotDemo.disparityScaling",disparityScaling);
    }

    /**
     * @return the tauMs
     */
    public float getTauArmMs (){
        return tauArmMs;
    }

    /**
     * The lowpass time constant of the arm.
     *
     * @param tauMs the tauMs to set
     */
    public void setTauArmMs (float tauMs){
        this.tauArmMs = tauMs;
        getPrefs().putFloat("StereoWavingHandRobotDemo.tauArmMs",tauMs);
        filter.setTauMs(tauArmMs);
    }
}