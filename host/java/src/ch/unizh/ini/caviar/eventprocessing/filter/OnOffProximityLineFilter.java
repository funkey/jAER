/*
 * RetinaBackgrondActivityFilter.java
 *
 * Created on October 21, 2005, 12:33 PM
 *
 * To change this template, choose Tools | Options and locate the template under
 * the Source Creation and Management node. Right-click the template and choose
 * Open. You can then make changes to the template in the Source Editor.
 */

package ch.unizh.ini.caviar.eventprocessing.filter;

import ch.unizh.ini.caviar.chip.*;
import ch.unizh.ini.caviar.event.*;
import ch.unizh.ini.caviar.event.EventPacket;
import ch.unizh.ini.caviar.eventprocessing.EventFilter2D;
import java.util.*;

/**
 * An AE filter that outputs only events that are supported by a nearby event of the opposite polarity
 in the neighborhood. The neighborhood is defined
 * by a subsampling bit shift.
 * @author tobi
 */
public class OnOffProximityLineFilter extends EventFilter2D implements Observer  {
    
    final int DEFAULT_TIMESTAMP=Integer.MIN_VALUE;
    
    /** the time in timestamp ticks (1us at present) that a spike
     * needs to be supported by a prior event in the neighborhood by to pass through
     */
    protected int dt=prefs.getInt("OnOffProximityLineFilter.dt",30000);
    {setPropertyTooltip("dt","Events with less than this delta time to neighbors pass through");}
    
    /** the amount to subsample x and y event location by in bit shifts when writing to past event times
     *map. This effectively increases the range of support. E.g. setting subSamplingShift to 1 quadruples range
     *because both x and y are shifted right by one bit */
    private int subsampleBy=prefs.getInt("OnOffProximityLineFilter.subsampleBy",0);
    {setPropertyTooltip("subsampleBy","Past events are subsampled by this many bits");}
    
    
    int[][][] lastTimestamps;
    
    public OnOffProximityLineFilter(AEChip chip){
        super(chip);
        chip.addObserver(this);
        initFilter();
        resetFilter();
    }
    
    void allocateMaps(AEChip chip){
        lastTimestamps=new int[chip.getSizeX()][chip.getSizeY()][2];
    }
    
    int ts=0; // used to reset filter
    
    /**
     * filters in to out. if filtering is enabled, the number of out may be less
     * than the number put in
     *@param in input events can be null or empty.
     *@return the processed events, may be fewer in number. filtering may occur in place in the in packet.
     */
    synchronized public EventPacket filterPacket(EventPacket in) {
        if(!filterEnabled) return in;
        if(in.getEventClass()!=PolarityEvent.class) {
            log.warning("can only process PolarityEvent, disabling filter");
            setFilterEnabled(false);
            return in;
        }
        checkOutputPacketEventType(in);
        if(lastTimestamps==null) allocateMaps(chip);
        // for each event only write it to the out buffers if it is within dt of the last time an event happened in neighborhood
        OutputEventIterator outItr=out.outputIterator();
        int sx=chip.getSizeX()-1;
        int sy=chip.getSizeY()-1;
        for(Object e:in){
            PolarityEvent i=(PolarityEvent)e;
            ts=i.timestamp;
            byte oppType=i.type==0? (byte)1:(byte)0;
            // subsample space part of address, check if delta t of current event to opposite polarity
            // is within dt. if so, output event. in either case write event to lastTimestamps map.
            short x=(short)(i.x>>>subsampleBy), y=(short)(i.y>>>subsampleBy);
            int lastt=lastTimestamps[x][y][oppType];
            int deltat=(ts-lastt);
            if(deltat<dt && lastt!=DEFAULT_TIMESTAMP){
                PolarityEvent o=(PolarityEvent)outItr.nextOutput();
                o.copyFrom(i);
            }
            lastTimestamps[x][y][i.type]=ts;
        }
        return out;
    }
    
    /**
     * gets the background allowed delay in us
     * @return delay allowed for spike since last in neighborhood to pass (us)
     */
    public int getDt() {
        return this.dt;
    }
    
    /**
     * sets the background delay in us
     * @see #getDt
     * @param dt delay in us
     */
    public void setDt(final int dt) {
        prefs.putInt("OnOffProximityLineFilter.dt",dt);
        support.firePropertyChange("dt",this.dt,dt);
        this.dt = dt;
    }
    
    public Object getFilterState() {
        return lastTimestamps;
    }
    
    void resetLastTimestamps(){
        for(int i=0;i<lastTimestamps.length;i++){
            for(int j=0;j<lastTimestamps[i].length;j++){
                Arrays.fill(lastTimestamps[i][j],DEFAULT_TIMESTAMP);
            }
        }
    }
    
    synchronized public void resetFilter() {
        // set all lastTimestamps to max value so that any event is soon enough, guarenteed to be less than it
        resetLastTimestamps();
    }
    
    
    public void update(Observable o, Object arg) {
//        if(!isFilterEnabled()) return;
        initFilter();
    }
    
    public void initFilter() {
        allocateMaps(chip);
    }
    
    public int getSubsampleBy() {
        return subsampleBy;
    }
    
    /** Sets the number of bits to subsample by when storing events into the map of past events.
     *Increasing this value will increase the number of events that pass through and will also allow
     *passing events from small sources that do not stimulate every pixel.
     *@param subsampleBy the number of bits, 0 means no subsampling, 1 means cut event time map resolution by a factor of two in x and in y
     **/
    public void setSubsampleBy(int subsampleBy) {
        if(subsampleBy<0) subsampleBy=0; else if(subsampleBy>4) subsampleBy=4;
        this.subsampleBy = subsampleBy;
        prefs.putInt("OnOffProximityLineFilter.subsampleBy",subsampleBy);
    }
    
    
}
