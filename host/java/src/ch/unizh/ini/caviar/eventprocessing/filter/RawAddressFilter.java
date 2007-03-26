/*
 * AbstractEventFilter.java
 *
 * Created on October 30, 2005, 4:58 PM
 *
 * To change this template, choose Tools | Options and locate the template under
 * the Source Creation and Management node. Right-click the template and choose
 * Open. You can then make changes to the template in the Source Editor.
 */

package ch.unizh.ini.caviar.eventprocessing.filter;

import ch.unizh.ini.caviar.aemonitor.*;
import ch.unizh.ini.caviar.eventprocessing.EventFilter;
import java.beans.*;
import java.util.prefs.*;

/**
 * A class that filters only events whose bits are set as in address.
  * @author tobi
 */
public class RawAddressFilter {
    
    protected Preferences prefs=Preferences.userNodeForPackage(EventFilter.class);
    protected PropertyChangeSupport support=new PropertyChangeSupport(this);
    
    /** default true */
    protected boolean filterEnabled=true;
    
    protected short address=0;
    
    AEPacketRaw out=null;
    
    /** Creates a new instance of AbstractEventFilter */
    public RawAddressFilter() {
    }
 
     /**
     * filters in to out. if filtering is enabled, the number of out may be less
     * than the number put in
     *@param in input events can be null or empty.
     *@return the processed events, may be fewer in number. filtering may occur in place in the in packet.
     */
    synchronized public AEPacketRaw filter(AEPacketRaw in) {
        if(!filterEnabled) return in;
        out=new AEPacketRaw();
        short[] a=in.getAddresses();
        int[] t=in.getTimestamps();
        int n=in.getNumEvents();
        for(int i=0;i<n;i++){
            if((a[i]&address)!=0){
                out.addEvent(new EventRaw(a[i],t[i]));
            }
        }
        return out;
    }
    

    /** @return true if filter is enabled */
    public boolean isFilterEnabled() {
        return filterEnabled;
    }

    /** @param enabled true to enable filter. false means output events are the same as input */
    synchronized public void setFilterEnabled(boolean enabled) {
        support.firePropertyChange("filterEnabled",new Boolean(this.filterEnabled),new Boolean(enabled));
        this.filterEnabled=enabled;
        System.out.println("AbstractEventFilter: setFilterEnabled to "+filterEnabled);
    }

    
    public PropertyChangeSupport getPropertyChangeSupport(){
        return support;
    }
    

    public short getAddress() {
        return this.address;
    }

    public void setAddress(final short address) {
        this.address = address;
    }
 }
