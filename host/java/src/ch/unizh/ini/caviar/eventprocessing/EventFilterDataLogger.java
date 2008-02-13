/*
 * EventFilterDataLogger.java
 *
 * Created on September 7, 2006, 5:02 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 *
 *
 *Copyright September 7, 2006 Tobi Delbruck, Inst. of Neuroinformatics, UNI-ETH Zurich
 */

package ch.unizh.ini.caviar.eventprocessing;

import java.io.*;
import java.util.logging.*;

/**
 * Eases writing data log files for a filter
 
 * @author tobi
 */
public class EventFilterDataLogger {
    
    static Logger log=Logger.getLogger("EventFilterDataLogger");
    protected PrintStream logStream;
    boolean logDataEnabled=false;
    EventFilter2D filter;
    String headerLine;
    
    /**
     * Creates a new instance of EventFilterDataLogger for a filter.
     *@param filter the filter (used to make filename)
     *@param headerLineComment a comment usuually specifying the contents and data fields
     */
    public EventFilterDataLogger(EventFilter2D filter, String headerLineComment) {
        this.filter=filter;
        this.headerLine=headerLineComment;
    }
    
    synchronized public void log(String s) {
        if(logStream!=null) {
            logStream.println(s);
            if(logStream.checkError()) log.warning("eroror logging data");
        }
    }
    
    public boolean isEnabled() {
        return logDataEnabled;
    }
    
    synchronized public void setEnabled(boolean logDataEnabled){
        this.logDataEnabled = logDataEnabled;
        if(!logDataEnabled) {
            logStream.flush();
            logStream.close();
            logStream=null;
        }else{
            try{
                logStream=new PrintStream(new BufferedOutputStream(new FileOutputStream(new File(filter.getClass().getSimpleName()+".txt"))));
                logStream.println(headerLine);
            }catch(Exception e){
                e.printStackTrace();
            }
        }
        
    }
    
}
