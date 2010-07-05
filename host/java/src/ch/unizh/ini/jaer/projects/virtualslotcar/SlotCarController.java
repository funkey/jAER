/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ch.unizh.ini.jaer.projects.virtualslotcar;

/**
 * The interface for a slot car controller.
 *
 * @author tobi
 *
 * This is part of jAER
<a href="http://jaer.wiki.sourceforge.net">jaer.wiki.sourceforge.net</a>,
licensed under the LGPL (<a href="http://en.wikipedia.org/wiki/GNU_Lesser_General_Public_License">http://en.wikipedia.org/wiki/GNU_Lesser_General_Public_License</a>.
 */
public interface SlotCarController {

    /** Computes the control signal given the car tracker and the track model.
     *
     * @param tracker
     * @param track
     * @return the throttle setting ranging from 0 to 1.
     */
    public float computeControl(CarTracker tracker, SlotcarTrack track);

    /** Returns the last computed throttle setting.
     *
     * @return the throttle setting.
     */
    public float getThrottle();

}