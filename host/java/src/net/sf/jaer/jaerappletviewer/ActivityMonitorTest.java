/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * ActivityMonitorTest.java
 *
 * Created on Jan 29, 2009, 6:58:28 PM
 */
package net.sf.jaer.jaerappletviewer;

import java.awt.BorderLayout;
import javax.swing.JApplet;
import javax.swing.JFrame;

/**
 * Test DVSActApplet locally in a JFrame.
 * @author tobi
 */
public class ActivityMonitorTest extends javax.swing.JFrame {

    JApplet applet;

    /** Creates new form ActivityMonitorTest */
    public ActivityMonitorTest(JApplet applet) {
        initComponents();
        this.applet = applet;
        setSize(800,300);
        getContentPane().setLayout(new BorderLayout());
        getContentPane().add(applet,BorderLayout.CENTER);
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 400, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 300, Short.MAX_VALUE)
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        JApplet applet = new DVSActApplet();
        JFrame frame = new ActivityMonitorTest(applet);
        applet.init();
        applet.start();
        frame.setVisible(true);
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    // End of variables declaration//GEN-END:variables
}