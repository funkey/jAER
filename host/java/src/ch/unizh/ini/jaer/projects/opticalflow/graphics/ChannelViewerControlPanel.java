/*
 * ChannelViewerControlPanel.java
 *
 * Created on December 17, 2006, 8:30 AM
 */

package ch.unizh.ini.jaer.projects.opticalflow.graphics;

import javax.swing.*;

/**
 *
 * @author  tobi
 */
public class ChannelViewerControlPanel extends javax.swing.JPanel {
    OpticalFlowDisplayMethod displayMethod=null;
    
//    /** Creates new form ChannelViewerControlPanel */
//    public ChannelViewerControlPanel() {
//        initComponents();
//    }
    
    /** Creates new form ChannelViewerControlPanel
     @param displayMethod the OpticalFlowDisplayMethod to control
     */
    public ChannelViewerControlPanel(OpticalFlowDisplayMethod displayMethod) {
        this.displayMethod=displayMethod;
        initComponents();
        buttonGroup1.add(jRadioButton1);
        buttonGroup1.add(jRadioButton2);
        buttonGroup1.add(jRadioButton3);
        buttonGroup1.add(jRadioButton4);
        buttonGroup1.add(jRadioButton5);
        buttonGroup1.add(jRadioButton6);




        //localOffsetSlider.setValue(intFrom(displayMethod.getLocalMotionOffset()));
        //localGainSlider.setValue(intFrom(displayMethod.getLocalMotionGain()));
    }
    
    // returns int 0-100 from float 0-1
    int intFrom(float f){
        return (int)(f*100);
    }
    
    // returns a float 0-1 range value from the event assuming it is from slider that has range 0-100
    float floatFrom(javax.swing.event.ChangeEvent e){
        float f=0.01f*((JSlider)e.getSource()).getValue();
        return f;
    }
    
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        buttonGroup1 = new javax.swing.ButtonGroup();
        localPanel = new javax.swing.JPanel();
        jPanel4 = new javax.swing.JPanel();
        localOffsetSlider = new javax.swing.JSlider();
        jPanel5 = new javax.swing.JPanel();
        localGainSlider = new javax.swing.JSlider();
        displayPanel = new javax.swing.JPanel();
        jRadioButton1 = new javax.swing.JRadioButton();
        jRadioButton2 = new javax.swing.JRadioButton();
        jRadioButton3 = new javax.swing.JRadioButton();
        jRadioButton4 = new javax.swing.JRadioButton();
        jRadioButton5 = new javax.swing.JRadioButton();
        jRadioButton6 = new javax.swing.JRadioButton();

        localPanel.setBorder(javax.swing.BorderFactory.createEtchedBorder());

        jPanel4.setBorder(javax.swing.BorderFactory.createTitledBorder("Offset"));

        localOffsetSlider.setMinimumSize(new java.awt.Dimension(36, 12));
        localOffsetSlider.setPreferredSize(new java.awt.Dimension(200, 15));
        localOffsetSlider.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                localOffsetSliderStateChanged(evt);
            }
        });

        org.jdesktop.layout.GroupLayout jPanel4Layout = new org.jdesktop.layout.GroupLayout(jPanel4);
        jPanel4.setLayout(jPanel4Layout);
        jPanel4Layout.setHorizontalGroup(
            jPanel4Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(org.jdesktop.layout.GroupLayout.TRAILING, jPanel4Layout.createSequentialGroup()
                .addContainerGap()
                .add(localOffsetSlider, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 149, Short.MAX_VALUE))
        );
        jPanel4Layout.setVerticalGroup(
            jPanel4Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(jPanel4Layout.createSequentialGroup()
                .add(localOffsetSlider, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        jPanel5.setBorder(javax.swing.BorderFactory.createTitledBorder("Gain"));

        localGainSlider.setMinimumSize(new java.awt.Dimension(36, 12));
        localGainSlider.setPreferredSize(new java.awt.Dimension(200, 15));
        localGainSlider.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                localGainSliderStateChanged(evt);
            }
        });

        org.jdesktop.layout.GroupLayout jPanel5Layout = new org.jdesktop.layout.GroupLayout(jPanel5);
        jPanel5.setLayout(jPanel5Layout);
        jPanel5Layout.setHorizontalGroup(
            jPanel5Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(localGainSlider, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 149, Short.MAX_VALUE)
        );
        jPanel5Layout.setVerticalGroup(
            jPanel5Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(jPanel5Layout.createSequentialGroup()
                .add(localGainSlider, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        org.jdesktop.layout.GroupLayout localPanelLayout = new org.jdesktop.layout.GroupLayout(localPanel);
        localPanel.setLayout(localPanelLayout);
        localPanelLayout.setHorizontalGroup(
            localPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(localPanelLayout.createSequentialGroup()
                .add(jPanel4, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(jPanel5, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        localPanelLayout.setVerticalGroup(
            localPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(localPanelLayout.createSequentialGroup()
                .add(jPanel4, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .add(69, 69, 69))
            .add(localPanelLayout.createSequentialGroup()
                .add(jPanel5, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );

        displayPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Display Control"));

        jRadioButton1.setText("chan1");
        jRadioButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jRadioButton1ActionPerformed(evt);
            }
        });

        jRadioButton2.setText("chan2");

        jRadioButton3.setText("chan3");
        jRadioButton3.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jRadioButton3ActionPerformed(evt);
            }
        });

        jRadioButton4.setText("chan4");

        jRadioButton5.setText("chan6");
        jRadioButton5.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jRadioButton5ActionPerformed(evt);
            }
        });

        jRadioButton6.setText("chan5");
        jRadioButton6.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jRadioButton6ActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout displayPanelLayout = new org.jdesktop.layout.GroupLayout(displayPanel);
        displayPanel.setLayout(displayPanelLayout);
        displayPanelLayout.setHorizontalGroup(
            displayPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(displayPanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(displayPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(jRadioButton2)
                    .add(jRadioButton1)
                    .add(jRadioButton3))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, 118, Short.MAX_VALUE)
                .add(displayPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(jRadioButton5)
                    .add(jRadioButton6)
                    .add(jRadioButton4))
                .add(106, 106, 106))
        );
        displayPanelLayout.setVerticalGroup(
            displayPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(displayPanelLayout.createSequentialGroup()
                .add(jRadioButton1)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(jRadioButton2)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(jRadioButton3))
            .add(displayPanelLayout.createSequentialGroup()
                .add(jRadioButton4)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(jRadioButton6)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(jRadioButton5))
        );

        org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(layout.createSequentialGroup()
                .addContainerGap()
                .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING, false)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, localPanel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, displayPanel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap(10, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(layout.createSequentialGroup()
                .add(localPanel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 61, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(displayPanel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );
    }// </editor-fold>//GEN-END:initComponents
                        
                        
    private void localGainSliderStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_localGainSliderStateChanged
        displayMethod.setLocalMotionGain(floatFrom(evt));
    }//GEN-LAST:event_localGainSliderStateChanged
    
    
    
    private void localOffsetSliderStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_localOffsetSliderStateChanged
        displayMethod.setLocalMotionOffset(floatFrom(evt));
    }//GEN-LAST:event_localOffsetSliderStateChanged

    private void jRadioButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jRadioButton1ActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_jRadioButton1ActionPerformed

    private void jRadioButton3ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jRadioButton3ActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_jRadioButton3ActionPerformed

    private void jRadioButton5ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jRadioButton5ActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_jRadioButton5ActionPerformed

    private void jRadioButton6ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jRadioButton6ActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_jRadioButton6ActionPerformed
    
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.ButtonGroup buttonGroup1;
    private javax.swing.JPanel displayPanel;
    private javax.swing.JPanel jPanel4;
    private javax.swing.JPanel jPanel5;
    private javax.swing.JRadioButton jRadioButton1;
    private javax.swing.JRadioButton jRadioButton2;
    private javax.swing.JRadioButton jRadioButton3;
    private javax.swing.JRadioButton jRadioButton4;
    private javax.swing.JRadioButton jRadioButton5;
    private javax.swing.JRadioButton jRadioButton6;
    private javax.swing.JSlider localGainSlider;
    private javax.swing.JSlider localOffsetSlider;
    private javax.swing.JPanel localPanel;
    // End of variables declaration//GEN-END:variables
    
}
