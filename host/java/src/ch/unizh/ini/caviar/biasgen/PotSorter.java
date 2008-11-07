/*
 * PotSorter.java
 *
 * Created on November 4, 2008, 7:57 AM
 */

package ch.unizh.ini.caviar.biasgen;

import java.util.ArrayList;
import javax.swing.JComponent;

/**
 *
 * @author  thkoch
 */
public class PotSorter extends javax.swing.JPanel {
ArrayList<JComponent> guiList;
ArrayList<Pot> pots;

    /** Creates new form PotSorter */
    public PotSorter(ArrayList<JComponent> guiList,ArrayList<Pot> pots) {
        initComponents();
        this.pots=pots;
        this.guiList=guiList;
    }

    private void filterBy(String s) {
        int i=0;
        for(Pot p:pots){
            String n=p.getName();
            n=n.toLowerCase();
            String t=p.getTooltipString().toLowerCase();
            if(s==null|| s.length()==0||n.contains(s)||(t!=null&&t.contains(s))){
                guiList.get(i).setVisible(true);
            }else{
                guiList.get(i).setVisible(false);
            }
            i++;
        }
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        filterPanel = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        filterTextField = new javax.swing.JTextField();
        headerPanel = new javax.swing.JPanel();
        nameLabel = new javax.swing.JLabel();
        jPanel4 = new javax.swing.JPanel();
        sexLabel = new javax.swing.JLabel();
        jPanel1 = new javax.swing.JPanel();
        typeLabel = new javax.swing.JLabel();
        sliderAndValuePanel = new javax.swing.JPanel();
        bitValueTextField = new javax.swing.JTextField();
        bitPatternTextField = new javax.swing.JTextField();
        jPanel2 = new javax.swing.JPanel();

        setMaximumSize(new java.awt.Dimension(2147483647, 50));
        setMinimumSize(new java.awt.Dimension(151, 50));
        setPreferredSize(new java.awt.Dimension(250, 50));
        setLayout(new javax.swing.BoxLayout(this, javax.swing.BoxLayout.Y_AXIS));

        filterPanel.setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT));

        jLabel1.setLabelFor(filterTextField);
        jLabel1.setText("Filter");
        filterPanel.add(jLabel1);

        filterTextField.setColumns(20);
        filterTextField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                filterTextFieldActionPerformed(evt);
            }
        });
        filterTextField.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                filterTextFieldKeyReleased(evt);
            }
            public void keyTyped(java.awt.event.KeyEvent evt) {
                filterTextFieldKeyTyped(evt);
            }
        });
        filterPanel.add(filterTextField);

        add(filterPanel);

        headerPanel.setBorder(javax.swing.BorderFactory.createEtchedBorder());
        headerPanel.setMinimumSize(new java.awt.Dimension(109, 25));
        headerPanel.setPreferredSize(new java.awt.Dimension(254, 25));
        headerPanel.setLayout(new javax.swing.BoxLayout(headerPanel, javax.swing.BoxLayout.X_AXIS));

        nameLabel.setFont(new java.awt.Font("Microsoft Sans Serif", 1, 12));
        nameLabel.setHorizontalAlignment(javax.swing.SwingConstants.LEFT);
        nameLabel.setText("name");
        nameLabel.setHorizontalTextPosition(javax.swing.SwingConstants.RIGHT);
        nameLabel.setMaximumSize(new java.awt.Dimension(100, 15));
        nameLabel.setMinimumSize(new java.awt.Dimension(17, 10));
        nameLabel.setPreferredSize(new java.awt.Dimension(85, 15));
        headerPanel.add(nameLabel);

        jPanel4.setPreferredSize(new java.awt.Dimension(3, 0));
        headerPanel.add(jPanel4);

        sexLabel.setText("sex");
        sexLabel.setToolTipText("Sex (N- or P-type)");
        sexLabel.setMinimumSize(new java.awt.Dimension(17, 10));
        headerPanel.add(sexLabel);

        jPanel1.setPreferredSize(new java.awt.Dimension(3, 0));
        headerPanel.add(jPanel1);

        typeLabel.setText("type");
        typeLabel.setToolTipText("Type (Normal or Cascode)");
        typeLabel.setMinimumSize(new java.awt.Dimension(17, 10));
        headerPanel.add(typeLabel);

        sliderAndValuePanel.setLayout(new java.awt.BorderLayout());
        headerPanel.add(sliderAndValuePanel);

        bitValueTextField.setColumns(8);
        bitValueTextField.setEditable(false);
        bitValueTextField.setFont(new java.awt.Font("Courier New", 0, 10));
        bitValueTextField.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
        bitValueTextField.setText("bitValue");
        bitValueTextField.setToolTipText("bit value as an int");
        bitValueTextField.setMaximumSize(new java.awt.Dimension(100, 2147483647));
        bitValueTextField.setMinimumSize(new java.awt.Dimension(17, 10));
        bitValueTextField.setPreferredSize(new java.awt.Dimension(59, 10));
        bitValueTextField.addMouseWheelListener(new java.awt.event.MouseWheelListener() {
            public void mouseWheelMoved(java.awt.event.MouseWheelEvent evt) {
                bitValueTextFieldMouseWheelMoved(evt);
            }
        });
        bitValueTextField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                bitValueTextFieldActionPerformed(evt);
            }
        });
        bitValueTextField.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyPressed(java.awt.event.KeyEvent evt) {
                bitValueTextFieldKeyPressed(evt);
            }
        });
        headerPanel.add(bitValueTextField);

        bitPatternTextField.setColumns(10);
        bitPatternTextField.setEditable(false);
        bitPatternTextField.setFont(new java.awt.Font("Monospaced", 0, 10));
        bitPatternTextField.setText("bitPattern");
        bitPatternTextField.setToolTipText("bit value as bits");
        bitPatternTextField.setMaximumSize(new java.awt.Dimension(100, 2147483647));
        bitPatternTextField.setMinimumSize(new java.awt.Dimension(17, 10));
        bitPatternTextField.setPreferredSize(new java.awt.Dimension(71, 10));
        headerPanel.add(bitPatternTextField);

        jPanel2.setMaximumSize(new java.awt.Dimension(0, 32767));
        jPanel2.setMinimumSize(new java.awt.Dimension(0, 10));
        jPanel2.setPreferredSize(new java.awt.Dimension(0, 10));
        jPanel2.setRequestFocusEnabled(false);
        headerPanel.add(jPanel2);

        add(headerPanel);
    }// </editor-fold>//GEN-END:initComponents

private void bitValueTextFieldKeyPressed(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_bitValueTextFieldKeyPressed

}//GEN-LAST:event_bitValueTextFieldKeyPressed

private void bitValueTextFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_bitValueTextFieldActionPerformed

}//GEN-LAST:event_bitValueTextFieldActionPerformed

private void bitValueTextFieldMouseWheelMoved(java.awt.event.MouseWheelEvent evt) {//GEN-FIRST:event_bitValueTextFieldMouseWheelMoved

}//GEN-LAST:event_bitValueTextFieldMouseWheelMoved

private void filterTextFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_filterTextFieldActionPerformed
    String s=filterTextField.getText().toLowerCase();
    filterBy(s);
}//GEN-LAST:event_filterTextFieldActionPerformed

private void filterTextFieldKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_filterTextFieldKeyReleased
// TODO add your handling code here:
}//GEN-LAST:event_filterTextFieldKeyReleased

private void filterTextFieldKeyTyped(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_filterTextFieldKeyTyped
    String s=filterTextField.getText().toLowerCase();
    filterBy(s);
}//GEN-LAST:event_filterTextFieldKeyTyped


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JTextField bitPatternTextField;
    private javax.swing.JTextField bitValueTextField;
    private javax.swing.JPanel filterPanel;
    private javax.swing.JTextField filterTextField;
    private javax.swing.JPanel headerPanel;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel4;
    private javax.swing.JLabel nameLabel;
    private javax.swing.JLabel sexLabel;
    private javax.swing.JPanel sliderAndValuePanel;
    private javax.swing.JLabel typeLabel;
    // End of variables declaration//GEN-END:variables

}