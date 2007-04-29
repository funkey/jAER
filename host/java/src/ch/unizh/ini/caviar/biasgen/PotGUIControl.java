/*
 * IPotSliderTextControl.java
 *
 * Created on September 21, 2005, 12:23 PM
 */

package ch.unizh.ini.caviar.biasgen;

import ch.unizh.ini.caviar.util.*;
import java.awt.Toolkit;
import java.awt.event.*;
import java.lang.reflect.*;
import java.util.*;
import java.util.prefs.*;
import javax.swing.*;
import javax.swing.JSlider;
import javax.swing.border.*;
import javax.swing.border.EtchedBorder;
import javax.swing.undo.*;

/**
 * A GUI control component for controlling a Pot.
 * It shows the name of the Pot, its attributes and provides fields for direct bit editing of the Pot value. Subclasses provide customized control
 of voltage or current biases via the sliderAndValuePanel contents.
 * @author  tobi
 */
public class PotGUIControl extends javax.swing.JPanel implements  Observer, StateEditable {
    // the IPot is the master; it is an Observable that notifies Observers when its value changes.
    // thus if the slider changes the pot value, the pot calls us back here to update the appearance of the slider and of the
    // text field. likewise, if code changes the pot, the appearance here will automagically be updated.
    
    static Preferences prefs=Preferences.userNodeForPackage(IPotSliderTextControl.class);
    
    Pot pot;
    StateEdit edit=null;
    UndoableEditSupport editSupport=new UndoableEditSupport();
    BiasgenFrame frame;

    public static boolean sliderEnabled=prefs.getBoolean("PotGUIControl.sliderEnabled",true);
    public static boolean valueEnabled=prefs.getBoolean("PotGUIControl.valueEnabled",true);
    public static boolean bitValueEnabled=prefs.getBoolean("PotGUIControl.bitValueEnabled",false);
    public static boolean bitViewEnabled=prefs.getBoolean("PotGUIControl.bitViewEnabled",false);
    public static boolean sexEnabled=prefs.getBoolean("PotGUIControl.sexEnabled",true);
    public static boolean typeEnabled=prefs.getBoolean("PotGUIControl.typeEnabled",true);
    
    
    // see java tuturial http://java.sun.com/docs/books/tutorial/uiswing/components/slider.html
    // and http://java.sun.com/docs/books/tutorial/uiswing/components/formattedtextfield.html
    
    /**
     * Creates new form IPotSliderTextControl
     */
    public PotGUIControl(Pot pot, BiasgenFrame frame) {
        this.frame=frame;
        this.pot=pot;
        initComponents(); // this has unfortunate byproduect of resetting pot value to 0... don't know how to prevent stateChanged event
        if(pot!=null){
            nameLabel.setText(pot.getName()); // the name of the bias
            nameLabel.setHorizontalAlignment(SwingConstants.LEFT);
            nameLabel.setBorder(null);
            if(pot.getTooltipString()!=null) nameLabel.setToolTipText(pot.getTooltipString());
            
            typeLabel.setText(pot.getType().toString());
            sexLabel.setText(pot.getSex().toString());
            bitPatternTextField.setColumns(pot.getNumBits()+1);
            
            sliderAndValuePanel.setVisible(true);
            pot.loadPreferedBitValue(); // to get around slider value change
            pot.addObserver(this); // when pot changes, so does this gui control view
        }
        updateAppearance();  // set controls up with values from ipot
        editSupport.addUndoableEditListener(frame);
        allInstances.add(this);
    }
    
    public String toString(){
        return "IPotGUIControl for pot "+pot.getName();
    }
    
    void rr(){
        revalidate();
        repaint();
    }
    
// updates the gui slider and text fields to match actual pot values
// neither of these trigger events
    protected void updateAppearance(){
        if(pot==null) return;
        if(typeLabel.isVisible()!=typeEnabled){ typeLabel.setVisible(typeEnabled); rr(); }
        if(sexLabel.isVisible()!=sexEnabled){ sexLabel.setVisible(sexEnabled); rr(); }
        if(bitValueTextField.isVisible()!=bitValueEnabled){ bitValueTextField.setVisible(bitValueEnabled); rr(); }
        if(bitPatternTextField.isVisible()!=bitViewEnabled){ bitPatternTextField.setVisible(bitViewEnabled); rr(); }
        
        bitPatternTextField.setText(pot.toBitPatternString());
        bitValueTextField.setText(Integer.toString(pot.getBitValue()));
    }
    
     
    /** called when Observable changes (pot changes) */
    public void update(Observable observable, Object obj) {
        if(observable instanceof Pot){
//            slider.setValueIsAdjusting(false); // try to prevent a new event from the slider
//            System.out.println("ipotguicontrol observer update");
            updateAppearance();
        }
    }
    
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc=" Generated Code ">//GEN-BEGIN:initComponents
    private void initComponents() {
        nameLabel = new javax.swing.JLabel();
        jPanel4 = new javax.swing.JPanel();
        sexLabel = new javax.swing.JLabel();
        jPanel1 = new javax.swing.JPanel();
        typeLabel = new javax.swing.JLabel();
        spacerPanel = new javax.swing.JPanel();
        sliderAndValuePanel = new javax.swing.JPanel();
        jPanel3 = new javax.swing.JPanel();
        bitValueTextField = new javax.swing.JTextField();
        bitPatternTextField = new javax.swing.JTextField();
        jPanel2 = new javax.swing.JPanel();

        setLayout(new javax.swing.BoxLayout(this, javax.swing.BoxLayout.X_AXIS));

        addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseEntered(java.awt.event.MouseEvent evt) {
                formMouseEntered(evt);
            }
            public void mouseExited(java.awt.event.MouseEvent evt) {
                formMouseExited(evt);
            }
        });

        nameLabel.setFont(new java.awt.Font("Microsoft Sans Serif", 1, 12));
        nameLabel.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);
        nameLabel.setText("name");
        nameLabel.setHorizontalTextPosition(javax.swing.SwingConstants.RIGHT);
        nameLabel.setMaximumSize(new java.awt.Dimension(75, 15));
        nameLabel.setMinimumSize(new java.awt.Dimension(50, 15));
        nameLabel.setPreferredSize(new java.awt.Dimension(70, 15));
        add(nameLabel);

        jPanel4.setPreferredSize(new java.awt.Dimension(3, 0));
        add(jPanel4);

        sexLabel.setText("sex");
        sexLabel.setToolTipText("Sex (N- or P-type)");
        add(sexLabel);

        jPanel1.setPreferredSize(new java.awt.Dimension(3, 0));
        add(jPanel1);

        typeLabel.setText("type");
        typeLabel.setToolTipText("Type (Normal or Cascode)");
        add(typeLabel);

        spacerPanel.setPreferredSize(new java.awt.Dimension(3, 0));
        add(spacerPanel);

        add(sliderAndValuePanel);

        jPanel3.setMinimumSize(new java.awt.Dimension(0, 10));
        jPanel3.setPreferredSize(new java.awt.Dimension(0, 10));
        add(jPanel3);

        bitValueTextField.setColumns(8);
        bitValueTextField.setFont(new java.awt.Font("Courier New", 0, 10));
        bitValueTextField.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
        bitValueTextField.setText("bitValue");
        bitValueTextField.setToolTipText("bit value as an int");
        bitValueTextField.setMaximumSize(new java.awt.Dimension(100, 2147483647));
        bitValueTextField.setMinimumSize(new java.awt.Dimension(11, 15));
        bitValueTextField.setPreferredSize(new java.awt.Dimension(59, 15));
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
        bitValueTextField.addMouseWheelListener(new java.awt.event.MouseWheelListener() {
            public void mouseWheelMoved(java.awt.event.MouseWheelEvent evt) {
                bitValueTextFieldMouseWheelMoved(evt);
            }
        });

        add(bitValueTextField);

        bitPatternTextField.setColumns(10);
        bitPatternTextField.setEditable(false);
        bitPatternTextField.setFont(new java.awt.Font("Monospaced", 0, 10));
        bitPatternTextField.setText("bitPattern");
        bitPatternTextField.setToolTipText("bit value as bits");
        bitPatternTextField.setMaximumSize(new java.awt.Dimension(100, 2147483647));
        bitPatternTextField.setMinimumSize(new java.awt.Dimension(11, 15));
        bitPatternTextField.setPreferredSize(new java.awt.Dimension(71, 15));
        add(bitPatternTextField);

        jPanel2.setMaximumSize(new java.awt.Dimension(0, 32767));
        jPanel2.setMinimumSize(new java.awt.Dimension(0, 20));
        jPanel2.setPreferredSize(new java.awt.Dimension(0, 20));
        jPanel2.setRequestFocusEnabled(false);
        add(jPanel2);

    }// </editor-fold>//GEN-END:initComponents
    Border selectedBorder=new EtchedBorder(), unselectedBorder=new EmptyBorder(1,1,1,1);
    
    private void formMouseExited(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_formMouseExited
//        setBorder(unselectedBorder); // TODO add your handling code here:
    }//GEN-LAST:event_formMouseExited
    
    
    private void formMouseEntered(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_formMouseEntered
//        setBorder(selectedBorder);
    }//GEN-LAST:event_formMouseEntered
    
    private void bitValueTextFieldMouseWheelMoved(java.awt.event.MouseWheelEvent evt) {//GEN-FIRST:event_bitValueTextFieldMouseWheelMoved
        int clicks=evt.getWheelRotation();
        startEdit();
        pot.setBitValue(pot.getBitValue()-clicks);
        endEdit();
    }//GEN-LAST:event_bitValueTextFieldMouseWheelMoved
    
    private void bitValueTextFieldKeyPressed(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_bitValueTextFieldKeyPressed
        // key pressed in text field
//        System.out.println("bit value text field key pressed event");
        int code=evt.getKeyCode();
        boolean shift=evt.isShiftDown();
        if(!shift){
            if(code==KeyEvent.VK_UP){
                startEdit();
                pot.incrementBitValue(); // appearance updated by observer event
                endEdit();
            }else if(code==KeyEvent.VK_DOWN){
                startEdit();
                pot.decrementBitValue();
                endEdit();
            }
        }else{ // shifted, do bit shift
            int v=pot.getBitValue();
            if(code==KeyEvent.VK_UP){
                v=v<<1;
            }else if(code==KeyEvent.VK_DOWN){
                v=v>>>1;
            }
            startEdit();
            pot.setBitValue(v);
            endEdit();
        }
    }//GEN-LAST:event_bitValueTextFieldKeyPressed
    
    private void bitValueTextFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_bitValueTextFieldActionPerformed
        try{
            int v=Integer.parseInt(bitValueTextField.getText());
            startEdit();
            pot.setBitValue(v);
            endEdit();
        }catch(NumberFormatException e){
            Toolkit.getDefaultToolkit().beep();
            bitValueTextField.selectAll();
        }
        
    }//GEN-LAST:event_bitValueTextFieldActionPerformed
    
    
     private int oldPotValue=0;
     
    /** when slider is moved, event is sent here. The slider is the 'master' of the value in the text field.
     * Slider is log scale, from pot min to pot max with caveat that zero position is zero current (no current splitter
     * outputs switched on) and rest of values are log scale from pot.getCurrentResolution to pot.getMaxCurrent
     * @param e the ChangeEvent
     */   
     void startEdit(){
//        System.out.println("ipot start edit "+pot);
         edit=new MyStateEdit(this, "pot change");
         oldPotValue=pot.getBitValue();
     }
     
     void endEdit(){
         if(oldPotValue==pot.getBitValue()){
//            System.out.println("no edit, because no change in "+pot);
             return;
         }
//        System.out.println("ipot endEdit "+pot);
         if(edit!=null) edit.end();
//        System.out.println("ipot "+pot+" postEdit");
         editSupport.postEdit(edit);
     }
     
     String STATE_KEY="pot state";
     
     public void restoreState(Hashtable<?,?> hashtable) {
//        System.out.println("restore state");
         if(hashtable==null) throw new RuntimeException("null hashtable");
         if(hashtable.get(STATE_KEY)==null) {
             System.err.println("pot "+pot+" not in hashtable "+hashtable+" with size="+hashtable.size());
//            Set s=hashtable.entrySet();
//            System.out.println("hashtable entries");
//            for(Iterator i=s.iterator();i.hasNext();){
//                Map.Entry me=(Map.Entry)i.next();
//                System.out.println(me);
//            }
             return;
         }
         int v=(Integer)hashtable.get(STATE_KEY);
         pot.setBitValue(v);
     }
     
     public void storeState(Hashtable<Object, Object> hashtable) {
//        System.out.println(" storeState "+pot);
         hashtable.put(STATE_KEY, new Integer(pot.getBitValue()));
     }
     
     class MyStateEdit extends StateEdit{
         public MyStateEdit(StateEditable o, String s){
             super(o,s);
         }
         protected void removeRedundantState(){}; // override this to actually get a state stored!!
     }
     
     
     
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JTextField bitPatternTextField;
    private javax.swing.JTextField bitValueTextField;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JPanel jPanel4;
    private javax.swing.JLabel nameLabel;
    private javax.swing.JLabel sexLabel;
    private javax.swing.JPanel sliderAndValuePanel;
    private javax.swing.JPanel spacerPanel;
    private javax.swing.JLabel typeLabel;
    // End of variables declaration//GEN-END:variables
    
    
    public JTextField getBitPatternTextField() {
        return this.bitPatternTextField;
    }
    
    public JTextField getBitValueTextField() {
        return this.bitValueTextField;
    }
    
    public static boolean isBitValueEnabled() {
        return PotGUIControl.bitValueEnabled;
    }
    
    public static void setBitValueEnabled(final boolean bitValueEnabled) {
        PotGUIControl.bitValueEnabled = bitValueEnabled;
        prefs.putBoolean("PotGUIControl.bitValueEnabled", bitValueEnabled);
    }
    
    public static boolean isBitViewEnabled() {
        return PotGUIControl.bitViewEnabled;
    }
    
    public static void setBitViewEnabled(final boolean bitViewEnabled) {
        PotGUIControl.bitViewEnabled = bitViewEnabled;
        prefs.putBoolean("PotGUIControl.bitViewEnabled", bitViewEnabled);
    }
    
    public static boolean isValueEnabled() {
        return PotGUIControl.valueEnabled;
    }
    
    public static void setValueEnabled(final boolean valueEnabled) {
        PotGUIControl.valueEnabled = valueEnabled;
        prefs.putBoolean("PotGUIControl.valueEnabled", valueEnabled);
    }
    
    public static boolean isSexEnabled() {
        return PotGUIControl.sexEnabled;
    }
    
    public static void setSexEnabled(final boolean sexEnabled) {
        PotGUIControl.sexEnabled = sexEnabled;
        prefs.putBoolean("PotGUIControl.sliderEnabled", sliderEnabled);
    }
    
    public static boolean isSliderEnabled() {
        return IPotSliderTextControl.sliderEnabled;
    }
    
    public static void setSliderEnabled(final boolean sliderEnabled) {
        PotGUIControl.sliderEnabled = sliderEnabled;
        prefs.putBoolean("PotGUIControl.sliderEnabled", sliderEnabled);
    }
    
    public static boolean isTypeEnabled() {
        return PotGUIControl.typeEnabled;
    }
    
    public static void setTypeEnabled(final boolean typeEnabled) {
        PotGUIControl.typeEnabled = typeEnabled;
        prefs.putBoolean("PotGUIControl.typeEnabled", typeEnabled);
    }
    
    static ArrayList<PotGUIControl> allInstances=new ArrayList<PotGUIControl>();
    
    public static void revalidateAllInstances(){
        for(PotGUIControl c:allInstances){
            c.updateAppearance();
            c.revalidate();
        }
    }
    
    
    static String[] controlNames={"Type","Sex","Slider","BitValue","BitView"};
    public static JMenu viewMenu;
    static {
        viewMenu=new JMenu("View options");
        viewMenu.setMnemonic('V');
        for(int i=0;i<controlNames.length;i++){
            viewMenu.add(new VisibleSetter(controlNames[i])); // add a menu item to enable view of this class of information
        }
    }
    
    /** this inner static class updates the appearance of all instances of the control 
     */
    static class VisibleSetter extends JCheckBoxMenuItem{
        public String myName;
        Method setMethod,isSetMethod;
        public VisibleSetter(String myName){
            super(myName);
            this.myName=myName;
            try{
                setMethod=PotGUIControl.class.getMethod("set"+myName+"Enabled", Boolean.TYPE);
                isSetMethod=PotGUIControl.class.getMethod("is"+myName+"Enabled");
                boolean isSel=(Boolean)isSetMethod.invoke(PotGUIControl.class);
                setSelected(isSel);
            }catch(Exception e){
                e.printStackTrace();
            }
            addActionListener(new ActionListener(){
                public void actionPerformed(ActionEvent e){
                    try{
                        setMethod.invoke(IPotSliderTextControl.class, new Boolean(isSelected()));
                        setSelected(isSelected());
                        ;
                    }catch(Exception e2){
                        e2.printStackTrace();
                    }
                    PotGUIControl.revalidateAllInstances();
                }
            });
        }
    }
    
    public javax.swing.JPanel getSliderAndValuePanel() {
        return sliderAndValuePanel;
    }
    
    public void setSliderAndValuePanel(javax.swing.JPanel sliderAndValuePanel) {
        this.sliderAndValuePanel = sliderAndValuePanel;
    }
    
}

