/*
 * DVS320ControlPanel.java
 *
 * Created on October 26, 2008, 9:37 PM
 */
package ch.unizh.ini.jaer.chip.dvs320;

import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.util.HashMap;
import java.util.Observable;
import java.util.Observer;
import java.util.logging.Logger;
import javax.swing.AbstractAction;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.border.TitledBorder;

/**
 * Control panel for cDVSTest10 diagnostic output configuration.
 * @author  tobi
 */
public class cDVSTest20ControlPanel extends javax.swing.JPanel {

    static Logger log = Logger.getLogger("cDVSTest20ControlPanel");
    class OutputSelectionAction extends AbstractAction implements Observer {

        cDVSTest20.cDVSTestBiasgen.OutputMux mux;
        int channel;
        JRadioButton button;

        OutputSelectionAction(cDVSTest20.cDVSTestBiasgen.OutputMux m, int i) {
            super(m.getChannelName(i));
            mux = m;
            channel = i;
            m.addObserver(this);
        }
        
        void setButton(JRadioButton b){
            button=b;
        }

        public void actionPerformed(ActionEvent e) {
             mux.selectWithoutNotify(channel);
           log.info("Selected " + mux);
        }

        public void update(Observable o, Object arg) {
            if(channel==mux.selectedChannel) button.setSelected(true);
        }
    }
    
    cDVSTest20 chip;
    boolean panelBuilt = false;

    /** Creates new form DVS320ControlPanel */
    public cDVSTest20ControlPanel(cDVSTest20 chip) {
        this.chip = chip;
        initComponents();
    }

    void buildPanel() {
        if (panelBuilt) {
            return;
        }
        panelBuilt = true;
        cDVSTest20.cDVSTestBiasgen biasgen = (cDVSTest20.cDVSTestBiasgen) chip.getBiasgen();
        cDVSTest20.cDVSTestBiasgen.AllMuxes muxes = biasgen.allMuxes;
        for (cDVSTest20.cDVSTestBiasgen.OutputMux m : muxes) {
            JPanel p = new JPanel();
            p.setAlignmentY(0);
            p.setBorder(new TitledBorder(m.getName()));
            p.setLayout(new BoxLayout(p, BoxLayout.Y_AXIS));
            ButtonGroup group = new ButtonGroup();
            final Insets insets=new Insets(0,0,0,0);
            for (int i = 0; i < m.nInputs; i++) {
                
                JRadioButton b=new JRadioButton();
                OutputSelectionAction action=new OutputSelectionAction(m, i);
                b.setAction(action);
                action.setButton(b);
                b.setSelected(i==m.selectedChannel);
                b.setFont(b.getFont().deriveFont(10f));
                b.setToolTipText(b.getText());
                b.setMargin(insets);
//                b.setMinimumSize(new Dimension(30, 14));
                group.add(b);
                p.add(b);
            }
            add(p);
        }
        validate();
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        addComponentListener(new java.awt.event.ComponentAdapter() {
            public void componentShown(java.awt.event.ComponentEvent evt) {
                formComponentShown(evt);
            }
        });
        setLayout(new java.awt.GridLayout(1, 0));
    }// </editor-fold>//GEN-END:initComponents

private void formComponentShown(java.awt.event.ComponentEvent evt) {//GEN-FIRST:event_formComponentShown
    buildPanel();
}//GEN-LAST:event_formComponentShown
    // Variables declaration - do not modify//GEN-BEGIN:variables
    // End of variables declaration//GEN-END:variables
}
