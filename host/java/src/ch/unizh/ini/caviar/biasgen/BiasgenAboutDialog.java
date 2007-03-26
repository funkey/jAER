/*
 $Id: BiasgenAboutDialog.java,v 1.20 2004/11/12 13:37:12 tobi Exp $
 
 Copyright 2002 Institute of Neuroinformatics, University and ETH Zurich, Switzerland
 
 This file is part of The Physiologist's Friend.
 
 The Physiologist's Friend is free software; you can redistribute it
 and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2 of
 the License, or (at your option) any later version.
 
 The Physiologist's Friend is distributed in the hope that it will be
 useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with The Physiologist's Friend; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
 * Created on September 15, 2002, 8:08 PM
 
 */

package ch.unizh.ini.caviar.biasgen;


/**
 * The About dialog.  It displays About information and latest CVS commit and build dates.
 * The CVS commit date is updated by automagic ant touch of BiasgenAboutDialog.java, and the build date comes
 * from a resource TSTAMP file that is echo'ed onto the root of the classes directory by the ant build.
 * This resource is loaded and the contents are displayed.
 * 
 * @author tobi
 * @version $Revision: 1.20 $
 */
public class BiasgenAboutDialog extends javax.swing.JDialog {
    String tagName="$Name:  $";
    String revision="$Date: 2004/11/12 13:37:12 $";
    String tStampFileName="TSTAMP";
    
    /**
     * Creates new form BiasgenAboutDialog
     */
    public BiasgenAboutDialog(java.awt.Frame parent, boolean modal) {
        super(parent, modal);
        initComponents();
        
//        // when running from webstart  we are not allowed to open a file on the local file system, but we can
//        // get a the contents of a resource, which in this case is the echo'ed date stamp written by ant on the last build
//        String dateModified=null;
//        ClassLoader cl=this.getClass().getClassLoader(); // get this class'es class loader
//        System.out.println("cl="+cl);
//        URL tStampURL=cl.getResource(tStampFileName); // get a URL to the time stamp file
//        System.out.println("URL="+tStampURL);
//        try{
//            Object urlContents=tStampURL.getContent();
//            System.out.println("contents="+urlContents);
////            JOptionPane.showMessageDialog(parent,"urlContents="+urlContents);
//            BufferedReader in=null;
//            if(urlContents instanceof InputStream){
//                in=new BufferedReader(new InputStreamReader((InputStream)urlContents));
////            }else if(urlContents instanceof ZipFile){
////                ZipFile zf=(ZipFile)urlContents;
//////                JOptionPane.showMessageDialog(parent,"zf="+zf);
//////                JOptionPane.showMessageDialog(parent,zf.size()+" entries");
////                Enumeration en=zf.entries();
////                in=new BufferedReader(new InputStreamReader(.getInputStream()));
//            }
//            if(in!=null) dateModified=in.readLine();
////            JOptionPane.showMessageDialog(parent,"dateModifed="+dateModified);
//            }catch(Exception e){
//                e.printStackTrace();
//                JOptionPane.showMessageDialog(parent,e);
//            }
//            
//            aboutLabel.setText(aboutLabel.getText() +"<center>CVS "+revision+"<p>Last built "+dateModified+"</center>");
            aboutLabel.setText(aboutLabel.getText());
            pack();
        }
        
        /** This method is called from within the constructor to
         * initialize the form.
         * WARNING: Do NOT modify this code. The content of this method is
         * always regenerated by the Form Editor.
         */
    private void initComponents() {//GEN-BEGIN:initComponents
        aboutLabel = new javax.swing.JLabel();
        jSeparator1 = new javax.swing.JSeparator();
        jPanel1 = new javax.swing.JPanel();
        jPanel2 = new javax.swing.JPanel();
        okButton = new javax.swing.JButton();
        jPanel3 = new javax.swing.JPanel();

        getContentPane().setLayout(new java.awt.FlowLayout());

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                closeDialog(evt);
            }
        });

        aboutLabel.setText("<html> <center> <h1> Biasgen IPot Controller </h1> <p> The Institute of Neuroinformatics <p> Uni/ETH Zurich, Switzerland <p> <em> <a href=\"http://www.ini.unizh.ch/~tobi/\">http://www.ini.unizh.ch/~tobi/</a> </em> </center>");
        getContentPane().add(aboutLabel);

        getContentPane().add(jSeparator1);

        jPanel1.add(jPanel2);

        okButton.setText("OK");
        okButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                okButtonActionPerformed(evt);
            }
        });

        jPanel1.add(okButton);

        jPanel1.add(jPanel3);

        getContentPane().add(jPanel1);

        pack();
    }//GEN-END:initComponents
    
    private void okButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_okButtonActionPerformed
        dispose();
        // Add your handling code here:
    }//GEN-LAST:event_okButtonActionPerformed
    
    /** Closes the dialog */
    private void closeDialog(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_closeDialog
        setVisible(false);
        dispose();
    }//GEN-LAST:event_closeDialog
    
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel aboutLabel;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JSeparator jSeparator1;
    private javax.swing.JButton okButton;
    // End of variables declaration//GEN-END:variables
    
    }
/*
 $Log: BiasgenAboutDialog.java,v $
 Revision 1.20  2004/11/12 13:37:12  tobi
 on shipping friend23 and updating manual. migrating to subversion now.

 Revision 1.19  2004/02/09 06:35:17  tobi
 cleanup with refactor

 Revision 1.18  2003/07/08 15:20:58  tobi
 added back splash screen

 Revision 1.17  2003/07/07 02:44:20  tobi

 added time constant setter/getter to hcell.  messed with javadoc.

 Revision 1.16  2003/07/06 05:22:01  tobi
 *** empty log message ***

 Revision 1.15  2003/07/03 16:54:25  tobi
 fixed a bunch of javadoc errors.
 made IntegrateFireCell gettter/setter methods for settings timeconstants and used those in simulation setup factory to set complex cell properties better. (need to move this inside complex cell factory method)

 made lowpass and highpass filters time constants settable.

 Revision 1.14  2003/06/26 00:33:40  tobi

 added simulation properties dialog and fixed simple and complex cells so that they work.
 simple cell had incomplete RF. complex cell had time constant that was too long.
 fiddled with audio input and output

 Revision 1.13  2003/06/23 11:30:15  tobi
 greatly improved recording display speed, capability

 added full screen exclusive display

 Revision 1.12  2003/06/16 07:46:27  tobi
 fixed javadoc

 added target to build windows installer

 Revision 1.11  2003/06/16 06:32:52  tobi
 added more mystery cells for daniel

 Revision 1.10  2003/06/06 08:42:16  tobi
 built

 Revision 1.9  2003/05/11 09:46:13  tobi
 moved cvs log to bottom

 Revision 1.8  2003/05/11 09:34:12  tobi
 successfully included build date from combination of ant echo onto TSTAMP file
  and BiasgenAboutDialog use of the resource.  this works in webstart as well.

 Revision 1.7  2003/05/10 17:27:42  jgyger
 Merge from color-branch
 
 Revision 1.6.2.2  2003/05/08 17:10:08  tobi
 added authors and CVS tags to about box
 
 Revision 1.6.2.1  2003/05/08 17:03:31  tobi
 added authors and CVS tags to about box
 
 Revision 1.6  2002/10/24 12:05:49  cmarti
 add GPL header
 
 Revision 1.5  2002/10/08 14:53:23  tobi
 minor changes to make preliminary release tag
 
 Revision 1.4  2002/10/01 16:16:52  cmarti
 change package and import names to new hierarchy
 
 Revision 1.3  2002/09/25 08:41:26  tobi
 fixed javadoc so no errors generated and added package descriptions.
 
 Revision 1.2  2002/09/21 20:34:14  tobi
 minor changes to javadoc and to activate some menus (stimuli) in FriendGUI
 */
 // Built November 9 2004 1118  // Built September 11 2005 1006