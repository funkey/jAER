package ch.unizh.ini.jaer.projects.holger;

import gnu.io.CommPort;
import gnu.io.CommPortIdentifier;
import gnu.io.SerialPort;
//import gnu.io.SerialPortEvent;
//import gnu.io.SerialPortEventListener;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.logging.Logger;

public class PanTiltControl
{
    private Logger log = Logger.getLogger("PanTiltControl");
    public OutputStream out;
    public InputStream in;
    private boolean moving = false;
    private byte[] buffer = new byte[1024];
    private boolean connected = false;
    private boolean waitingForStarResponse = false;
    private int panPosTransformed = 0;
    int OldMinPanPos = -800;
    int OldMaxPanPos = 800;
    int NewMinPanPos = -1800;
    int NewMaxPanPos = 1800;
    boolean invert = false;

    public PanTiltControl() {
        super();
    }

    void connect ( String portName ) throws Exception
    {
        CommPortIdentifier portIdentifier = CommPortIdentifier.getPortIdentifier(portName);
        if ( portIdentifier.isCurrentlyOwned() )
        {
            log.warning("Error: Port for Pan-Tilt-Communication is currently in use");
        }
        else
        {
            CommPort commPort = portIdentifier.open(this.getClass().getName(),2000);

            if ( commPort instanceof SerialPort )
            {
                SerialPort serialPort = (SerialPort) commPort;
                serialPort.setSerialPortParams(9600,SerialPort.DATABITS_8,SerialPort.STOPBITS_1,SerialPort.PARITY_NONE);

                in = serialPort.getInputStream();
                out = serialPort.getOutputStream();

                //(new Thread(new SerialWriter(out))).start();
                //serialPort.addEventListener(new SerialReader(in));
                serialPort.notifyOnDataAvailable(true);
                connected = true;
                
                log.info("Connected to Pan-Tilt-Unit!");

            }
            else
            {
                log.warning("Error: Cannot connect to Pan-Tilt-Unit!");
            }
        }
    }

    public void setPanPos(int pos) {
        String strPanTilt = "PP" + pos + "\nA\n";
        try {
            this.out.write(strPanTilt.getBytes());
            this.moving = true;
        } catch (IOException ex) {
            log.warning("In setPanPos(position) caught IOexception " + ex);
        }
    }

    public void setPanPosTransformed(int pos) { //pos between -800 to 800
        int newpos = 0;

        if (invert==true)
            newpos=(-pos-OldMinPanPos)*(NewMaxPanPos-NewMinPanPos)/(OldMaxPanPos-OldMinPanPos)+NewMinPanPos;
        else
            newpos=(pos-OldMinPanPos)*(NewMaxPanPos-NewMinPanPos)/(OldMaxPanPos-OldMinPanPos)+NewMinPanPos;
        //log.info("newpos: " + newpos);
        panPosTransformed = pos;
        setPanPos(newpos);
    }

    public void setTransformation(int OldMinPanPos, int OldMaxPanPos, int NewMinPanPos, int NewMaxPanPos, boolean invert) {
        this.OldMinPanPos = OldMinPanPos;
        this.OldMaxPanPos = OldMaxPanPos;
        this.NewMinPanPos = NewMinPanPos;
        this.NewMaxPanPos = NewMaxPanPos;
        this.invert = invert;
    }

    public void setPanSpeed(int speed) {
        String strSpeed = "PS" + speed + "\n";
        try {
            this.out.write(strSpeed.getBytes());
            this.moving = false;
        } catch (IOException ex) {
            log.warning("In setPanSpeed caught IOexception " + ex);
        }
    }

    public void executeCommand(String command) {
        String strCommand = command + "\n";
        try {
            this.out.write(strCommand.getBytes());
        } catch (IOException ex) {
            log.warning("In executeCommand caught IOexception " + ex);
        }
    }

    public void halt() {
        String strHalt = "H\n";
        try {
            this.out.write(strHalt.getBytes());
        } catch (IOException ex) {
            log.warning("In halt() caught IOexception " + ex);
        }
    }

    public String getResponse() {
        int data;
        String response = "";
        try {
            int len = 0;
            while ((data = in.read()) > -1) {
                if (data == '\n') {
                    break;
                }
                buffer[len++] = (byte) data;
            }
            response = new String(buffer, 0, len);
            if (waitingForStarResponse == true && response.contains("*")) {
                this.moving = false;
                this.waitingForStarResponse = false;
                //log.info("Movement is done!");
            }
            else if (response.equalsIgnoreCase("A")) {
                waitingForStarResponse = true;
                //log.info("Waiting for the next star response!");
            }
            
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }
        return response;
    }

    /**
     * @return the isMoving
     */
    public boolean isMoving() {
        return moving;
    }

    /**
     * @return the isConnected
     */
    public boolean isConnected() {
        return connected;
    }

    /**
     * @return the panPosTransformed
     */
    public int getPanPosTransformed() {
        return panPosTransformed;
    }

}