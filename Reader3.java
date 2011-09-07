/*
 * Reader3
 * Written by: Scott Clark, Oregon State University, University of California: Davis
 *
 * This program reads in a *.pdb file for a left or right handed beta helix and computes various statistics about it including
 * The volume (total both approximated with a cylinder and Monte Carlo integration)
 * The volume of internal side chains
 * The hydrophobisity of the side chains that point in, as well as those that point out
 * As well as the hydrophobisity of the burried (non top, bottom, or loop residues)
 * The Ramachandran Plot of the helix (Phi-Psi angles)
 * The Ramachandran statistics (how many beta etc)
 * The Stacking information (what residues lie directly above and below eachother)
 * Possible sites of Hydrogen bonding
 *
 * It writes all of this information out to the terminal
 * And then makes a *.dat file called stats.dat with the information in it
 * along with the orientation of the individual resiudes
 * And in RAM.dat there is the plot if the phi-psi angles in the form:
 * residueNumber PHI PSI
 */

import java.io.*; // needed for input and output
import java.util.Random; // needed for Monte Carlo integration

public class Reader3
{
    public static void main(String args[]) throws
    IOException,
    FileNotFoundException
    {
    
        // ****************************************************************
        // *** THE FOLLOWING FOUR LINES ARE VERY IMPORTANT TO DEFINE!!! ***
    // ****************************************************************
    
    int START = 1; // Where the beta helix starts (residue number)
    int END = 61; // Where the beta helix ends (residue number)
    boolean isLH = true; // Is the helix left handed? (true = yes)
    
    // Define the *.pdb file here, only change what is in the quotes
    String ins = "DI_v2.pdb";
    
    // ***************************************************************
    // *** MAKE SURE TO DEFINE THE LINES ABOVE FOR EACH PROTIEN!!! ***
    // ***************************************************************
    
    // ********************************************************************
    // *** SWITCHES FOR EACH OF THE OUTPUT MODES, ANY COMBINATION WORKS ***
    // ********************************************************************
    
    boolean MCon = false; // turns on or off the Monte Carlo integration (off will run faster)
    boolean VolumeStats = true; // turns on or off the volume (packing) statistics
    boolean HydroStats = true; // turns on or off the hydrophobisity statistics
    boolean RamPlot = false; // turns on or off the Ramachandran plot
    boolean RamRot = false; // turns on or off the rotation of the Ramachandran plot
    boolean RamStats = true; // turns on or off the Ramachandran statistics
    boolean TransformPDB = false; // turns on or off a function that will reflect the coordinated of the *.pdb file and re-export them, more options in function.
    boolean StackStats = true; // turns on or off the Stacking statistics
    boolean HydroOn = false; // turns on or off the hydrogen bond finding
    boolean TotStats = true; // turns on or off the total statistic including orientation of individual residues
    
    // ******************** END OF SWITCHES ********************************
    
    // Variable declaration for input and output in main
    BufferedReader fin = new BufferedReader(new FileReader("" + ins + ""));
    BufferedReader fin2 = new BufferedReader(new FileReader("" + ins + ""));
    BufferedReader fin3 = new BufferedReader(new FileReader("" + ins + ""));
    PrintWriter fout = new PrintWriter(new FileOutputStream("stats_" + ins + ".dat"), true); // Where basic statistics are stored
    PrintWriter debug = new PrintWriter(new FileOutputStream("RAM.dat"), true); // Where the Ramachandran plot is stored
        String line = null;
        String line2 = null;
        String data = null;
        String data2 = null;
        String data3 = null;
        int proLen = 0;
        double CBpos[][] = new double[END-START][3]; // x = 0, y = 1, z = 2
        double CApos[][] = new double[END-START][3]; // x = 0, y = 1, z = 2
    double CPpos[][] = new double[END-START][3]; // etc.
    double Npos[][] = new double[END-START][3];
    double Opos[][] = new double[END-START][3];
    double NdonPos[][] = new double[END-START][3];
    double OdonPos[][] = new double[END-START][3];
    double NaccPos[][] = new double[END-START][3];
    double OaccPos[][] = new double[END-START][3];
    double NePos[][] = new double[END-START][3];
        double volume = 0.0, volumeTot = 0.0, volumeApprox = 0.0;
        double hydro = 0.0;
        int i = 0,j,ci=0;
        boolean stop = false;
        boolean isCoord = false;
        
    /*
     * INPUT
     */
       
    // Reads in lines and discards them until it gets to a line that starts with "SEQRES"   
        while( stop != true ){
            line = fin.readLine();
            //System.out.println("line: " + line);
            stop = (line.substring(0,4).compareTo("ATOM")==0);
            if(stop == true){
                    isCoord = true;
            }else{
            stop = (line.substring(0,6).compareTo("SEQRES")==0);
            ci++;
            }
        }
        //System.out.println("Out of loop 1: " + isCoord);
        
        
    if(isCoord == false){     
        proLen = Integer.parseInt(line.substring(14,17).trim());
    }else{
        stop = false;
        while( stop != true ){
            data2 = fin2.readLine();
            stop = (data2.substring(0,3).compareTo("END")==0);
            if(stop == false){
                if(data2.substring(0,3).compareTo("TER")==0){
                    stop = true;
                    //end = Integer.parseInt(line.substring(23,26).trim());
                }
            }
            if(stop != false){
                System.out.println(data3);
                proLen = Integer.parseInt(data3.substring(23,26).trim());
            }
            data3 = data2.toString();
        }
        stop = false;
    }
    //System.out.println("Defined seq[]: " + proLen);
        char seq[] = new char[proLen];
    if(isCoord == false){    
        data = line.substring(19);
      
    // Reads in the entire aminoacid sequence and converts it to the single letter form and then stores it in an array
        while( stop == true ){
            for(j = 0; j < 13; j++){
                if(i < proLen){
                    seq[i] = seqres(data.substring(0+(4*j),3+(4*j)));
                    //System.out.println(data.substring(0+(4*j),3+(4*j)) + " " + seqres(data.substring(0+(4*j),3+(4*j))));
                    i++;
                }
            }
            line = fin.readLine();
        // This makes sure we are still looking at the sequence and not some other part of the pdb file
            stop = (line.substring(0,6).compareTo("SEQRES")==0);
        // This makes sure we are only looking at the 'A' form of the protien
        if(stop == true){
            stop = (line.substring(11,12).compareTo("A")==0);
        if(stop == false){
            stop = (line.substring(13,14).compareTo(" ")==0);
        }
        }
            if(stop == true){
                data = line.substring(19);
            }
        }
    }else{
        int resp = 1;
        while( stop != true ){
            for(i = 0; i < START-1; i++){seq[i] = seqres("GLY");}
            for(i = 0; i <= ci; i++){line2 = fin3.readLine();}
            resp = Integer.parseInt(line2.substring(23,26).trim());
            //System.out.println("resp: " + resp);
            seq[resp-1] = seqres(line2.substring(17,20));
            if(resp == proLen){
                stop = true;
            }
            line2 = fin3.readLine();
        }
    }
    //for(i = 0; i < proLen; i++){System.out.println("Seq[" + i + "]: " + seq[i]);}
    //System.out.println("Seq initialized");
    //System.out.println("line: " + line);
    
            
        stop = false;
        
    // Reads through and discards data until it finds the individual atomic information
        while( stop != true ){
            stop = (line.substring(0,4).compareTo("ATOM")==0);
            if(stop != true){
                line = fin.readLine();
            }
        }
        System.out.println("line: " + line);
        stop = false;
        
    // Reads through and discards data until it finds the amino acid that we wanted to start with (START)
        while( stop !=true ){
            if( Integer.parseInt(line.substring(23,26).trim()) == (START) ){
                stop = true;
            }
        if(stop == false){
            line = fin.readLine();
        }
        }
//System.out.println("line: " + line);
        stop = false;
        
        for(i = 0; i < (END-START); i++){
        // Finds the position of the backbone N for a particular amino acid
        while( stop != true ){
            //System.out.println("line: " + line);
                data = line.substring(13,15);
            //System.out.println("data: " + data);
        data = data.trim();
        //System.out.println("data: " + data);
                if(data.compareTo("N") == 0){
                    //System.out.println("IN N: " + i);
                    //System.out.println("line: " + line);
                    stop = true;
                }
        if(stop == false){
            line = fin.readLine();
        }
            }
            Npos[i][0] = Double.parseDouble(line.substring(31,38).trim());
            Npos[i][1] = Double.parseDouble(line.substring(39,46).trim());
            Npos[i][2] = Double.parseDouble(line.substring(47,54).trim());
            
            stop = false;
        
        // finds the coords of the C alpha
            while( stop != true ){
                line = fin.readLine();
                data = line.substring(13,15);
                if(data.compareTo("CA") == 0){
                    stop = true;
                }
            }
            CApos[i][0] = Double.parseDouble(line.substring(31,38).trim());
            CApos[i][1] = Double.parseDouble(line.substring(39,46).trim());
            CApos[i][2] = Double.parseDouble(line.substring(47,54).trim());
            
            stop = false;
        
        // finds the coords for the C'
        while( stop != true ){
                line = fin.readLine();
                data = line.substring(13,15);
        data = data.trim();
                if(data.compareTo("C") == 0){
                    stop = true;
                }
            }
            CPpos[i][0] = Double.parseDouble(line.substring(31,38).trim());
            CPpos[i][1] = Double.parseDouble(line.substring(39,46).trim());
            CPpos[i][2] = Double.parseDouble(line.substring(47,54).trim());
            
            stop = false;
            
        // Finds the coords of the backbone O
        while( stop != true ){
                line = fin.readLine();
                data = line.substring(13,15);
        data = data.trim();
                if(data.compareTo("O") == 0){
                    stop = true;
                }
            }
            Opos[i][0] = Double.parseDouble(line.substring(31,38).trim());
            Opos[i][1] = Double.parseDouble(line.substring(39,46).trim());
            Opos[i][2] = Double.parseDouble(line.substring(47,54).trim());
        
        stop = false;
        
        /*
         * This finds the coords of the C Beta (and thus the orientation of the side chain
         * unless the residue is GLY at which point it finds the coords of the oxygen.
         */
            if(seq[START+i-1] != 'G'){
                while( stop != true ){
                    line = fin.readLine();
                    System.out.println(line);
                    data = line.substring(13,15);
                    if(data.compareTo("CB") == 0){
                        //System.out.println("IN CB: " + i);
                        //System.out.println("line: " + line);
                        stop = true;
                    }
                }
                CBpos[i][0] = Double.parseDouble(line.substring(31,38).trim());
                CBpos[i][1] = Double.parseDouble(line.substring(39,46).trim());
                CBpos[i][2] = Double.parseDouble(line.substring(47,54).trim());
            }else{
                CBpos[i][0] = Opos[i][0];
                CBpos[i][1] = Opos[i][1];
                CBpos[i][2] = Opos[i][2];
            }
            
            stop = false;
        
        /*
         * If the residue is a GLN or ASN then finds the coords of the
         * H bond accepting Oxygen and the
         * H bond donating Nitrogen
         */
        if(seq[START+i-1] == 'Q' || seq[START+i-1] == 'N'){
            while( stop != true ){
                    data = line.substring(13,15);
            data = data.trim();
            if(seq[START+i-1] == 'N'){
                        if(data.compareTo("CG") == 0){
                            stop = true;
                        }
            }else{
                if(data.compareTo("CD") == 0){
                            stop = true;
                        }
            }
            if(stop == false){
                line = fin.readLine();
            }
                }
        NePos[i][0] = Double.parseDouble(line.substring(31,38).trim());
                NePos[i][1] = Double.parseDouble(line.substring(39,46).trim());
                NePos[i][2] = Double.parseDouble(line.substring(47,54).trim());
        
        stop = false;
        
        while( stop != true ){
                    data = line.substring(13,14);
            data = data.trim();
                    if(data.compareTo("O") == 0){
                        stop = true;
                    }
            if(stop == false){
                line = fin.readLine();
            }
                }
        OaccPos[i][0] = Double.parseDouble(line.substring(31,38).trim());
                OaccPos[i][1] = Double.parseDouble(line.substring(39,46).trim());
                OaccPos[i][2] = Double.parseDouble(line.substring(47,54).trim());
        
        stop = false;
        
        while( stop != true ){
                    data = line.substring(13,14);
            data = data.trim();
                    if(data.compareTo("N") == 0){
                        stop = true;
                    }
            if(stop == false){
                line = fin.readLine();
            }
                }
        NdonPos[i][0] = Double.parseDouble(line.substring(31,38).trim());
                NdonPos[i][1] = Double.parseDouble(line.substring(39,46).trim());
                NdonPos[i][2] = Double.parseDouble(line.substring(47,54).trim());
        }
        
        stop = false;
        // If the residue is HIS, finds the coords of the H bond accepting Nitrogen
        if(seq[START+i-1] == 'H'){
            //System.out.println("H STOP: "+ seq[START+i-1]);
            while( stop != true ){
                //System.out.println("line: " + line);
                    data = line.substring(13,15);
            //data = data.trim();
                    if(data.compareTo("ND") == 0){
                        stop = true;
                        //System.out.println("ND STOP: "+ seq[START+i-1]);
                    }
            if(stop == false){
                line = fin.readLine();
            }
                }
        NaccPos[i][0] = Double.parseDouble(line.substring(31,38).trim());
                NaccPos[i][1] = Double.parseDouble(line.substring(39,46).trim());
                NaccPos[i][2] = Double.parseDouble(line.substring(47,54).trim());
        
        stop = false;
        
        while( stop != true ){
                    data = line.substring(13,15);
            data = data.trim();
                    if(data.compareTo("NE") == 0){
                        stop = true;
                    }
            if(stop == false){
                line = fin.readLine();
            }
                }
        NePos[i][0] = Double.parseDouble(line.substring(31,38).trim());
                NePos[i][1] = Double.parseDouble(line.substring(39,46).trim());
                NePos[i][2] = Double.parseDouble(line.substring(47,54).trim());
        }
        
        stop = false;
        
        /*
         * If the residue is SER, THR, or TYR; finds the coords of
         * the H bond donnor Oxygen
         * and the H bond acceptor Oxygen
         */
        if(seq[START+i-1] == 'S' || seq[START+i-1] == 'T' || seq[START+i-1] == 'Y'){
            while( stop != true ){ 
                    line = fin.readLine();
                    data = line.substring(13,14);
                    if(data.compareTo("O") == 0){
                        stop = true;
                    }
            if(stop == false){
                line = fin.readLine();
            }
                }
                OdonPos[i][0] = Double.parseDouble(line.substring(31,38).trim());
                OdonPos[i][1] = Double.parseDouble(line.substring(39,46).trim());
                OdonPos[i][2] = Double.parseDouble(line.substring(47,54).trim());
        OaccPos[i][0] = Double.parseDouble(line.substring(31,38).trim());
                OaccPos[i][1] = Double.parseDouble(line.substring(39,46).trim());
                OaccPos[i][2] = Double.parseDouble(line.substring(47,54).trim());
        }
        
        line = fin.readLine();
        stop = false;
            //System.out.println("CA: " + CApos[i][0] + " " + CApos[i][1] + " " + CApos[i][2]);
            //System.out.println("CB: " + CBpos[i][0] + " " + CBpos[i][1] + " " + CBpos[i][2]);
        }
    
    fin.close();
    
        /*
     * OUTPUT
     */
         
    //Calculates the volume of the helix using Monte Carlo techniques
    if(MCon){
            volumeTot = findVolume(CApos,END-START,true,isLH);
    }
        
    // VOLUME STATS
    if(VolumeStats){
        //Calculates the volume of the helix using a triangular cylinder approximation
        volumeApprox = findVolApp(CApos,START,END);
    
        //finds the volume of all side chains pointing in
            for(i = 0; i < (END-START); i++){
                if(CBisInside(i,CBpos[i][0],CBpos[i][1],CBpos[i][2],CApos,END-START,false,isLH)){
                    volume += getVol(seq[START+i-1]);
                }
            }
        //Output of data and various statistics
        System.out.println("Volume: " + volume);
            System.out.println("Volume total (of helix) (MC): " + volumeTot);
        System.out.println("Volume ratio V/Vt(MC): " + (volume/volumeTot));
        if(isLH){
            System.out.println("Volume total (of helix) (Ap): " + volumeApprox);
                System.out.println("Volume ratio V/Vt(Ap): " + (volume/volumeApprox));
        }
        fout.println("Volume: " + volume);
            fout.println("Volume total (of helix) (MC): " + volumeTot);
        fout.println("Volume ratio V/Vt(MC): " + (volume/volumeTot));
        if(isLH){
            fout.println("Volume total (of helix) (Ap): " + volumeApprox);
                fout.println("Volume ratio V/Vt: " + (volume/volumeApprox));
        }
    }
        
    // HYDROPHOBISITY
    if(HydroStats){
        int hydroStati[] = new int[5];
        int hydroStato[] = new int[5];
        int hydroOut = 0;
            for(i = 0; i < (END-START); i++){
                if(CBisInside(i,CBpos[i][0],CBpos[i][1],CBpos[i][2],CApos,END-START,false,isLH)){ //Checks whether or not the amino acid is inside or outside
                    hydro += getHydro(seq[START+i-1]); // adds to the running hydro count for inside
            hydroStati[getHydroType(seq[START+i-1])]++; //keeps track of what type of hydrophobisity it has: acid, base, hydrophobe, hydrophile
            // getHydroType: 0 = hydrophobe, 1 = hydrophile, 2 = base, 3 = acid
            if(seq[START+i-1] == 'G'){
                hydroStati[4]++; // keeps track of the number of GLY that point in
            }
                }else{
                hydroOut += getHydro(seq[START+i-1]); // adds to the running count of the outside hydrophobisity
            hydroStato[getHydroType(seq[START+i-1])]++; // keeps track of the types that are facing out
            // getHydroType: 0 = hydrophobe, 1 = hydrophile, 2 = base, 3 = acid
            if(seq[START+i-1] == 'G'){
                hydroStato[4]++; // keeps track of the number of GLY
            }
            }
            }
        // Output of data and various statistics
        System.out.println("Hydrophobisity:");
        System.out.println("IN: " + hydroStati[0] + " Hydrophobic residue(s) " + hydroStati[1] + " Hydrophilic residue(s) including " + hydroStati[4] + " GLY residue(s)");
        System.out.println("IN: " + hydroStati[2] + " basic residue(s) " + hydroStati[3] + " acidic residue(s)");
        System.out.println("IN: Total hydrophobisity: " + hydro);
        System.out.println("OUT: " + hydroStato[0] + " Hydrophobic residue(s) " + hydroStato[1] + " Hydrophilic residue(s) including " + hydroStato[4] + " GLY residue(s)");
        System.out.println("OUT: " + hydroStato[2] + " basic residue(s) " + hydroStato[3] + " acidic residue(s)");
        System.out.println("OUT: Total hydrophobisity: " + hydroOut);
        fout.println("Hydrophobisity:");
        fout.println("IN: " + hydroStati[0] + " Hydrophobic residue(s) " + hydroStati[1] + " Hydrophilic residue(s) including " + hydroStati[4] + " GLY residue(s)");
        fout.println("IN: " + hydroStati[2] + " basic residue(s) " + hydroStati[3] + " acidic residue(s)");
        fout.println("IN: Total hydrophobisity: " + hydro);
        fout.println("OUT: " + hydroStato[0] + " Hydrophobic residue(s) " + hydroStato[1] + " Hydrophilic residue(s) including " + hydroStato[4] + " GLY residue(s)");
        fout.println("OUT: " + hydroStato[2] + " basic residue(s) " + hydroStato[3] + " acidic residue(s)");
        fout.println("OUT: Total hydrophobisity: " + hydroOut);
        
        // Burried Hydrophobisity
        // This looks at the hydrdophobisity of the residues NOT on the top, bottom or loop of the helix
        for(i = 0; i < 5; i++){
            hydroStati[i] = 0;
            hydroStato[i] = 0;
        }
        hydroOut = 0;
            for(i = 0; i < (END-START); i++){
            if(getNplus(i, 3, 0, 9, Math.min(30,(END-START) - i), CApos) !=0 && getNminus(i, 3, 0, 9, Math.min(30,i), CApos) != 0){
                    if(CBisInside(i,CBpos[i][0],CBpos[i][1],CBpos[i][2],CApos,END-START,false,isLH)){ //Checks whether or not the amino acid is inside or outside
                        hydro += getHydro(seq[START+i-1]); // adds to the running hydro count for inside
                hydroStati[getHydroType(seq[START+i-1])]++; //keeps track of what type of hydrophobisity it has: acid, base, hydrophobe, hydrophile
                // getHydroType: 0 = hydrophobe, 1 = hydrophile, 2 = base, 3 = acid
                if(seq[START+i-1] == 'G'){
                    hydroStati[4]++; // keeps track of the number of GLY that point in
                }
                    }else{
                    hydroOut += getHydro(seq[START+i-1]); // adds to the running count of the outside hydrophobisity
                hydroStato[getHydroType(seq[START+i-1])]++; // keeps track of the types that are facing out
                // getHydroType: 0 = hydrophobe, 1 = hydrophile, 2 = base, 3 = acid
                if(seq[START+i-1] == 'G'){
                    hydroStato[4]++; // keeps track of the number of GLY
                }
                }
        }
            }
        // Output of data and various statistics
        System.out.println("Burried Hydrophobisity:");
        System.out.println("IN: " + hydroStati[0] + " Hydrophobic residue(s) " + hydroStati[1] + " Hydrophilic residue(s) including " + hydroStati[4] + " GLY residue(s)");
        System.out.println("IN: " + hydroStati[2] + " basic residue(s) " + hydroStati[3] + " acidic residue(s)");
        System.out.println("IN: Total hydrophobisity: " + hydro);
        System.out.println("OUT: " + hydroStato[0] + " Hydrophobic residue(s) " + hydroStato[1] + " Hydrophilic residue(s) including " + hydroStato[4] + " GLY residue(s)");
        System.out.println("OUT: " + hydroStato[2] + " basic residue(s) " + hydroStato[3] + " acidic residue(s)");
        System.out.println("OUT: Total hydrophobisity: " + hydroOut);
        fout.println("Burried Hydrophobisity:");
        fout.println("IN: " + hydroStati[0] + " Hydrophobic residue(s) " + hydroStati[1] + " Hydrophilic residue(s) including " + hydroStati[4] + " GLY residue(s)");
        fout.println("IN: " + hydroStati[2] + " basic residue(s) " + hydroStati[3] + " acidic residue(s)");
        fout.println("IN: Total hydrophobisity: " + hydro);
        fout.println("OUT: " + hydroStato[0] + " Hydrophobic residue(s) " + hydroStato[1] + " Hydrophilic residue(s) including " + hydroStato[4] + " GLY residue(s)");
        fout.println("OUT: " + hydroStato[2] + " basic residue(s) " + hydroStato[3] + " acidic residue(s)");
        fout.println("OUT: Total hydrophobisity: " + hydroOut);
        
    }
    
    // Ramachandran plot
    // Plots the phi vs. psi angles of each residue on the helix
    if(RamPlot){
        boolean includeGLY = false;
        for(i = 0; i < (END-START); i++){
            if(!includeGLY){
            if(seq[START+i-1] != 'G'){
                    debug.println(i + " " + getPhi(i, CApos, Npos, CPpos, END-START) + " " + getPsi(i, CApos, Npos, CPpos, END-START));
                }
        }else{
            debug.println(i + " " + getPhi(i, CApos, Npos, CPpos, END-START) + " " + getPsi(i, CApos, Npos, CPpos, END-START));
        }
        }
    }
    
    if(RamRot){
        double theta = 3.14159/4;
            PrintWriter foutRR = new PrintWriter(new FileOutputStream("rotRAM.dat"), true);
        double phi, psi, temp;
        int space1, space2;
        i = 0;
        stop = false;
        line = null;
        data = null;
        int sector;
        double vphi[] = new double[2];
        double vpsi[] = new double[2];
        double vrot[] = new double[2];
        double nvrot[] = new double[2];
        double pvrot[] = new double[2];
        double npvrot[] = new double[2];
        double beta = 0.0;
        double dot, magPhiPsi;
    
        vrot[0] = Math.cos(theta);
        vrot[1] = -1.0*Math.sin(theta);
    
        nvrot[0] = -1.0*vrot[0];http://www.pricegrabber.com/
        nvrot[1] = -1.0*vrot[1];
    
        pvrot[0] = Math.sin(theta);
        pvrot[1] = Math.cos(theta);
    
        npvrot[0] = -1.0*pvrot[0];
        npvrot[1] = -1.0*pvrot[1];
    
        while(stop != true){
            line = fin.readLine();
            stop = (line.substring(0,3).compareTo("END")==0);
            if(stop != true){
                space1 = line.indexOf(' ');
            data = line.substring(space1+1);
            space2 = data.indexOf(' ');
            space2 = space1 + space2 + 2;
                phi = Double.parseDouble(line.substring(space1+1,space2).trim());
            psi = Double.parseDouble(line.substring(space2).trim());
        
            magPhiPsi = Math.sqrt(phi*phi + psi*psi);
        
            sector = getSector(phi, psi, vrot[0], vrot[1], pvrot[0], pvrot[1]);
        
            switch(sector){
                case 1:
                    dot = phi*nvrot[0] + psi*nvrot[1];
                beta = Math.abs(Math.acos(dot/magPhiPsi));
                beta = 2.0*beta;
                temp = phi;
                phi = Math.cos(beta)*phi + Math.sin(beta)*psi;
                psi = -1.0*Math.sin(beta)*temp + Math.cos(beta)*psi;
                break;
                case 2:
                    dot = phi*nvrot[0] + psi*nvrot[1];
                beta = Math.abs(Math.acos(dot/magPhiPsi));
                beta = 2.0*beta;
                temp = phi;
                phi = Math.cos(beta)*phi - Math.sin(beta)*psi;
                psi = Math.sin(beta)*temp + Math.cos(beta)*psi;
                break;
                case 3:
                    dot = phi*vrot[0] + psi*vrot[1];
                beta = Math.abs(Math.acos(dot/magPhiPsi));
                beta = 2.0*beta;
                temp = phi;
                phi = Math.cos(beta)*phi + Math.sin(beta)*psi;
                psi = -1.0*Math.sin(beta)*temp + Math.cos(beta)*psi;
                break;
                case 4:
                    dot = phi*vrot[0] + psi*vrot[1];
                beta = Math.abs(Math.acos(dot/magPhiPsi));
                beta = 2.0*beta;
                temp = phi;
                phi = Math.cos(beta)*phi - Math.sin(beta)*psi;
                psi = Math.sin(beta)*temp + Math.cos(beta)*psi;
                break;
                default: break;
            }
            foutRR.println("residue["+i+"].setDihedrals(" + phi + "," + psi + ");");
                i++;
            }
        }
    }
    
    // Ramachandran Statistics
    if(RamStats){
        int resType[] = new int[END-START];
        int resTypeT[] = new int[END-START];
        boolean counted[] = new boolean[END-START];
        int typeTot[] = new int[3];
        int typeT[] = new int[3];
        int tot = 0, tottot = 0;
        int typeCount = 1;
        for(i = 0; i < (END-START); i++){
            counted[i] = false;
        }
        for(i = 0; i < (END-START); i++){
            if(i != 0 && i != (END-START-1)){
                if(onCorner(i, CApos)){
                //System.out.println("CORNER AT: " + (START+i));
                    for(j = -1; j < 1; j++){
                if(counted[i+j] == false){
                    if(seq[START+i+j-1] != 'G'){
                            resTypeT[i+j] = classifyPhiPsi(getPhi((i+j), CApos, Npos, CPpos, END-START),getPsi((i+j), CApos, Npos, CPpos, END-START));
                                typeT[resTypeT[i+j]]++;
                            tot++;
                    counted[i+j] = true;
                }
                    }
                    }
            }
            }
            if(seq[START+i-1] != 'G'){
                resType[i] = classifyPhiPsi(getPhi((i), CApos, Npos, CPpos, END-START),getPsi((i), CApos, Npos, CPpos, END-START));
                typeTot[resType[i]]++;
                tottot++;
            }
        }
        System.out.println("Ramachandran Plot stored in RAM.dat");
        System.out.println("There were " + typeTot[0] + " (" + ((double)typeTot[0]/(double)tottot) + ") Beta residues, " + typeTot[1] + " (" + ((double)typeTot[1]/(double)tottot) + ") RH alpha and " + typeTot[2] + " (" + ((double)typeTot[2]/(double)tottot) + ") LH alpha of " + tottot + " non-GLY residues.");
        System.out.println("On the turns: " + typeT[0] + " (" + ((double)typeT[0]/(double)tot) + ") were Beta, " + typeT[1] + " (" + ((double)typeT[1]/(double)tot) + ") RH alpha and " + typeT[2] + " (" + ((double)typeT[2]/(double)tot) + ") LH alpha of " + tot + " non-GLY residues");
        fout.println("There were " + typeTot[0] + " (" + ((double)typeTot[0]/(double)tottot) + ") Beta residues, " + typeTot[1] + " (" + ((double)typeTot[1]/(double)tottot) + ") RH alpha and " + typeTot[2] + " (" + ((double)typeTot[2]/(double)tottot) + ") LH alpha of " + tottot + " non-GLY residues.");
        fout.println("On the turns: " + typeT[0] + " (" + ((double)typeT[0]/(double)tot) + ") were Beta, " + typeT[1] + " (" + ((double)typeT[1]/(double)tot) + ") RH alpha and " + typeT[2] + " (" + ((double)typeT[2]/(double)tot) + ") LH alpha of " + tot + " non-GLY residues");
        for(i = 0; i < (END-START); i++){
            if(resType[i] != 0){
            for(j = 1; j < Math.min((typeTot[1] + typeTot[2]), (END-START) - i - 1); j++){
                if(resType[i+j] != 0){
                typeCount++;
            }else{
                i = i + j - 1;
                j = (typeTot[1] + typeTot[2]);
            }
            }
            if(typeCount > 1){
                System.out.println("There was a sequence of " + typeCount + " alpha pieces in a row (" + (START + i) + "-" + (START + i + typeCount - 1) + ").");
                fout.println("There was a sequence of " + typeCount + " alpha pieces in a row (" + (START + i) + "-" + (START + i + typeCount - 1) + ").");
            typeCount = 1;
            }
        }
        }
    }
    
    // Transformation of *.pbd coordinate file
    if(TransformPDB){
        boolean rotX = false;
        boolean rotY = false;
        boolean rotZ = false;
        transPDB(rotX, rotY, rotZ, ins);
    }
    
    int Np;
    //Stacking info
    if(StackStats){
        stop = false;
        int temp = -1;
        int stack[] = new int[20];
        int dStack[][] = new int[25][3];
        int dStackc = 0;
        int pStack[][] = new int[50][2];
        int numStack = 0;
        // loops through entire helix
        for(i = 0; i < (END-START); i++){
            Np = getNplus(i, 3, 0, 9, Math.min(30,(END-START) - i), CApos); // finds how far away the residue directly above it in 3D is
            if(stop == false){
                if(i < (END-START-1)){
                    temp = getNplus(i+1, 3, 0, 9, Math.min(30,(END-START) - i), CApos); // finds how far away the residue above the residue next to it is
                }
                if(Np == temp){ // checks to see if they are the same distance away
                    stop = true; // if they are then it alows the next algorithm to search for a possible stack
                // this is put in to avoid multiple residues thinking that the same residue is directly above them
                // which would result in stacks recorded that were not actual stacks
            }
            }
            //System.out.println((i+START)+ " " + (i+START+Np) + " " + seq[i+START-1] + " " + seq[i+Np+START-1]);
            if(Np != 0){ // makes sure that the residue will not check similarity with itself, but rather the residue above it
                if(stop == true){ // makes sure the basic error checking above was satisfied
                if(seq[i+START-1] != 'G'){ // we are ignoring GLY because it does not have enough of a side chain to make stacking relevent
                        if(seq[i+START-1] == seq[i+Np+START-1]){ // this actually checks whether or not the residues are the same (the two residues directly next to eachother vertically
                    // This next line makes sure that either both residues are pointing in or pointing out or pointing in (stacking would be irrelevent if they were on opposite sides of the helical wall
                    if(CBisInside(i-1,CBpos[i-1][0],CBpos[i-1][1],CBpos[i-1][2],CApos,(END-START),false,isLH)==CBisInside(i+Np-1,CBpos[i+Np-1][0],CBpos[i+Np-1][1],CBpos[i+Np-1][2],CApos,(END-START),false,isLH)){
                        // this records information about where the stack was for use in testing for a "double stack"
                    pStack[numStack][0] = i+Np+START;
                        pStack[numStack][1] = Np;
                    // this searches through and looks for double stacks and records where they are
                        for(j = 0; j < numStack; j++){
                            if((i+START) == pStack[j][0]){
                            dStack[dStackc][0] = i+START-pStack[j][1];
                            dStack[dStackc][1] = i+START;
                            dStack[dStackc][2] = i+START+Np;
                            dStackc++;
                        }
                            }
                        numStack++;
                    // this records where the stacks were in the dat file as well as on the screen
                    // It also records whether the stack was inside or outside
                        if((START + Np + i) != END){
                        if(CBisInside(i-1,CBpos[i][0],CBpos[i][1],CBpos[i][2],CApos,(END-START),false,isLH)){
                                fout.println("Inward Stack at " + (START + i) + " and " + (START + Np + i) + " with " + seq[START + i-1]);
                                    System.out.println("Inward Stack at " + (START + i) + " and " + (START + Np + i) + " with " + seq[START + i-1]);
                        }else{
                            fout.println("Outward Stack at " + (START + i) + " and " + (START + Np + i) + " with " + seq[START + i-1]);
                                    System.out.println("Outward Stack at " + (START + i) + " and " + (START + Np + i) + " with " + seq[START + i-1]);
                        }
                        }else{ // a * is appended for stacks involving the final residue because sometimes they can be suspect (not perfectly aligned)
                            if(CBisInside(i-1,CBpos[i-1][0],CBpos[i-1][1],CBpos[i-1][2],CApos,(END-START),false,isLH)){
                                fout.println("Inward Stack at " + (START + i) + " and " + (START + Np + i) + " with " + seq[START + i-1] + "*");
                                    System.out.println("Inward Stack at " + (START + i) + " and " + (START + Np + i) + " with " + seq[START + i-1] + "*");
                        }else{
                            fout.println("Outward Stack at " + (START + i) + " and " + (START + Np + i) + " with " + seq[START + i-1] + "*");
                                    System.out.println("Outward Stack at " + (START + i) + " and " + (START + Np + i) + " with " + seq[START + i-1] + "*");
                        }
                        }
                        stack[getNum(seq[START + i-1])]++; // This records what type of residue was stacked
                }
                }
                }
            }
            }
        }
    
        // Spits out the info on the double stacks that were found
        for(i = 0; i < dStackc; i++){
            if(CBisInside((dStack[i][0]-START),CBpos[dStack[i][0]-START][0],CBpos[dStack[i][0]-START][1],CBpos[dStack[i][0]-START][2],CApos,(END-START),false,isLH)){
                System.out.println("Inward double stack of " + (char)seq[dStack[i][0]-1] + " at " + dStack[i][0] + "-" + dStack[i][1] + "-" + dStack[i][2]);
                fout.println("Inward double stack of " + (char)seq[dStack[i][0]-1] + " at " + dStack[i][0] + "-" + dStack[i][1] + "-" + dStack[i][2]);
            }else{
                System.out.println("Outward double stack of " + (char)seq[dStack[i][0]-1] + " at " + dStack[i][0] + "-" + dStack[i][1] + "-" + dStack[i][2]);
                fout.println("Outward double stack of " + (char)seq[dStack[i][0]-1] + " at " + dStack[i][0] + "-" + dStack[i][1] + "-" + dStack[i][2]);
            }
        }
        // Spits out the info on what residues stacked and at what frequency
        for(i = 0; i < 20; i++){
            if(stack[i] != 0){
                System.out.println((char)getName(i) + " stacked " + stack[i] + " time(s)");
            fout.println((char)getName(i) + " stacked " + stack[i] + " time(s)");
            }
        }
    }
    
    // Hydrogen Bonding
    if(HydroOn){
        //PrintWriter hout = new PrintWriter(new FileOutputStream("Reader3Hbond.dat"), true);
        for(i = 0; i < (END-START); i++){
            getHbond(i,NdonPos,NaccPos,Npos,OdonPos,OaccPos,Opos,CPpos,CBpos,NePos,START,END,seq, fout);
        }
        //hout.close();
    }
    
    // Total Stats
    if(TotStats){
        // Gives information about the dimensions of the helix that was worked on
        double turns = 0;
        Np = 0;
        // loops through until it cannot find another residue directly above the current one, counting the rings as it goes
        for(i = 0; i < (END-START); i=i){
            Np = getNplus(i, 3, 0, 9, Math.min(30,(END-START) - i), CApos);
            i += Np;
            if(Np != 0){
                turns = turns + 1.0;
            }else{
                i = END-START;
                if((i+30)<END){
                        turns = 0.0;
            }
            }
        }
        // if it failed (due to complicated side loop), or if wanted, this gives a 18 residue/ring approximation for the turns. turns = 0; will force this approximation
        if(turns == 0){
            if(isLH){
                turns = ((END-START)/18.0);
        }else{
            turns = ((END-START)/22.0);
        }
        }      
            System.out.println("There were " + (END-START) + " residues (" + START + "-" + END + "), and approximatly " + turns + " turns.");
        fout.println("There were " + (END-START) + " residues (" + START + "-" + END + "), and approximatly " + turns + " turns.");
    
    
        // Sequence and oreintation (true = in, false = out)
            for(i = 0; i < (END-START) ; i++){
                fout.println(""+ (START+i)+ " " + seq[START+i-1] + ": " + CBisInside(i,CBpos[i][0],CBpos[i][1],CBpos[i][2],CApos,END-START,false,isLH));
            }
    
        System.out.println("This data and the chain information including reisdue name, number, and whether it is fliped in (true) or out (false) in stats.dat");
    }
    }
    
    /*
     * getHbond
     * takes int, double[][], double[][], double[][], double[][], double[][], double[][], double[][], int, int ,char
     * returns int (dummy placeholder for now)
     * Input: (residue to test for possible H bonds on, Nitrogen donor position array, Nitrogen acceptor position array
     *      Nitrogen backbone position array, Oxygen donor position array, Oxygen acceptor position array,
     *      Oxygen backbone position array, the C' position arrau, the C beta position array,
     *          the bond angle atom postion array (see comments below), start residue, end residue, sequence array)
     *
     * This function takes in information regarding the positions of the atoms that play a role in hydrogen bonding,
     * it determines whether or not there is a possible hydrogen bond by finding the distance between two atoms that
     * could participate in a hydrogen bond with eachother. The distances were taken from Whitford's Proteins. The distance
     * of the covalent bond between the atom and its hydrogen is assumed to be 0.96 A (Wikipedia). If there is a hydrogen bond
     * the type of H bond as well as the asociated atom and residues are outputed to the screen. This information is then also
     * outputed to a file Reader3Hbond.dat
     *
     * To determine the bond angles of the H bonds we need to have another atom take another atom into account to determine
     * the angles at which the hydrogen bonds are formed and throw out the ones at too sharp of an angle.
     * To do this we look at the next heavy atom next to the acceptor or donor and find the angle between the vector that
     * runs from heavy atom to heavy atom along the covalent bond and the vector that runs from the heavy atom all the way
     * across the H bond to the donor or accepto atom. This will define the angle of our H bond and we can put restrictions on it
     * In Protiens page 56 it is stated that H bonds can have up to a +/- 40 degree shift from the horizontal. In the 3 spatial
     * dimensions that we are working in this equates to a cone of possible H bonding for every donor/acceptor. The heavy atoms
     * that we are using to compute the H bond angles are:
     * Amide-carbonyl: C'
     * Amide-hydroxyl: C beta
     * Amide-imidozole: NE2
     * Hydroxyl-carbonyl: C'
     * Hydroxyl-Hydroxyl: no restriction
     *
     * an example of the angle is (amide-carbonyl):
     *
     *      \              /
     *       N - H -- O = C'
     *      /              \
     *
     * The angle between the vector that points from C' to O and the vector that points from C' to N.
     */
    private static int getHbond(int res, double Nd[][], double Na[][], double N[][], double Od[][], double Oa[][], double O[][], double cp[][], double cb[][], double ne[][], int start, int end, char seq[], PrintWriter hout) throws
    IOException,
    FileNotFoundException
    {
        
    double dist = 0.0;
    int END = end, START = start;
    int i;
    if(res == 0){
        System.out.println("Side chain hydrogen bonds");
        hout.println("Side chain hydrogen bonds");
    }
    switch(seq[res+start-1]){
        case 'Q':
            for(i = (res - Math.min(30,res)); i < Math.min(30,(END-START-res)); i++){
            dist = Math.sqrt((Nd[res][0] - O[res + i][0])*(Nd[res][0] - O[res + i][0]) + (Nd[res][1] - O[res + i][1])*(Nd[res][1] - O[res + i][1]) + (Nd[res][2] - O[res + i][2])*(Nd[res][2] - O[res + i][2]));
            if(dist < 3.86 && i != 0){
                if(angleWithin(Nd[res],cp[res+i],O[res+i],.698132)){
                    System.out.println("Amide - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and backbone O " + (i+start+res));
                hout.println("Amide - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and backbone O " + (i+start+res));
            }
            }
            if(seq[start + res + i - 1] == 'S'){
                dist = Math.sqrt((Nd[res][0] - Oa[res + i][0])*(Nd[res][0] - Oa[res + i][0]) + (Nd[res][1] - Oa[res + i][1])*(Nd[res][1] - Oa[res + i][1]) + (Nd[res][2] - Oa[res + i][2])*(Nd[res][2] - Oa[res + i][2]));
            if(dist < 3.96 && i != 0){
                if(angleWithin(Nd[res],cb[res+i],Oa[res+i],.698132)){
                        System.out.println("Amide - hydroxyl H bond at " + (res+start) + " " + seq[res+start-1] + " (N) to " + (start+res+i) + " (O)");
                    hout.println("Amide - hydroxyl H bond at " + (res+start) + " " + seq[res+start-1] + " (N) to " + (start+res+i) + " (O)");
                    }
            }
            }
            if(seq[start + res + i - 1] == 'H'){
                dist = Math.sqrt((Nd[res][0] - Na[res + i][0])*(Nd[res][0] - Na[res + i][0]) + (Nd[res][1] - Na[res + i][1])*(Nd[res][1] - Na[res + i][1]) + (Nd[res][2] - Na[res + i][2])*(Nd[res][2] - Na[res + i][2]));
            if(dist < 4.06 && i != 0){
                if(angleWithin(Nd[res],ne[res+i],Na[res+i],.698132)){
                        System.out.println("Amide - imidazole H bond at " + (res+start) + " " + seq[res+start-1] + " (N) to " + (start+res+i) + " (N)");
                    hout.println("Amide - imidazole H bond at " + (res+start) + " " + seq[res+start-1] + " (N) to " + (start+res+i) + " (N)");
                }
                }
            }
        }
        break;
            case 'N':
            for(i = (res - Math.min(30,res)); i < Math.min(30,(END-START-res)); i++){
            dist = Math.sqrt((Nd[res][0] - O[res + i][0])*(Nd[res][0] - O[res + i][0]) + (Nd[res][1] - O[res + i][1])*(Nd[res][1] - O[res + i][1]) + (Nd[res][2] - O[res + i][2])*(Nd[res][2] - O[res + i][2]));
            if(dist < 3.86 && i != 0){
                if(angleWithin(Nd[res],cp[res+i],O[res+i],.698132)){
                    System.out.println("Amide - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and backbone O " + (i+start+res));
                hout.println("Amide - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and backbone O " + (i+start+res));
                }
            }
            if(seq[start + res + i - 1] == 'S'){
                dist = Math.sqrt((Nd[res][0] - Oa[res + i][0])*(Nd[res][0] - Oa[res + i][0]) + (Nd[res][1] - Oa[res + i][1])*(Nd[res][1] - Oa[res + i][1]) + (Nd[res][2] - Oa[res + i][2])*(Nd[res][2] - Oa[res + i][2]));
            if(dist < 3.96 && i != 0){
                if(angleWithin(Nd[res],cb[res+i],Oa[res+i],.698132)){
                        System.out.println("Amide - hydroxyl H bond at " + (res+start) + " " + seq[res+start-1] + " (N) to " + (start+res+i) + " (O)");
                    hout.println("Amide - hydroxyl H bond at " + (res+start) + " " + seq[res+start-1] + " (N) to " + (start+res+i) + " (O)");
                    }
            }
            }
            if(seq[start + res + i - 1] == 'H'){
                dist = Math.sqrt((Nd[res][0] - Na[res + i][0])*(Nd[res][0] - Na[res + i][0]) + (Nd[res][1] - Na[res + i][1])*(Nd[res][1] - Na[res + i][1]) + (Nd[res][2] - Na[res + i][2])*(Nd[res][2] - Na[res + i][2]));
            if(dist < 4.06 && i != 0){
                if(angleWithin(Nd[res],ne[res+i],Na[res+i],.698132)){
                        System.out.println("Amide - imidazole H bond at " + (res+start) + " " + seq[res+start-1] + " (N) to " + (start+res+i) + " (N)");
                    hout.println("Amide - imidazole H bond at " + (res+start) + " " + seq[res+start-1] + " (N) to " + (start+res+i) + " (N)");
                }
                }
            }
        }
        break;
        case 'S':
            for(i = (res - Math.min(30,res)); i < Math.min(30,(END-START-res)); i++){
            dist = Math.sqrt((Od[res][0] - O[res + i][0])*(Od[res][0] - O[res + i][0]) + (Od[res][1] - O[res + i][1])*(Od[res][1] - O[res + i][1]) + (Od[res][2] - O[res + i][2])*(Od[res][2] - O[res + i][2]));
            if(dist < 3.76 && i != 0){
                if(angleWithin(Od[res],cp[res+i],Od[res+i],.698132)){
                    System.out.println("Hydroxyl - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and backbone O " + (i+start+res));
                hout.println("Hydroxyl - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and backbone O " + (i+start+res));
            }
            }
            dist = Math.sqrt((Od[res][0] - N[res + i][0])*(Od[res][0] - N[res + i][0]) + (Od[res][1] - N[res + i][1])*(Od[res][1] - N[res + i][1]) + (Od[res][2] - N[res + i][2])*(Od[res][2] - N[res + i][2]));
            if(dist < 3.96 && i != 0){
                if(angleWithin(N[res+i],Od[res],cb[res],.698132)){
                    System.out.println("Amide - hydroxyl H bond at " + (res+start) + " " + seq[res+start-1] + " and backbone N " + (i+start+res));
                hout.println("Amide - hydroxyl H bond at " + (res+start) + " " + seq[res+start-1] + " and backbone N " + (i+start+res));
            }
            }
            if(seq[start + res + i -1] == 'S' || seq[start + res + i -1] == 'T' || seq[start + res + i -1] == 'Y'){
                dist = Math.sqrt((Od[res][0] - Oa[res + i][0])*(Od[res][0] - Oa[res + i][0]) + (Od[res][1] - Oa[res + i][1])*(Od[res][1] - Oa[res + i][1]) + (Od[res][2] - Oa[res + i][2])*(Od[res][2] - Oa[res + i][2]));
                if(dist < 3.76 && i != 0){
                System.out.println("Hydroxyl - hydroxyl H bond at " + (res+start) + " " + seq[res+start-1] + " and " + (res+start+i) + " " + seq[res+start+i-1]);
                hout.println("Hydroxyl - hydroxyl H bond at " + (res+start) + " " + seq[res+start-1] + " and " + (res+start+i) + " " + seq[res+start+i-1]);
            }
            }
            if(seq[start + res + i - 1] == 'Q' || seq[start + res + i - 1] == 'N'){
                dist = Math.sqrt((Od[res][0] - Oa[res + i][0])*(Od[res][0] - Oa[res + i][0]) + (Od[res][1] - Oa[res + i][1])*(Od[res][1] - Oa[res + i][1]) + (Od[res][2] - Oa[res + i][2])*(Od[res][2] - Oa[res + i][2]));
            if(dist < 4.06 && i != 0){
                if(angleWithin(Od[res],ne[res+i],Oa[res+i],.698132)){
                    System.out.println("Hydroxyl - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and " + (res+start+i) + " " + seq[res+start+i-1]);
                    hout.println("Hydroxyl - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and " + (res+start+i) + " " + seq[res+start+i-1]);
                }
            }
            }
        }
        break;
        case 'T':
            for(i = (res - Math.min(30,res)); i < Math.min(30,(END-START-res)); i++){
            dist = Math.sqrt((Od[res][0] - O[res + i][0])*(Od[res][0] - O[res + i][0]) + (Od[res][1] - O[res + i][1])*(Od[res][1] - O[res + i][1]) + (Od[res][2] - O[res + i][2])*(Od[res][2] - O[res + i][2]));
            if(dist < 3.76 && i != 0){
                if(angleWithin(Od[res],cp[res+i],Od[res+i],.698132)){
                    System.out.println("Hydroxyl - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and backbone O " + (i+start+res));
                hout.println("Hydroxyl - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and backbone O " + (i+start+res));
            }
            }
            if(seq[start + res + i -1] == 'S' || seq[start + res + i -1] == 'T' || seq[start + res + i -1] == 'Y'){
                dist = Math.sqrt((Od[res][0] - Oa[res + i][0])*(Od[res][0] - Oa[res + i][0]) + (Od[res][1] - Oa[res + i][1])*(Od[res][1] - Oa[res + i][1]) + (Od[res][2] - Oa[res + i][2])*(Od[res][2] - Oa[res + i][2]));
                if(dist < 3.76 && i != 0){
                System.out.println("Hydroxyl - hydroxyl H bond at " + (res+start) + " " + seq[res+start-1] + " and " + (res+start+i) + " " + seq[res+start+i-1]);
                hout.println("Hydroxyl - hydroxyl H bond at " + (res+start) + " " + seq[res+start-1] + " and " + (res+start+i) + " " + seq[res+start+i-1]);
            }
            }
            if(seq[start + res + i - 1] == 'Q' || seq[start + res + i - 1] == 'N'){
                dist = Math.sqrt((Od[res][0] - Oa[res + i][0])*(Od[res][0] - Oa[res + i][0]) + (Od[res][1] - Oa[res + i][1])*(Od[res][1] - Oa[res + i][1]) + (Od[res][2] - Oa[res + i][2])*(Od[res][2] - Oa[res + i][2]));
            if(dist < 4.06 && i != 0){
                if(angleWithin(Od[res],ne[res+i],Oa[res+i],.698132)){
                    System.out.println("Hydroxyl - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and " + (res+start+i) + " " + seq[res+start+i-1]);
                    hout.println("Hydroxyl - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and " + (res+start+i) + " " + seq[res+start+i-1]);
                }
            }
            }
        }
        break;
        case 'Y':
            for(i = (res - Math.min(30,res)); i < Math.min(30,(END-START-res)); i++){
            dist = Math.sqrt((Od[res][0] - O[res + i][0])*(Od[res][0] - O[res + i][0]) + (Od[res][1] - O[res + i][1])*(Od[res][1] - O[res + i][1]) + (Od[res][2] - O[res + i][2])*(Od[res][2] - O[res + i][2]));
            if(dist < 3.76 && i != 0){
                if(angleWithin(Od[res],cp[res+i],Od[res+i],.698132)){
                    System.out.println("Hydroxyl - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and backbone O " + (i+start+res));
                hout.println("Hydroxyl - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and backbone O " + (i+start+res));
            }
            }
            if(seq[start + res + i -1] == 'S' || seq[start + res + i -1] == 'T' || seq[start + res + i -1] == 'Y'){
                dist = Math.sqrt((Od[res][0] - Oa[res + i][0])*(Od[res][0] - Oa[res + i][0]) + (Od[res][1] - Oa[res + i][1])*(Od[res][1] - Oa[res + i][1]) + (Od[res][2] - Oa[res + i][2])*(Od[res][2] - Oa[res + i][2]));
                if(dist < 3.76 && i != 0){
                System.out.println("Hydroxyl - hydroxyl H bond at " + (res+start) + " " + seq[res+start-1] + " and " + (res+start+i) + " " + seq[res+start+i-1]);
                hout.println("Hydroxyl - hydroxyl H bond at " + (res+start) + " " + seq[res+start-1] + " and " + (res+start+i) + " " + seq[res+start+i-1]);
            }
            }
            if(seq[start + res + i - 1] == 'Q' || seq[start + res + i - 1] == 'N'){
                dist = Math.sqrt((Od[res][0] - Oa[res + i][0])*(Od[res][0] - Oa[res + i][0]) + (Od[res][1] - Oa[res + i][1])*(Od[res][1] - Oa[res + i][1]) + (Od[res][2] - Oa[res + i][2])*(Od[res][2] - Oa[res + i][2]));
            if(dist < 4.06 && i != 0){
                if(angleWithin(Od[res],ne[res+i],Oa[res+i],.698132)){
                    System.out.println("Hydroxyl - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and " + (res+start+i) + " " + seq[res+start+i-1]);
                    hout.println("Hydroxyl - carbonyl H bond at " + (res+start) + " " + seq[res+start-1] + " and " + (res+start+i) + " " + seq[res+start+i-1]);
                }
            }
            }
        }
        break;
        case 'H':
            for(i = (res - Math.min(30,res)); i < Math.min(30,(END-START-res)); i++){
            dist = Math.sqrt((N[res][0] - Na[res + i][0])*(N[res][0] - Na[res + i][0]) + (N[res][1] - Na[res + i][1])*(N[res][1] - Na[res + i][1]) + (N[res][2] - Na[res + i][2])*(N[res][2] - Na[res + i][2]));
            if(dist < 3.96 && i !=0){
                if(angleWithin(N[res],ne[res],Na[res+i],.698132)){
                    System.out.println("Amide - imidazole H bond at " + (res+start) + " " + seq[res+start-1] + " and backbone N " + (res+start+i));
                hout.println("Amide - imidazole H bond at " + (res+start) + " " + seq[res+start-1] + " and backbone N " + (res+start+i));
            }
            }
        }
        break;
        default: //System.out.println("No H bond on " + (res+start));
    }
    
    return 0;
    }
    
    /* angleWithin
     * takes: double[], double[], double[], double
     * returns boolean
     * input: (donor position, heavy atom angle refrence position, acceptor position, angle it must be less than (radians)
     *
     * This function returns true if the vector from b to a makes an angle less than a specified angle
     * with the vector that points from b to c, and false otherwise
     */
    private static boolean angleWithin(double a[], double b[], double c[], double angle){
        double vect[][] = new double [2][3];
    double dot = 0.0;
    double mag1, mag2;
    double actAng;
    int i;
    for(i = 1; i < 3; i++){
        vect[0][i] = a[i] - b[i];
        vect[1][i] = c[i] - b[i];
    }
    mag1 = Math.sqrt(vect[0][0]*vect[0][0] + vect[0][1]*vect[0][1] + vect[0][2]*vect[0][2]);
    mag2 = Math.sqrt(vect[1][0]*vect[1][0] + vect[1][1]*vect[1][1] + vect[1][2]*vect[1][2]);
    dot = vect[0][0]*vect[1][0] + vect[0][1]*vect[1][1] + vect[0][2]*vect[1][2];
    actAng = Math.acos(dot/(mag1*mag2));
    if( Math.abs(actAng) < angle){
        return true;
    }else{
        return false;
    }
    }
    
    
    /*
     * findVolApp
     * takes: double[][], int, int
     * returns double
     * input: (C alpha position array, starting residue, ending residue)
     *
     * This will find an approximation to the volume by assuming that the helix is a triangular cylinder
     * it finds the value of the longest side of the (effectivly) equiladiral triangle
     * and then compute the volume by using (3)^(1/2)/4*L^2*H
     * where L is the length of a side of the equiladiral triangle and
     * where H is the hieght of the structure
     */
    private static double findVolApp(double CAgrid[][], int START, int END){
        int i, j = 0, Np;
    double MaxDist = 0.0;
    double dist, turns = 0.0;
    // finds the distance away in 1D of the closest vertical residue in 3D to the first residue we do this to know how far in 1D we have to go to make a ring in 3D
    Np = getNplus(0,3,0,9,Math.min(30,END),CAgrid);
    for(i = 0; i < Np; i++){ // loops over the ring described above
        for(j = 0; j < Np; j++){ // searches through again (so we find the max dist between ANY two points on the ring
            // computes the euclidian distace between the two points
            dist = Math.sqrt((CAgrid[i][0] - CAgrid[j][0])*(CAgrid[i][0] - CAgrid[j][0]) + (CAgrid[i][1] - CAgrid[j][1])*(CAgrid[i][1] - CAgrid[j][1]) + (CAgrid[i][2] - CAgrid[j][2])*(CAgrid[i][2] - CAgrid[j][2]));
            if(dist > MaxDist){ // if its bigger than the current maximum distance then it replaces it as the maximum distance
            MaxDist = dist; // this allows the maximum distance between any two points to "bubble sort" itself to the top
            }
        }
    }
    turns = ((END-START)/18.0); // approximates 18 turns per amino acid
    System.out.println("Edge length: " + MaxDist); // output for debuging pourposes
    return Math.sqrt(3)/4.0*MaxDist*MaxDist*(turns*4.8); // returns the approximated volume
    }
    
    /*
     * getPhi
     * takes: int, double[][], double[][], double[][], int
     * returns double
     * input: (residue that you want the phi angle for, C alpha position array, N position array, C' position array, length of the 1D string of residues that make up the helix (END-START)
     *
     * This finds the phi angle defined by the angle between the two planes that are defined as
     * plane 1: the plane in 3-space consisting of the positions of the N, C' (from the previous residue) and CA atoms
     * plane 2: the plane in 3-space consisting of the positions of the C alpha, N and C' atoms
     *
     * for more information look up information regarding a Ramachandran Plot (page 48 of Protiens)
     */
    private static double getPhi(int res, double CAgrid[][], double Ngrid[][], double CPgrid[][], int len){
    double CpN[] = new double[3]; // the vector that points from N to C'
    double CAN[] = new double[3]; // the vector that points from N to CA
    double cp1[] = new double[3]; // the cross product of the above two vectors
    double NCA[] = new double[3]; // the vector that points from CA to N
    double CPCA[] = new double[3]; // the vector that points from CA to C'
    double cp2[] = new double[3]; // the cross product of the above two vectors
    double dot = 0.0;
    double cp1mag = 0.0, cp2mag = 0.0;
    double phi = 0.0;
    if(res != 0){
        // Assign vectors based off of the grids
        CpN[0] = CPgrid[res-1][0] - Ngrid[res][0];
        CpN[1] = CPgrid[res-1][1] - Ngrid[res][1];
        CpN[2] = CPgrid[res-1][2] - Ngrid[res][2];
        CAN[0] = CAgrid[res][0] - Ngrid[res][0];
        CAN[1] = CAgrid[res][1] - Ngrid[res][1];
        CAN[2] = CAgrid[res][2] - Ngrid[res][2];
        NCA[0] = Ngrid[res][0] - CAgrid[res][0];
        NCA[1] = Ngrid[res][1] - CAgrid[res][1];
        NCA[2] = Ngrid[res][2] - CAgrid[res][2];
        CPCA[0] = CPgrid[res][0] - CAgrid[res][0];
        CPCA[1] = CPgrid[res][1] - CAgrid[res][1];
        CPCA[2] = CPgrid[res][2] - CAgrid[res][2];
        // compute cross products CAN x CpN
        cp1[0] = CAN[1]*CpN[2] - CpN[1]*CAN[2];
        cp1[1] = CAN[2]*CpN[0] - CpN[2]*CAN[0];
        cp1[2] = CAN[0]*CpN[1] - CpN[0]*CAN[1];
        // CPCA x NCA
        cp2[0] = CPCA[1]*NCA[2] - NCA[1]*CPCA[2];
        cp2[1] = CPCA[2]*NCA[0] - NCA[2]*CPCA[0];
        cp2[2] = CPCA[0]*NCA[1] - NCA[0]*CPCA[1];
        // compute dot product of two cross products
        dot = cp1[0]*cp2[0] + cp1[1]*cp2[1] + cp1[2]*cp2[2];
        // compute magnitude of two cross products
        cp1mag = Math.sqrt((cp1[0])*(cp1[0]) + (cp1[1])*(cp1[1]) + (cp1[2])*(cp1[2]));
        cp2mag = Math.sqrt((cp2[0])*(cp2[0]) + (cp2[1])*(cp2[1]) + (cp2[2])*(cp2[2]));
        // compute angles using dot product rule cp1*cp2 = |cp1||cp2|cos(theta)
        phi = Math.acos(dot/(cp1mag*cp2mag));
        /*
         * These next lines make phi into (-pi,pi) form this is done by
         * checking whether or not the angle that it sweeps out is past
         * 6 O'clock. To do this we dot the normal vector of the plane that
         * is stationary (which will give us a vector at 3 O'clock) with the
         * vector that extends from the C alpha to the C'. When this is < 0
         * then it is past 6 O'clock and we must transform phi.
         */
        if((cp1[0]*CPCA[0] + cp1[1]*CPCA[1] + cp1[2]*CPCA[2]) < 0){
            phi = (-1.0)*phi;
        }
        // return answer
        return phi;  
    }else{
        return 0.0;
    }
    
    }
    
    /*
     * getPsi
     * takes: int, double[][], double[][], double[][], int
     * returns double
     * input: (residue that you want the psi angle for, C alpha position array, N position array, C' position array, length of the 1D string of residues that make up the helix (END-START)
     *
     * This finds the psi angle defined by the angle between the two planes that are defined as
     * plane 1: the plane in 3-space consisting of the positions of the N (from the following residue), C' and CA atoms
     * plane 2: the plane in 3-space consisting of the positions of the C alpha, N and C' atoms
     *
     * for more information look up information regarding a Ramachandran Plot (Page 48 of Protiens)
     */
    private static double getPsi(int res, double CAgrid[][], double Ngrid[][], double CPgrid[][], int len){
    double NCA[] = new double[3]; // the vector that points from CA to N
    double CPCA[] = new double[3]; // the vector that points from CA to C'
    double cp2[] = new double[3]; // the cross product of the above two vectors
    double NCp[] = new double[3]; // the vector that points from C' to N
    double CACP[] = new double[3]; // the vector that points from C' to CA
    double cp1[] = new double[3]; // the cross product of the above two vectors
    double dot = 0.0;
    double cp1mag = 0.0, cp2mag = 0.0;
    double psi = 0.0;
    len--;
    if(res != len){
        // Assign vectors based off of the grids
        NCp[0] = Ngrid[res+1][0] - CPgrid[res][0];
        NCp[1] = Ngrid[res+1][1] - CPgrid[res][1];
        NCp[2] = Ngrid[res+1][2] - CPgrid[res][2];
        CACP[0] = CAgrid[res][0] - CPgrid[res][0];
        CACP[1] = CAgrid[res][1] - CPgrid[res][1];
        CACP[2] = CAgrid[res][2] - CPgrid[res][2];
        NCA[0] = Ngrid[res][0] - CAgrid[res][0];
        NCA[1] = Ngrid[res][1] - CAgrid[res][1];
        NCA[2] = Ngrid[res][2] - CAgrid[res][2];
        CPCA[0] = CPgrid[res][0] - CAgrid[res][0];
        CPCA[1] = CPgrid[res][1] - CAgrid[res][1];
        CPCA[2] = CPgrid[res][2] - CAgrid[res][2];
        // compute cross products NCp x CACP
        cp1[0] = NCp[1]*CACP[2] - CACP[1]*NCp[2];
        cp1[1] = NCp[2]*CACP[0] - CACP[2]*NCp[0];
        cp1[2] = NCp[0]*CACP[1] - CACP[0]*NCp[1];
        // CPCA x NCA
        cp2[0] = CPCA[1]*NCA[2] - NCA[1]*CPCA[2];
        cp2[1] = CPCA[2]*NCA[0] - NCA[2]*CPCA[0];
        cp2[2] = CPCA[0]*NCA[1] - NCA[0]*CPCA[1];
        // compute dot product of two cross products
        dot = cp1[0]*cp2[0] + cp1[1]*cp2[1] + cp1[2]*cp2[2];
        // compute magnitude of two cross products
        cp1mag = Math.sqrt((cp1[0])*(cp1[0]) + (cp1[1])*(cp1[1]) + (cp1[2])*(cp1[2]));
        cp2mag = Math.sqrt((cp2[0])*(cp2[0]) + (cp2[1])*(cp2[1]) + (cp2[2])*(cp2[2]));
        // compute angles using dot product rule cp1*cp2 = |cp1||cp2|cos(theta)
        psi = Math.acos(dot/(cp1mag*cp2mag));
        /*
         * These next lines make psi into (-pi,pi) form this is done by
         * checking whether or not the angle that it sweeps out is past
         * 6 O'clock. To do this we dot the normal vector of the plane that
         * is stationary (which will give us a vector at 3 O'clock) with the
         * vector that extends from the C' to the N. When this is < 0
         * then it is past 6 O'clock and we must transform psi.
         */
        if((cp2[0]*NCp[0] + cp2[1]*NCp[1] + cp2[2]*NCp[2]) < 0){
            psi = (-1.0)*psi;
        }
        // return answer
        return psi;  
    }else{
        return 0.0;
    }
    
    }
    
    /*
     * findVolume
     * takes: double[][], int, boolean, boolean
     * returns double
     * input (C alpha position array, length of 1D helix (END-START), strictness of the in/out test: true includes planar caps to the helix, is the helix left handed?)
     *
     * This function computes the volume of the helix using Monte Carlo integration
     * It finds the range of space the the helix is in and generates a ranfom test point
     * it then finds the closest C alpha to that test point and tests whether or not it is inside the helix
     * by pretending that the point is that C alphas corresponding C beta point
     * If it is inside the helix then it is added to the count
     *
     * at the end we use the monte carlo integration formula that V = (in/total)*(space of possible values)
     */
    private static double findVolume(double CAgrid[][], int len, boolean strict, boolean isLH){
    double in = 0;
        int closest = 0;
        double test[] = new double[3];
        double max[] = new double[3];
    double min[] = new double[3];
        double volume = 0.0;
        max[0] = -9000.0;
        max[1] = -9000.0;
        max[2] = -9000.0;
        int i, j;
        double minDist = 1000000.0, dist = 0.0;
    double minMin = -1000000.0;
        
    // Searches for the maximum range in all three orthogonal directions that the C alpha position array includes
        for(i = 0; i < len; i++){
            if(CAgrid[i][0] > max[0]){
                max[0] = CAgrid[i][0];
            }
            if(CAgrid[i][1] > max[1]){
                max[1] = CAgrid[i][1];
            }
            if(CAgrid[i][2] > max[2]){
                max[2] = CAgrid[i][2];
            }
        }
    
    // Searches for the minimum range in all three directions that the C alpha position array includes
    min[0] = max[0];
    min[1] = max[1];
    min[2] = max[2];
    for(i = 0; i < len; i++){
        if(CAgrid[i][0] < min[0]){
            min[0] = CAgrid[i][0];
        }
        if(CAgrid[i][1] < min[1]){
            min[1] = CAgrid[i][1];
        }
        if(CAgrid[i][2] < min[2]){
            min[2] = CAgrid[i][2];
        }
    }
    // We now have the maximum and minimum values for x,y,z
    
    // Creates a new random number generator
        Random rand1 = new Random();
        
    // loops over 500000 Monte Carlo points
        for(i = 0; i < 500000; i++){
            test[0] = rand1.nextDouble()*(max[0]-min[0]) + min[0]; // gets test x value
            test[1] = rand1.nextDouble()*(max[1]-min[1]) + min[1]; // gets test y value
            test[2] = rand1.nextDouble()*(max[2]-min[2]) + min[2]; // gets test z value
        
        // finds the closest C alpha point to the test value using a "bubble sort"
            for(j = 0; j < len; j++){
                dist = Math.sqrt((test[0]-CAgrid[j][0])*(test[0]-CAgrid[j][0]) + (test[1]-CAgrid[j][1])*(test[1]-CAgrid[j][1]) + (test[2]-CAgrid[j][2])*(test[2]-CAgrid[j][2]));
                if(dist < minDist){
                    minDist = dist;
                    closest = j;
                }
            }
        // determines whether or not the test value is inside the helix by pretending that it is a C beta to the closest C alpha (for info on strict see CBisINside)
            if(CBisInside(closest, test[0], test[1], test[2], CAgrid, len, strict,isLH)){
                in = in + 1.0;
            }
            minDist = 1000000.0;
        }
    // the monte carlo integration formula that V = (in/total)*(space of possible values)
        volume = (in/500000.0)*(max[0]-min[0])*(max[1]-min[1])*(max[2]-min[2]);
        return volume;
    }
    
    /*
     * getNplus
     * takes: int, int, int, int, int, double[][]  // That is 5 ints
     * returns int
     * input: (residue that you want Nplus for, minimum distance away it will evaluate against (3), dummy variable (0), where it starts when looking for the residue (9), how deep into the 1D chain it will search (must be less than the number of residues left), C alpha postion array)
     *
     * This function searches for and returns the number of residues in between it and the residue directly above it in the helix
     * It starts by looking at some residue a minumum distance away and compares it to the distance between it and another residue which is relivly close
     * if it is closer than the criteria then it becomes the new closeset and it tries to find a residue which is closer yet
     * if there is no residue above it (if its at the top of the helix) then it returns a value of zero
     */
    private static int getNplus(int res, int min, int minsearch, int minstart, int searchDomP, double CAgrid[][]){
        int j = res, i;
    double dist = 1000.0;
    double minDist = 1000.0;
    int Nplus = 0, Nplusp = 0;
    boolean npOK = false;
    //System.out.println("TEST");
    //System.out.println(res + " " + min + " " + minsearch + " " + minstart + " " + searchDomP);
    for(minsearch = minstart; minsearch < searchDomP; minsearch++){
            Nplusp = Nplus;
            //if(searchDomP > minsearch){
                for(i = minsearch; i < searchDomP; i++){
                    dist = Math.sqrt((CAgrid[j][0]-CAgrid[j+i][0])*(CAgrid[j][0]-CAgrid[j+i][0])+(CAgrid[j][1]-CAgrid[j+i][1])*(CAgrid[j][1]-CAgrid[j+i][1])+(CAgrid[j][2]-CAgrid[j+i][2])*(CAgrid[j][2]-CAgrid[j+i][2]));
                    //System.out.println(dist);
            if(dist < minDist){
                        minDist = dist;
                        Nplus = i;
            //System.out.println("Lowest " + Nplus);
                    }
                }
                if(Math.sqrt((CAgrid[j][0]-CAgrid[j+min-1][0])*(CAgrid[j][0]-CAgrid[j+min-1][0])+(CAgrid[j][1]-CAgrid[j+min-1][1])*(CAgrid[j][1]-CAgrid[j+min-1][1])+(CAgrid[j][2]-CAgrid[j+min-1][2])*(CAgrid[j][2]-CAgrid[j+min-1][2]))<minDist){
                    Nplus = 0;
                }
            //}else{
            //    Nplus = 0;
            //}
        
            if(minsearch != minstart){
                if(Nplusp == Nplus){
                    npOK = true;
            //i = searchDomP;
                    minsearch = searchDomP;
                }else{
                    npOK = false;
                }
            }
            minDist = 1000.0;
        }
        if(npOK == false){
            Nplus = 0;
        }
    //System.out.println("Lowest " + Nplus);
    return Nplus;
    }
    
    /*
     * getNminus
     * takes int, int, int, int, int, double[][]
     * returns int
     * input: (residue that you want Nplus for, minimum distance away it will evaluate against (3), dummy variable (0), where it starts when looking for the residue (9), how deep into the 1D chain it will search (must be less than the number of residues left), C alpha postion array)
     *
     * This works in exactly the same way as getNplus only in the opposite direction
     */
    private static int getNminus(int res, int min, int minsearch, int minstart, int searchDomM, double CAgrid[][]){
        int j = res, i;
    double dist = 1000.0;
    double minDist = 1000.0;
    int Nminus = 0, Nminusp = 0;
    boolean nmOK = false;
    //System.out.println("TESTL");
    for(minsearch = minstart; minsearch < searchDomM; minsearch++){ 
            Nminusp = Nminus;
            //if(searchDomM > minsearch){
                for(i = minsearch; i < searchDomM; i++){
                    dist = Math.sqrt((CAgrid[j][0]-CAgrid[j-i][0])*(CAgrid[j][0]-CAgrid[j-i][0])+(CAgrid[j][1]-CAgrid[j-i][1])*(CAgrid[j][1]-CAgrid[j-i][1])+(CAgrid[j][2]-CAgrid[j-i][2])*(CAgrid[j][2]-CAgrid[j-i][2]));
                    //System.out.println(dist);
            if(dist < minDist){
                        minDist = dist;
                        Nminus = i;
            //System.out.println("Lowest");
                    }
                }
                if(Math.sqrt((CAgrid[j][0]-CAgrid[j-min-1][0])*(CAgrid[j][0]-CAgrid[j-min-1][0])+(CAgrid[j][1]-CAgrid[j-min-1][1])*(CAgrid[j][1]-CAgrid[j-min-1][1])+(CAgrid[j][2]-CAgrid[j-min-1][2])*(CAgrid[j][2]-CAgrid[j-min-1][2]))<minDist){
                    Nminus = 0;
                }
            //}else{
            //    Nminus = 0;
            //}
            if(minsearch != minstart){
                if(Nminusp == Nminus){
                    nmOK = true;
            //i = searchDomM;
                    minsearch = searchDomM;
                }else{
                    nmOK = false;
                }
            }
            minDist = 1000.0;
        } 
        if(nmOK == false){
            Nminus = 0;
        }
    return Nminus;
    }
    
    /*
     * CBisInside
     * takes int, double, double, double, double[][], int, boolean, boolean
     * returns boolean
     * input (the residue that we are comparing the test point to (usually the closest c alpha),
     *      the x coord of the test point, the y coord of the test point, the z coord of the test point,
     *      the position array of the C alpha atoms, the 1D length of helix (END-START),
     *      the strictness of the algorithm: true will include planar caps to the helix,
     *      whether or not the helix is left handed (true) or right handed (false))
     *
     * This function will test whether a single point in 3-space lies inside the helix
     * to do this 5 refrence points are found, the closest C alpha, the C alphas on either side of that,
     * and the C alpha right above and right below that closest C alpha
     * vectors are then made to point from that closest C alpha towards the other 4 points
     * these vectors make up 4 planes that are an approximation to the side of the helix
     * the cross product of these vectors (if taken in the right order) will produce a
     * vector that will always point inwards (due to the consistant left handed nature of the helix)
     * we then calculate the dot product of the vector that points from the closest C alpha to the C beta
     * (or any point that we wish to test) and the cross product vector from one of the planes
     * if this dot product is positive then the vector points in and the point (or C beta and therefore
     * the side chain) is inside the helix. We then repeat this step for the other (up to 3) sides.
     * The vector must pass a majority (>50%) of these tests for it to be considered "in."
     * There are some circumstances where there is no upper C alpha or no C alpha to the right etc.
     * for these situations only the possible tests are done.
     *
     * If strict == true then the function takes another step and puts another restriction on the test
     * planar caps are formed at the top and bottom of the helix in order to determine whether or not
     * the point is oriented inwards but is too far "above" the helix to be counted "in."
     * The caps are constructed by finding points on the top and bottom rungs of the helix that are in
     * the corners. A equation of a plane is then found using these points. We then test a point that is
     * known to be on the correct (in) side of the plane and make sure that all subsequent points follow
     * the same inequality. If either of these tests fail then the particle imeadiatly fails. This strict
     * parameter does not need to be invoked when simply checking whether a C beta is inside or not, becuase
     * one can be reasonably certain that orientation will dictate whether the side chain is in or not.
     * When performing the Monte Carlo integration on the other hand, points are not restricted to the same
     * geometry as teh side chains and without the strictness points above and below the helix (enough so
     * to be deemed out logically) will still be counted "in" based on orientation alone if the planar caps
     * are not inforced.
     *
     * If the helix is left handed then it will orient the cross products so that they point in and check the
     * test vector against them in a positive way (see if they DO point along the same line). If the helix is
     * right handed then it will orient the cross products so that the test vector wants to point in the opposite
     * direction to be counted in. Left handed = true, Right handed = false.
     */
    private static boolean CBisInside(int res, double x, double y, double z, double CAgrid[][], int len, boolean strict, boolean isLH){
        int i = 0, j = res, k;
        int searchDomP = Math.min(50,len - res); // This determines how far away the program will search for the C alpha above
    //_ the closest one, minimum of 50 residues except when length does not permit, in which case as far as possible
        int searchDomM = Math.min(50,res); // same as above only for the C alpha below
        double CACAp1[] = new double[3]; // Vector that points from the close C alpha to the next C alpha in the series
        double CACAm1[] = new double[3]; // Vector that points from the close C alpha to the previous one in the series
        double CACApN[] = new double[3]; // Vector that points from the close C alpha to the one directly above it
        double CACAmN[] = new double[3]; // Vector that points from the close C alpha to the one directly below it
        double CACB[] = new double[3]; // Vector that points from the close C alpha to the test point
    double cp[] = new double[3]; // Generic cross product vector
        double dot = 0.0;
    
    // Variables and conditions for searching for the C alpha above and below
    int Nplus, Nminus;
        int tests = 0, in = 0;
        int min = 3, minsearch = 0, minstart;
        minstart = min + 6;
        len--;
        
    // Defining the vector that points from the closest C alpha to the test point
        CACB[0] = x - CAgrid[j][0];
        CACB[1] = y - CAgrid[j][1];
        CACB[2] = z - CAgrid[j][2];
        
    // finds how many residues away the C alpha directly above is
        Nplus = getNplus(res, min, minsearch, minstart, searchDomP, CAgrid);
    // and then defines the vector from the closest C alpha to it
        if(Nplus != 0){ // but only if it exists
            CACApN[0] = CAgrid[j+Nplus][0] - CAgrid[j][0];
            CACApN[1] = CAgrid[j+Nplus][1] - CAgrid[j][1];
            CACApN[2] = CAgrid[j+Nplus][2] - CAgrid[j][2];
        }

    // finds how many residues away the C alpha directly below the closest is
        Nminus = getNminus(res, min, minsearch, minstart, searchDomM, CAgrid);
    // and defines the vector if it exists
        if(Nminus != 0){
            CACAmN[0] = CAgrid[j-Nminus][0] - CAgrid[j][0];
            CACAmN[1] = CAgrid[j-Nminus][1] - CAgrid[j][1];
            CACAmN[2] = CAgrid[j-Nminus][2] - CAgrid[j][2];
        }
        
    // if the closest residue is not the last one in the sequence
    // then it defines the vector that points from the closest C alpha to the one in front of it
        if(res != len){
            CACAp1[0] = CAgrid[j+1][0] - CAgrid[j][0];
            CACAp1[1] = CAgrid[j+1][1] - CAgrid[j][1];
            CACAp1[2] = CAgrid[j+1][2] - CAgrid[j][2];
        }
        
    // if the closest residue is not the first one in the sequence
    // then it defines the vector that points from the closest C alpha to the one behind it
        if(res != 0){
            CACAm1[0] = CAgrid[j-1][0] - CAgrid[j][0];
            CACAm1[1] = CAgrid[j-1][1] - CAgrid[j][1];
            CACAm1[2] = CAgrid[j-1][2] - CAgrid[j][2];

        }
    /* debuging outputs NOT EXECUTED WHILE COMENTED
        System.out.println("CACAp1: " + CACAp1[0] + " " + CACAp1[1] + " " + CACAp1[2]);
        System.out.println("CACApN: " + CACApN[0] + " " + CACApN[1] + " " + CACApN[2]);
        System.out.println("pN: " + Nplus);
        System.out.println("CACAm1: " + CACAm1[0] + " " + CACAm1[1] + " " + CACAm1[2]);
        System.out.println("CACAmN: " + CACAmN[0] + " " + CACAmN[1] + " " + CACAmN[2]);
        System.out.println("mN: " + Nminus);
        System.out.println("CACB " + CACB[0] + " " + CACB[1] + " " + CACB[2]);
        System.out.println("");*/
        //System.out.println(Nplus + " " + Nminus);
    
        //upper right test (CACAp1 and CACApN)
        if(Nplus != 0 && res != 0){
            tests++; // determines that a test has been done (for >50% check)
        // Computes the cross product of the two vectors
            cp[0] = CACAp1[1]*CACApN[2] - CACApN[1]*CACAp1[2];
            cp[1] = CACAp1[2]*CACApN[0] - CACApN[2]*CACAp1[0];
            cp[2] = CACAp1[0]*CACApN[1] - CACApN[0]*CACAp1[1];
        // Computes the dot product of the two vectors
            dot = cp[0]*CACB[0] + cp[1]*CACB[1] + cp[2]*CACB[2];
        // Determines whether or not it is in
            if(dot > 0){
            if(isLH){
                    in++; // and counts it if it is in
        }
            }else{
            if(!isLH){
            in++;
        }
            }
        }
        //upper left
        if(Nplus != 0 && res != len){
            tests++; // determines that a test has been done (for >50% check)
        // Computes the cross product of the two vectors
            cp[0] = CACApN[1]*CACAm1[2] - CACAm1[1]*CACApN[2];
            cp[1] = CACApN[2]*CACAm1[0] - CACAm1[2]*CACApN[0];
            cp[2] = CACApN[0]*CACAm1[1] - CACAm1[0]*CACApN[1];
            // Computes the dot product of the two vectors
            dot = cp[0]*CACB[0] + cp[1]*CACB[1] + cp[2]*CACB[2];
            // Determines whether or not it is in
            if(dot > 0){
            if(isLH){
                    in++; // and counts it if it is in
        }
            }else{
            if(!isLH){
            in++;
        }
            }
        }
        //lower left
        if(Nminus != 0 && res != 0){
            tests++; // determines that a test has been done (for >50% check)
        // Computes the cross product of the two vectors
            cp[0] = CACAm1[1]*CACAmN[2] - CACAmN[1]*CACAm1[2];
            cp[1] = CACAm1[2]*CACAmN[0] - CACAmN[2]*CACAm1[0];
            cp[2] = CACAm1[0]*CACAmN[1] - CACAmN[0]*CACAm1[1];
            // Computes the dot product of the two vectors
            dot = cp[0]*CACB[0] + cp[1]*CACB[1] + cp[2]*CACB[2];
            // Determines whether or not it is in
            if(dot > 0){
            if(isLH){
                    in++; // and counts it if it is in
        }
            }else{
            if(!isLH){
            in++;
        }
            }
        }
        //lower right
        if(Nminus != 0 && res != len){
            tests++; // determines that a test has been done (for >50% check)
        // Computes the cross product of the two vectors
            cp[0] = CACAmN[1]*CACAp1[2] - CACAp1[1]*CACAmN[2];
            cp[1] = CACAmN[2]*CACAp1[0] - CACAp1[2]*CACAmN[0];
            cp[2] = CACAmN[0]*CACAp1[1] - CACAp1[0]*CACAmN[1];
            // Computes the dot product of the two vectors
            dot = cp[0]*CACB[0] + cp[1]*CACB[1] + cp[2]*CACB[2];
            // Determines whether or not it is in
            if(dot > 0){
            if(isLH){
                    in++; // and counts it if it is in
        }
            }else{
            if(!isLH){
            in++;
        }
            }
        }
        
    // The following code will test to see if the chain is in a complicated side loop
    // If there is no C alpha directly above and below but the residue is in the middle of
    // the sequence then it must be in a side loop and therefore cannot be "in"
        if((Nplus == 0 || Nminus == 0) && (j+25 < len && j-25 > 0)){
            return false;
        }
    
    // Strictness (puts planar caps on the helix)
    if(strict == true){
        double Corner[][] = new double[3][3];
        double P1[] = new double[3];
        double d = 0.0, v = 0.0;
        k = 0;
        
        // BOTTOM CAP
        // finds the first 3 corners
        for(i = 1; i < len; i++){
            if(onCorner(i,CAgrid)){
            Corner[k][0] = CAgrid[i][0] - CAgrid[0][0];
            Corner[k][1] = CAgrid[i][1] - CAgrid[0][1];
            Corner[k][2] = CAgrid[i][2] - CAgrid[0][2];
            k++;
        }
        if(k == 3){
            i = len;
        }
        }
        // finds the orthogonal vector to the plane definded by the first two corners and the origin
        P1[0] = Corner[1][0]*Corner[0][1] - Corner[0][0]*Corner[1][1];
        P1[1] = Corner[1][1]*Corner[0][2] - Corner[0][1]*Corner[1][2];
        P1[2] = Corner[1][2]*Corner[0][0] - Corner[0][2]*Corner[1][0];
        // determines the constant d for the plane function
        d = (-1.0)*CAgrid[0][0]*P1[0] - CAgrid[0][1]*P1[1] - CAgrid[0][2]*P1[2];
        // orients what is meant by "in" by the plane inequality
        v = P1[0]*Corner[2][0] + P1[1]*Corner[2][1] + P1[2]*Corner[2][2] + d;
        // searches for points that fall "below" the bottom cap
        if( v < 0 ){
            if ( 0 < P1[0]*x + P1[1]*y + P1[2]*z + d ){
            return false;
        }
        }
        if( v >= 0 ){
            if ( 0 > P1[0]*x + P1[1]*y + P1[2]*z + d ){
            return false;
        }
        }
        
        // TOP CAP
        k = 0;
        // finds the last 3 corners
        for(i = len-2; i > 0; i--){
            if(onCorner(i,CAgrid)){
            Corner[k][0] = CAgrid[i][0] - CAgrid[len-1][0];
            Corner[k][1] = CAgrid[i][1] - CAgrid[len-1][1];
            Corner[k][2] = CAgrid[i][2] - CAgrid[len-1][2];
            k++;
        }
        if(k == 3){
            i = 0;
        }
        }
        // finds the orthogonal vector to the plane defined by the last two corners and the end
        P1[0] = Corner[1][0]*Corner[0][1] - Corner[0][0]*Corner[1][1];
        P1[1] = Corner[1][1]*Corner[0][2] - Corner[0][1]*Corner[1][2];
        P1[2] = Corner[1][2]*Corner[0][0] - Corner[0][2]*Corner[1][0];
        // determines the constant d for the plane function
        d = (-1.0)*CAgrid[0][0]*P1[0] - CAgrid[0][1]*P1[1] - CAgrid[0][2]*P1[2];
        // orients what is meant by "in" by the plane inequality
        v = P1[0]*Corner[2][0] + P1[1]*Corner[2][1] + P1[2]*Corner[2][2] + d;
        // searches for the points that fall "above" the top cap
        if( v < 0 ){
            if ( 0 < P1[0]*x + P1[1]*y + P1[2]*z + d ){
            return false;
        }
        }
        if( v >= 0 ){
            if ( 0 > P1[0]*x + P1[1]*y + P1[2]*z + d ){
            return false;
        }
        }
    }
        
        if(((double)in)/((double)tests) > .5){
            return true;
        }else{
            return false;
        }           
    }
    
    /*
     * classifyPhiPsi
     * takes double, double
     * returns int
     * input: (phi angle, psi angle)
     *
     * This function takes in the phi and psi angle for a residue and returns
     * what type of secondary structure that the residue is in
     * 0: Beta
     * 1: Right handed alpha
     * 2: Left handed alpha
     *
     * It determines this by using the plot on page 49 of Whitford's Protiens
     * Beta: phi < 0, psi > 1 or psi < -2
     * RH Alpha: phi < 0, -2 < psi < 1
     * LH Alpha: phi > 0
     */
    private static int classifyPhiPsi(double phi, double psi){
        if(phi > 0){
        return 2;
    }
    if(psi > -2.0){
        if(psi < 1.0){
            return 1;
        }
    }
    return 0;
    }
          
    /*
     * onCorner
     * takes int, double[][]
     * returns boolean
     * input: (residue that you want to check, the C alpha position array)
     *
     * This function checks to see if a residue is on the corner of the helix
     * it does this by making vector that point from the suspect residue to
     * the residue before and after it. Then by using the dot product rule:
     *      u*v = |u||v|cos(theta)
     * if the angle between them is less than 110 deg or 1.920 radians then
     * it is determined that the residue is on a corner and the value
     * true is returned.
     */
    private static boolean onCorner(int res, double CAgrid[][]){
        double CE[] = new double[3];
    double CS[] = new double[3];
    double dot = 0.0;
    double CEm = 0.0, CSm = 0.0;
    CE[0] = CAgrid[res][0] - CAgrid[res-1][0];
    CE[1] = CAgrid[res][1] - CAgrid[res-1][1];
    CE[2] = CAgrid[res][2] - CAgrid[res-1][2];
    CS[0] = CAgrid[res][0] - CAgrid[res+1][0];
    CS[1] = CAgrid[res][1] - CAgrid[res+1][1];
    CS[2] = CAgrid[res][2] - CAgrid[res+1][2];
    dot = CE[0]*CS[0] + CE[1]*CS[1] + CE[2]*CS[2];
    CEm = Math.sqrt(CE[0]*CE[0] + CE[1]*CE[1] + CE[2]*CE[2]);
    CSm = Math.sqrt(CS[0]*CS[0] + CS[1]*CS[1] + CS[2]*CS[2]);
    if(Math.abs(Math.acos(dot/(CEm*CSm))) < 1.920){
        return true;
    }else{
        return false;
    }
    }
    
    /* getVol
     * takes char
     * returns double
     * input: (The character code for the amino acid that you want the volume for)
     *
     * This function will take in the letter code for an amino acid and return
     * the VDW volume of the amino acid as defined by Whitford in Protiens on page 18-22
     * This includes the backbone (GLY as a non-zero volume)
     */
    private static double getVol(char r){
        switch (r) {
            case 'A': return 67.0;
            case 'R': return 167.0;
            case 'N': return 148.0;
            case 'D': return 67.0;
            case 'C': return 86.0;
            case 'Q': return 114.0;
            case 'E': return 109.0;
            case 'G': return 48.0;
            case 'H': return 118.0;
            case 'I': return 124.0;
            case 'L': return 124.0;
            case 'K': return 135.0;
            case 'M': return 124.0;
            case 'F': return 135.0;
            case 'P': return 90.0;
            case 'S': return 73.0;
            case 'T': return 93.0;
            case 'W': return 163.0;
            case 'Y': return 141.0;
            case 'V': return 105.0;
            default: System.out.println("Error in getVol");
        }
        return 0.0;
    }
    
    /*
     * getNum
     * takes: char
     * returns: int
     * input: (The character code for an amino acid)
     *
     * This function simply orders the amino acids in the same way that Whitford does in Protiens
     * and returns the number value of thier position starting from 0
     * This has no scientific use, it is merely a trick to make some arrays more compact and simple
     */
    private static int getNum(char r){
        switch (r) {
            case 'A': return 0;
            case 'R': return 1;
            case 'N': return 2;
            case 'D': return 3;
            case 'C': return 4;
            case 'Q': return 5;
            case 'E': return 6;
            case 'G': return 7;
            case 'H': return 8;
            case 'I': return 9;
            case 'L': return 10;
            case 'K': return 11;
            case 'M': return 12;
            case 'F': return 13;
            case 'P': return 14;
            case 'S': return 15;
            case 'T': return 16;
            case 'W': return 17;
            case 'Y': return 18;
            case 'V': return 19;
            default: System.out.println("Error in getNum");
        }
        return 0;
    }
    
    /*
     * getName
     * takes: int
     * returns: char
     * input: (The number value of an amino acid)
     *
     * This function basically undoes getNum for when the data needs to be read
     * logically from the array. It does nothing without first using getNum.
     */
    private static int getName(int i){
        switch(i){
        case 0: return 'A';
        case 1: return 'R';
        case 2: return 'N';
        case 3: return 'D';
        case 4: return 'C';
        case 5: return 'Q';
        case 6: return 'E';
        case 7: return 'G';
        case 8: return 'H';
        case 9: return 'I';
        case 10: return 'L';
        case 11: return 'K';
        case 12: return 'M';
        case 13: return 'F';
        case 14: return 'P';
        case 15: return 'S';
        case 16: return 'T';
        case 17: return 'W';
        case 18: return 'Y';
        case 19: return 'V';
        default: System.out.println("Error in getName");
    }
    return 'X';
    }
    
    /*
     * getHydro
     * takes: char
     * returns: double
     * input: (the char value for an amino acid)
     *
     * This function returns the hydrophobisity of an amino acid
     * for a hydrophobe this value is 1
     * for anything else the current value is -1
     * This can be changed on a grouping basis (4 groups are already defined:
     * hydrophobes, hydrophiles, acids, and bases)
     * This has the effect of giveing larger positive numbers to hydrophobic things
     * and larger (in magnitude) positive numbers to hydrophilic things
     */
    private static double getHydro(char r){
        double hydrophobe = 1.0;
        double hydrophile = -1.0;
        double base = -1.0;
        double acid = -1.0;
        switch (r) {
            case 'A': return hydrophobe;
            case 'R': return base;
            case 'N': return hydrophile;
            case 'D': return acid;
            case 'C': return hydrophile;
            case 'Q': return hydrophile;
            case 'E': return acid;
            case 'G': return hydrophile;
            case 'H': return base;
            case 'I': return hydrophobe;
            case 'L': return hydrophobe;
            case 'K': return base;
            case 'M': return hydrophobe;
            case 'F': return hydrophobe;
            case 'P': return hydrophobe;
            case 'S': return hydrophile;
            case 'T': return hydrophile;
            case 'W': return hydrophobe;
            case 'Y': return hydrophile;
            case 'V': return hydrophobe;
            default: System.out.println("Error in getHydro: " + r);
        }
        return 0.0;
    }
    
    /*
     * getHydroType
     * takes char
     * returns int
     * input: (character code for an amino acid)
     *
     * This function returns what hydrophobic group an amino acid belongs to
     * 0. Hydrophobes
     * 1. Hydrophiles
     * 2. bases
     * 3. acids
     * This can be easily expanded or contracted to make more gorups or change groups
     */
    private static int getHydroType(char r){
        int hydrophobe = 0;
        int hydrophile = 1;
        int base = 2;
        int acid = 3;
        switch (r) {
            case 'A': return hydrophobe;
            case 'R': return base;
            case 'N': return hydrophile;
            case 'D': return acid;
            case 'C': return hydrophile;
            case 'Q': return hydrophile;
            case 'E': return acid;
            case 'G': return hydrophile;
            case 'H': return base;
            case 'I': return hydrophobe;
            case 'L': return hydrophobe;
            case 'K': return base;
            case 'M': return hydrophobe;
            case 'F': return hydrophobe;
            case 'P': return hydrophobe;
            case 'S': return hydrophile;
            case 'T': return hydrophile;
            case 'W': return hydrophobe;
            case 'Y': return hydrophile;
            case 'V': return hydrophobe;
            default: System.out.println("Error in getHydroType: " + r);
        }
        return 0;
    }
    
    /*
     * seqres
     * takes String
     * returns char
     * input: (The three letter code for any amino acid)
     *
     * This function basically turns the three letter code for any amino acid into
     * its one letter code and returns it.
     */
    private static char seqres(String res){
        if(res.compareTo("ALA")==0){return 'A';}
        if(res.compareTo("ARG")==0){return 'R';}
        if(res.compareTo("ASN")==0){return 'N';}
        if(res.compareTo("ASP")==0){return 'D';}
        if(res.compareTo("CYS")==0){return 'C';}
        if(res.compareTo("GLN")==0){return 'Q';}
        if(res.compareTo("GLU")==0){return 'E';}
        if(res.compareTo("GLY")==0){return 'G';}
        if(res.compareTo("HIS")==0){return 'H';}
        if(res.compareTo("ILE")==0){return 'I';}
        if(res.compareTo("LEU")==0){return 'L';}
        if(res.compareTo("LYS")==0){return 'K';}
        if(res.compareTo("MET")==0){return 'M';}
        if(res.compareTo("PHE")==0){return 'F';}
        if(res.compareTo("PRO")==0){return 'P';}
        if(res.compareTo("SER")==0){return 'S';}
        if(res.compareTo("THR")==0){return 'T';}
        if(res.compareTo("TRP")==0){return 'W';}
        if(res.compareTo("TYR")==0){return 'Y';}
        if(res.compareTo("VAL")==0){return 'V';}
        System.out.println("Error in seqres: " + res);
        return 'X';
    }
    
    /*
     * getSector
     * takes double, double, double, double, double, double
     * returns int
     * input: (x coord of test point, y coord of test point,
     *         x coord of rotation vector, y coord of rotation vector
     *         x coord of vector orthogonal to rotaion vector, y coord of vector orthogonal to rotation vector)
     *
     * This function returns what sector a test point is in, in relation to the rotation vector
     * the sectors are defined as follows:
     *
     *    \  II | II  / <- Vector orthogonal to rotation vector
     *     \    |    /
     *      \   |   /
     *    I  \  |  /  III
     *        \ | /
     *    _____\|/_____
     *         /|\
     *    I   / | \   III
     *       /  |  \
     *      /   |   \
     *     / IV | IV \ <- Rotation vector
     */
    public static int getSector(double x, double y, double rotx, double roty, double protx, double proty){
        double dot = x*rotx + y*roty;
    double pdot = x*protx + y*proty;
    if(dot > 0){
        if(pdot > 0){
            return 3;
        }else{
            return 4;
        }
    }else{
        if(pdot > 0){
            return 2;
        }else{
            return 1;
        }
    }
    }
    
    /*
     * transPDB
     * takes boolean, boolean, boolean, String
     * returns void (outputs to a file)
     * input: (reflection about x axis?, y axis?, z axis?, input file)
     * the input is defined on/around line 625
     *
     * This function will reflect the *.pdb file about the x, y, and/or z axis and export the new
     * coordinate file to trans.pdb
     *
     * This can be used to artificially force a protien to assume a right handed helix instead of
     * a left handed helix, or in combination with the ramRot function can generate almost entirely
     * new protien configurations from an already stable *.pdb file
     */
    public static void transPDB(boolean rotX, boolean rotY, boolean rotZ, String ins) throws
    IOException,
    FileNotFoundException
    {
        BufferedReader fin = new BufferedReader(new FileReader("" + ins + ""));
        PrintWriter fout = new PrintWriter(new FileOutputStream("trans.pdb"), true);
        int i, len, toNeg = 0;
        double xTemp = 0.0, yTemp = 0.0, zTemp = 0.0;
        boolean stop = false;
        String line = null, dataPRE, dataX, dataY, dataZ, dataEND;
    
        while(stop != true){
            line = fin.readLine();
        stop = (line.substring(0,4).compareTo("ATOM")==0);
        if(stop != true){
            fout.println(line);
        }
        }
    
        stop = false;
    
        while(stop != true){
            dataPRE = line.substring(0,31);
            dataX = line.substring(31,38).trim();
        dataY = line.substring(39,46).trim();
        dataZ = line.substring(47,54).trim();
        dataEND = line.substring(55);
    
        if(rotX){
            len = dataX.length();
            xTemp = Double.parseDouble(dataX);
            if(xTemp > 0){
                toNeg = 1;
            }
            xTemp = (-1.0)*xTemp;
            dataX = Double.toString(xTemp);
            for(i = 0; i < 7; i++){
                if(dataX.length() != (len+toNeg)){
                    dataX = (dataX + "0");
            }
            }
            len = 6 - len;
            for(i = 0; i < len; i++){
                dataX = (" " + dataX);
            }
        }
        dataPRE = dataPRE.concat(dataX);
        
        if(rotY){
            len = dataY.length();
            yTemp = Double.parseDouble(dataY);
            if(yTemp > 0){
                toNeg = 1;
            }
            yTemp = (-1.0)*yTemp;
            dataY = Double.toString(yTemp);
            for(i = 0; i < 7; i++){
                if(dataY.length() != (len+toNeg)){
                    dataY = (dataY + "0");
            }
            }
            len = 6 - len;
            for(i = 0; i < len; i++){
                dataY = (" " + dataY);
            }
        }
        dataPRE = dataPRE.concat(dataY);
    
        if(rotZ){
            len = dataZ.length();
            zTemp = Double.parseDouble(dataZ);
            if(zTemp > 0){
                toNeg = 1;
            }
            zTemp = (-1.0)*zTemp;
            dataZ = Double.toString(zTemp);
            for(i = 0; i < 7; i++){
                if(dataZ.length() != (len+toNeg)){
                    dataZ = (dataZ + "0");
            }
            }
            len = 6 - len;
            for(i = 0; i < len; i++){
                dataZ = (" " + dataZ);
            }
        }
        dataPRE = dataPRE.concat(dataZ);
        
        dataPRE = dataPRE.concat(dataEND);
    
        fout.println(dataPRE);
    
        line = fin.readLine();
        stop = (line.substring(0,4).compareTo("ATOM")!=0);
        }
    }
}

// END of program
