import java.io.*; // needed for input and output

public class pdbConvert
{
    public static void main(String args[]) throws
    IOException,
    FileNotFoundException
    {
        
        // PUT THE COORDINATE FILE TO BE CONVERTED HERE
        String ins = "C3.pdb"; // <----------------
        
        BufferedReader fin = new BufferedReader(new FileReader("" + ins + ""));
        BufferedReader fin2 = new BufferedReader(new FileReader("" + ins + ""));
        BufferedReader fin3 = new BufferedReader(new FileReader("" + ins + ""));
        PrintWriter fout = new PrintWriter(new FileOutputStream("f_" + ins), true);
        
        String amino = null;
        String line = null;
        int res = 0;
        int resT = 0;
        int lines = 0;
        int i = 0;
        int start = 1;
        boolean stop = false;
        boolean stop2 = false;
        int end = 2;
        boolean CR = true;
        int totallines = 0;
        int toATOM = 0;
        int toTER = 0;
        
        fout.println("Start of fixed f_*.pdb file");
        
        while( stop != true ){
            line = fin2.readLine();
            stop = (line.substring(0,3).compareTo("END")==0);
            if(stop == false){
                if(line.substring(0,3).compareTo("TER")==0){
                    stop = true;
                    end = Integer.parseInt(line.substring(23,26).trim());
                }
            }
            if(stop == false){
                end = Integer.parseInt(line.substring(23,26).trim());
            }
            toTER++;
        }
        stop = false;
        //end = Integer.parseInt(line.substring(23,26).trim());
        //System.out.println("ERROR TEST: end: " + end);
        
        while( stop != true ){
            line = fin.readLine();
            stop = (line.substring(0,4).compareTo("ATOM")==0);
            toATOM++;
        }
        stop = false;
        
        totallines = toTER - toATOM;
        
        //start = Integer.parseInt(line.substring(23,26).trim());
        
        while( stop != true ){
            if(CR == true){
                CR = false;
                lines++;
                if(lines < 10){
                    if(end < 100){
                        fout.print("SEQRES   " + lines + " A  " + end + "  ");
                    }else{
                        fout.print("SEQRES   " + lines + " A  " + end + " ");
                    }
                }else{
                    if(end < 100){
                        fout.print("SEQRES  " + lines + " A  " + end + "  ");
                    }else{
                        fout.print("SEQRES  " + lines + " A  " + end + " ");
                    }
                }
            }
            amino = line.substring(17,20);
            res = Integer.parseInt(line.substring(23,26).trim());
            //System.out.println("ERROR TEST: amino: " + amino + ", res: " + res);
            //System.out.println("ERROR TEST: start: " + start);
            //stop = true;
            if(res > start){
                for(i = 0; i < (res-start); i++){
                    //System.out.println("ERROR TEST: in GLY loop");
                    if(i == 13){
                        //System.out.println("ERROR TEST: in GLY loop2");
                        fout.println(" GLY");
                        start = i + start;
                        i = res - start;
                        CR = true;
                    }else{
                        //System.out.println("ERROR TEST: in GLY loop3");
                        fout.print(" GLY" + "");
                    }
                }
            }
            if(CR == false){
                start = end;
                if(res%13==0){
                    fout.println(" " + amino);
                    CR = true;
                }else{
                    fout.print(" " + amino);
                }
                if(res == end){
                    stop = true;
                    stop2 = true;
                }
                while(stop2 != true){
                    line = fin.readLine();
                    resT = Integer.parseInt(line.substring(23,26).trim());
                    if(resT != res){
                        stop2 = true;
                    }
                }
                stop2 = false;
            }
        }
        
        fout.println("\nBuffer line in f_*.pdb");
        
        for(i = 0; i < toATOM-1; i++){
            line = fin3.readLine();
        }
        for(i = 0; i < totallines; i++){
            line = fin3.readLine();
            fout.println(line);
        }
        fout.println("END OF f_*.pdb");
        
        fout.close();
    }
}