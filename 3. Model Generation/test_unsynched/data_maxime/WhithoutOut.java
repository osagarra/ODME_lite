/*
 *  Copyright (C) 2010 Cemagref
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math.stat.inference.ChiSquareTestImpl;
import static java.lang.Math.*;

/**
 *
 * @author Maxime Lenormand <maxime.lenormand@cemagref.fr>
 */
public class WhithoutOut {

    static String wd = new File(System.getProperty("user.dir")) + File.separator;
    static ChiSquareTestImpl chiSquareTestImpl = new ChiSquareTestImpl();

    public static void main(String[] args) throws FileNotFoundException, IOException {

        // load observed data
        int nbb = 0;
        int NbComb = 0;
        //ID, IN et OUT
        List<int[]> INb = new ArrayList<int[]>();
        List<int[]> OUTb = new ArrayList<int[]>();
        Scanner scanner = new Scanner(new File(wd + "NetworkFRc0.csv"));
        scanner.nextLine();
        while (scanner.hasNextLine()) {
            String[] cols = scanner.nextLine().split(";");
            int[] in = new int[2];
            int[] out = new int[2];
            in[0] = Integer.parseInt(cols[0]);
            in[1] = Integer.parseInt(cols[1]);
            out[0] = Integer.parseInt(cols[0]);
            out[1] = Integer.parseInt(cols[2]);

            INb.add(in);
            OUTb.add(out);
            NbComb++;
            nbb = nbb + Integer.parseInt(cols[1]);
        }

        //Enlever les 0
        int v = 0;
        int u = INb.size();
        int y = 0;
        while (y < u) {
            if (INb.get(v)[1] == 0) {
                INb.remove(v);
            } else {
                v++;
            }
            y++;
        }
        v = 0;
        u = OUTb.size();
        y = 0;
        while (y < u) {
            if (OUTb.get(v)[1] == 0) {
                OUTb.remove(v);
            } else {
                v++;
            }
            y++;
        }

        //Distance
        double[][] Distb = new double[NbComb][NbComb];
        scanner = new Scanner(new File(wd + "DistanceFRc0.txt"));
        scanner.nextLine();
        int k = 0;
        while (scanner.hasNextLine()) {
            String[] cols = scanner.nextLine().split(" ");
            for (int i = 1; i < cols.length; i++) {
                Distb[k][(i - 1)] = Double.parseDouble(cols[i]) * 1000;
            }
            k++;
        }

        File reg = new File("C:/FRc/FRc0");
        reg.mkdirs();

        for (int beta = 0; beta <= 20; beta++) {
            if (beta != 4) {
                File bet = new File("C:/FRc/FRc0/Beta" + (double) (61 + beta * 2));
                bet.mkdirs();
                for (int t = 0; t < 1; t++) {
                    PrintWriter writer = new PrintWriter(new File("C:/FRc/FRc0/Beta" + (double) (61 + beta * 2) + "/CommutingNumber_r" + t + ".dat"));
                    int[][] com = Commuting((double) (61 + beta * 2), NbComb, nbb, INb, OUTb, Distb);
                    for (int i = 0; i < com.length; i++) {
                        for (int j = 0; j < com.length; j++) {
                            writer.print(com[i][j]);
                            writer.print(" ");
                        }
                        writer.println();
                    }
                    writer.close();
                }
            }
        }
    }

    static int[][] Commuting(double beta, int NbComb, int nbb, List<int[]> INb, List<int[]> OUTb, double[][] Dist) throws FileNotFoundException, IOException {
        int NbCom = NbComb;
        int nb = nbb;
        List<int[]> IN = new ArrayList<int[]>();
        List<int[]> OUT = new ArrayList<int[]>();
        for (int j = 0; j < INb.size(); j++) {
            IN.add(new int[2]);
            IN.get(j)[0] = INb.get(j)[0];
            IN.get(j)[1] = INb.get(j)[1];
        }
        for (int j = 0; j < OUTb.size(); j++) {
            OUT.add(new int[2]);
            OUT.get(j)[0] = OUTb.get(j)[0];
            OUT.get(j)[1] = OUTb.get(j)[1];
        }
        //Resultat
        int[][] MatCom = new int[NbCom][NbCom];
        for (int i = 0; i < NbCom; i++) {
            for (int j = 0; j < NbCom; j++) {
                MatCom[i][j] = 0;
            }
        }

        while (nb > 0) {
            int Nout = (int) round((double) (random() * (OUT.size() - 1)));
            OUT.get(Nout)[1]--;
            int Nin = 0;
            List<Double> prob1 = new ArrayList<Double>();
            List<Integer> prob0 = new ArrayList<Integer>();
            double sum = 0;
            for (int t = 0; t < IN.size(); t++) {
                if (IN.get(t)[0] != OUT.get(Nout)[0]) {
                    //sum = sum + (double) ((double) IN.get(t)[1] * pow(Dist[IN.get(t)[0]][OUT.get(Nout)[0]], -beta));
                    sum = sum + (double) ((double) IN.get(t)[1] * exp((Dist[OUT.get(Nout)[0]][IN.get(t)[0]]) * (-beta)));
                }
            }

            int y = 0;
            prob1.add((double) 0);
            for (int t = 0; t < IN.size(); t++) {
                if (IN.get(t)[0] != OUT.get(Nout)[0]) {
                    //prob1.add((double) ((double) IN.get(t)[1] * pow(Dist[IN.get(t)[0]][OUT.get(Nout)[0]], -beta)) / sum + prob1.get(y));
                    prob1.add((double) ((double) IN.get(t)[1] * exp((Dist[OUT.get(Nout)[0]][IN.get(t)[0]]) * (-beta))) / sum + prob1.get(y));
                    prob0.add(t);
                    y++;
                }
            }

            double r = random();
            for (int t = 0; t < (prob1.size() - 1); t++) {
                if (r >= prob1.get(t) & r < prob1.get(t + 1)) {
                    Nin = prob0.get(t);
                }
            }

            IN.get(Nin)[1]--;
            MatCom[OUT.get(Nout)[0]][IN.get(Nin)[0]]++;

            if (OUT.get(Nout)[1] == 0) {
                OUT.remove(Nout);
            }
            if (IN.get(Nin)[1] == 0) {
                IN.remove(Nin);
            }
            nb--;
        }
        return MatCom;
    }
}
