//
// UCESTest.cpp -- implementation of the namespace "UCESTest".
//
//    This file is part of the featsel program
//    Copyright (C) 2010  Marcelo S. Reis
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "UCESTest.h"

namespace UCESTest
{

    //always run nroLocalMinima before get_steps
    //always run nroLocalMinima before get_steps
    //always run nroLocalMinima before get_steps
    //always run nroLocalMinima before get_steps
//
//	bool check_number_of_minima() { //past tests, doesnt work now
//		ElementSet set1 ("set", "input/CostFunctionTest3ElementsFile.xml");
//		ElementSet set2 ("set", "input/UCSTest9ElementsFileA.xml");
//		ElementSet set3 ("set", "input/UCSTest7ElementsFileA.xml");
//		UCES uces1;
//		UCES uces2;
//		UCES uces3;
//		MeanAbsSum c1 (&set1);
//		MeanAbsSum c2 (&set2);
//		HammingDistance c3 (&set3);
//		uces1.set_parameters (&c1, &set1, false);
//		uces1.get_minima_list (1);
////		cout << uces1.print_list_of_minima () << "\n";
////		cout << uces1.print_search_report () << "\n";
////		cout << uces1.print_list_of_visited_subsets () << "\n";
//        uces2.set_parameters (&c2, &set2, false);
//		uces2.get_minima_list (3);
////		cout << uces2.print_list_of_minima () << "\n";
////		cout << uces2.print_search_report () << "\n";
////		cout << uces2.print_list_of_visited_subsets () << "\n";
//		uces3.set_parameters (&c3, &set3, false);
//		uces3.get_minima_list (1);
////		cout << uces3.print_list_of_minima () << "\n";
////		cout << uces3.print_search_report () << "\n";
////		cout << uces3.print_list_of_visited_subsets () << "\n";
//        return true;
//	}
//
//    bool check_results() {
//
//        int n = 25;
//
//        string arr[] = {
//        "input/UCS2bTest5ElementsFileA.xml", 
//        "input/UCS2bTest5ElementsFileB.xml",
//        "input/UCS2bTest6ElementsFileA.xml",
//        "input/UCS2bTest6ElementsFileB.xml",
//        "input/UCS2bTest9ElementsFileA.xml",
//        "input/UCSTest4ElementsFileA.xml",
//        "input/UCSTest4ElementsFileB.xml", 
//        "input/UCSTest5ElementsFileA.xml",
//        "input/UCSTest6ElementsFileA.xml",
//        "input/UCSTest6ElementsFileB.xml",
//        "input/UCSTest7ElementsFileA.xml",
//        "input/UCSTest7ElementsFileB.xml",
//        "input/UCSTest7ElementsFileC.xml",
//        "input/UCSTest7ElementsFileD.xml",
//        "input/UCSTest9ElementsFileA.xml",
//        "input/UCSTest9ElementsFileB.xml",
//        "input/UCSTest10ElementsFileA.xml",
//        "input/UCSTest11ElementsFileA.xml",
//        "input/UCSTest13ElementsFileA.xml",
//        "input/UCSTest13ElementsFileB.xml",
//        "input/UCSTest14ElementsFileA.xml",
//        "input/UCSTest15ElementsFileA.xml",
//        "input/UCSTest17ElementsFileA.xml",
//        "input/UCSTest19ElementsFileA.xml",
//        "input/UCSTest21ElementsFileA.xml"
//        };
//        vector <string> S(arr, arr+n);
//
//        vector <int> positions(5), results;
//        positions[0] = 0;
//        positions[1] = 1;
//        positions[2] = 2;
//        positions[3] = 4;
//        positions[4] = 9;
//
//
//        UCES uces;
//        for (int i = 0; i < n; i++) {
//            cout << S[i] << endl;
//            ElementSet pset("set", S[i]);
//            MeanAbsSum pc (&pset);
//            uces.set_parameters (&pc, &pset, false);
//            int nset = pset.get_set_cardinality ();
//            cout << "size of the set = " << nset << endl;
//            vector <int> Costs;
//            int nMinima = uces.nLocalMinima(Costs);
//            cout << "number of Local Minima = " << nMinima << endl;
//            double media = 0;
//            int iterations = 20, nresults = positions.size();
//            vector <double> mediaResults(nresults, 0.0);
//            for (int j = 0; j < iterations; j++) {
//                media += uces.get_steps(Costs, results, positions, 1000);
//                for (int k = 0; k < nresults; k++) {
//                    mediaResults[k] += results[k];
//                }
//            }
//
//            for (int k = 0; k < nresults; k++) {
//                mediaResults[k] /= iterations;
////                cout << "Media of results position = " << positions[k]+1 << " result = " << mediaResults[k] << endl;
//            }
//
//
//            double mediaSearchingMinima = media*1.0/iterations;
//            cout << "media of searching minima " << mediaSearchingMinima << endl;
//            for (int j = 0; j < nresults; j++) {
//                double muT = nMinima*1.0/(1<<nset);
//                double epsilon = (positions[j]+1)*1.0/(1<<nset);
//                double log1d = log(1.0/0.001); //delta 0.1%, 1-delta = 99.9%
//                double lowerbound = muT/epsilon*log1d;
//                double upperbound = 1.0/epsilon*log1d;
//                cout << "Expected values for position=" << positions[j]+1 << 
//                    ", lower bound= " << lowerbound << 
//                    ", upper bound= " << upperbound << 
//                    ", media real value = " << mediaResults[j] << 
//                    ", media of nodes visited = " << mediaResults[j]*mediaSearchingMinima << endl;
//            }
//        }
//        
//        return true;
//    }

    bool check_results() {

        int n = 9;

        string arr[] = {
        "input/subset_sum/Test_03_A.xml", 
        "input/subset_sum/Test_04_A.xml", 
        "input/subset_sum/Test_04_B.xml", 
        "input/subset_sum/Test_07_A.xml", 
        "input/subset_sum/Test_07_B.xml", 
        "input/subset_sum/Test_07_C.xml", 
        "input/subset_sum/Test_09_A.xml", 
        "input/subset_sum/Test_09_B.xml", 
        "input/subset_sum/Test_13_A.xml" 
        };
        vector <string> S(arr, arr+n);

        vector <int> positions(5), results;
        positions[0] = 0;
        positions[1] = 1;
        positions[2] = 2;
        positions[3] = 4;
        positions[4] = 7;


        UCES uces;
        for (int i = 0; i < n; i++) {
            cout << S[i] << endl;
            ElementSet pset("set", S[i]);
            SubsetSum pc (&pset);
            uces.set_parameters (&pc, &pset, false);
            int nset = pset.get_set_cardinality ();
            cout << "size of the set = " << nset << endl;
            vector <int> Costs;
            int nMinima = uces.nLocalMinima(Costs);
//            for (int j = 0; j < Costs.size(); j++) {
//                cout << " costs " << j << " = " << Costs[j] << endl;
//            }
            cout << "number of Local Minima = " << nMinima << endl;
            double media = 0;
            int iterations = 20, nresults = positions.size();
            vector <double> mediaResults(nresults, 0.0);
            for (int j = 0; j < iterations; j++) {
                media += uces.get_steps(Costs, results, positions, 1000);
                for (int k = 0; k < nresults; k++) {
                    mediaResults[k] += results[k];
//                    cout << "Media of results position = " << positions[k]+1 << " result = " << results[k] << endl;
                }
            }

            for (int k = 0; k < nresults; k++) {
                mediaResults[k] /= iterations;
//                cout << "Media of results position = " << positions[k]+1 << " result = " << mediaResults[k] << endl;
            }


            double mediaSearchingMinima = media*1.0/iterations;
            cout << "media of searching minima " << mediaSearchingMinima << endl;
            for (int j = 0; j < nresults; j++) {
                double muT = nMinima*1.0/(1<<nset);
                double epsilon = (positions[j]+1)*1.0/(1<<nset);
                double log1d = log(1.0/0.001); //delta 0.1%, 1-delta = 99.9%
                double lowerbound = muT/epsilon*log1d;
                double upperbound = 1.0/epsilon*log1d;
                cout << "Expected values for position=" << positions[j]+1 << 
                    ", lower bound= " << lowerbound << 
                    ", upper bound= " << upperbound << 
                    ", media real value = " << mediaResults[j] << 
                    ", media of nodes visited = " << mediaResults[j]*mediaSearchingMinima << endl;
            }
        }
        
        return true;
    }

    bool check_results2() { //works but it is slow

        int n = 10;

        string arr[] = {
        "input/tmp/Test_020_0001.xml",
        "input/tmp/Test_020_0002.xml",
        "input/tmp/Test_020_0003.xml",
        "input/tmp/Test_020_0004.xml",
        "input/tmp/Test_020_0005.xml",
        "input/tmp/Test_020_0006.xml",
        "input/tmp/Test_020_0007.xml",
        "input/tmp/Test_020_0008.xml",
        "input/tmp/Test_020_0009.xml",
        "input/tmp/Test_020_0010.xml"
        };
        vector <string> S(arr, arr+n);

        vector <double> epsilon(5), delta(5), results;
        epsilon[0] = 0.00001;
        epsilon[1] = 0.00002;
        epsilon[2] = 0.00005;
        epsilon[3] = 0.0001;
        epsilon[4] = 0.001;
        delta[0] = 0.1;
        delta[1] = 0.05;
        delta[2] = 0.01;
        delta[3] = 0.005;
        delta[4] = 0.001;

        UCES uces;
        for (int i = 0; i < n; i++) {
            cout << S[i] << endl;
            ElementSet pset("set", S[i]);
            SubsetSum pc (&pset);
            uces.set_parameters (&pc, &pset, false);
            int nset = pset.get_set_cardinality ();
            cout << "size of the set = " << nset << endl;
            vector <int> Costs;
//            int nMinima = uces.nLocalMinima(Costs); //returns number of Minima
//            cout << nMinima << " p real " << nMinima*1.0/pow(2.0, nset) << endl;
            uces.nCosts(Costs); //returns only costs
//            for (int j = 0; j < 300; j++) {
//                cout << " costs " << j << " = " << Costs[j] << endl;
//            }
            //cout << "number of Local Minima = " << nMinima << endl;
            int iterations = 20, nepsilon = epsilon.size(), ndelta = delta.size();
            for (int j = 0; j < iterations; j++) {
                for (int e = 0; e < nepsilon; e++) {
                    for (int d = 0; d < ndelta; d++) {
                        int position = max(pow(2.0, nset)*epsilon[e]-1, 0.0);
                        double valueFound = uces.get_minima_eps_delta(1, epsilon[e], delta[e]);
                        cout << "value found = " <<  valueFound << " position = " << position << " Cost[position] " << Costs[position] << endl;
                        if (valueFound <= Costs[position]) {
                            cout << "true" << endl;
                        }
                        else {
                            cout << "false" << endl;
                        }
                    }
                }
            }
        }
        
        return true;
    }

    bool check_results3() { //works faster than check_results2, epsilon[0] and delta[0] minimum values!!

        int n = 10;

        string arr[] = {
        "input/tmp/Test_020_0001.xml",
        "input/tmp/Test_020_0002.xml",
        "input/tmp/Test_020_0003.xml",
        "input/tmp/Test_020_0004.xml",
        "input/tmp/Test_020_0005.xml",
        "input/tmp/Test_020_0006.xml",
        "input/tmp/Test_020_0007.xml",
        "input/tmp/Test_020_0008.xml",
        "input/tmp/Test_020_0009.xml",
        "input/tmp/Test_020_0010.xml"
        };
        vector <string> S(arr, arr+n);

        vector <double> epsilon(5), delta(5); // epsilon[0] and delta[0] minimum values!!
        epsilon[0] = 0.00001;
        epsilon[1] = 0.00002;
        epsilon[2] = 0.00005;
        epsilon[3] = 0.0001;
        epsilon[4] = 0.001;
        delta[0] = 0.001;
        delta[1] = 0.005;
        delta[2] = 0.01;
        delta[3] = 0.05;
        delta[4] = 0.1;

        int nepsilon = epsilon.size(), ndelta = delta.size();


        UCES uces;
        for (int i = 0; i < n; i++) {
            cout << S[i] << endl;
            ElementSet pset("set", S[i]);
            SubsetSum pc (&pset);
            uces.set_parameters (&pc, &pset, false);
            int nset = pset.get_set_cardinality ();
            cout << "size of the set = " << nset << endl;
            vector <int> Costs;
//            int nMinima = uces.nLocalMinima(Costs); //returns number of Minima
//            cout << "number of Local Minima = " << nMinima << " p real " << nMinima*1.0/pow(2.0, nset) << endl;
            uces.nCosts(Costs); //returns only costs
//            for (int j = 0; j < 300; j++) {
//                cout << " costs " << j << " = " << Costs[j] << endl;
//            }
            //cout << "number of Local Minima = " << nMinima << endl;
            vector < vector <int> > Results( nepsilon, vector <int>(ndelta, 0));
            int TotalIterations = 100;
            for (int j = 0; j < TotalIterations; j++) {
                cout << j << endl;
                vector <double> values;
                double pHat;
                values = uces.get_minima_list_eps_delta(1, epsilon[0], delta[0], pHat);
                
                //iterationsPhat must be the same that UCES.cpp
//                int iterationsPhat = 13572; //phat error = 0.01 with 99% of confidence, p <= p
//                if (pHat*iterationsPhat < 10 ||  (1-pHat)*iterationsPhat < 10) {
//                    cout << "error estimation of p hat, np < 10 or n(1 − p) < 10" << endl;
//                }

                for (int e = 0; e < nepsilon; e++) {
                    for (int d = 0; d < ndelta; d++) {

                        int iterations = pHat/epsilon[e]*log(1.0/delta[d]);

                        int position = max(pow(2.0, nset)*epsilon[e]-1, 0.0);
                        double valueFound = values[iterations-1];
//                        cout << "epsilon = " << epsilon[e] << " delta = " << delta[d] << " value found = " <<  valueFound << " position = " << position << " Cost[position] " << Costs[position] << endl;
                        if (valueFound <= Costs[position]) {
                            Results[e][d]++;
//                            cout << "true" << endl;
                        }
                        else {
//                            cout << "false" << endl;
                        }
                    }
                }
            }

            for (int e = 0; e < nepsilon; e++) {
                for (int d = 0; d < ndelta; d++) {
                    cout << Results[e][d] << " ";
                }
                cout << endl;
            }
        }


        
        return true;
    }


} // end of namespace
