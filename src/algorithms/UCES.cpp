//
// ExhaustiveSearch.cpp -- implementation of the class "ExhaustiveSearch".
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

#include "UCES.h"


UCES::UCES ()
{
    list_of_visited_subsets = new Collection ();
	cost_function = NULL;
}


UCES::~UCES ()
{
  if (list_of_visited_subsets != NULL)
    delete list_of_visited_subsets;
}

// get_minima_list return the local minima, and we run this algorithm n times
// this returns the number of steps to get to a local minima
int UCES::dfs (ElementSubset * X) {
//    cout << "start dfs" << endl;
    X->cost = cost_function->cost (X);  
//    cout << "X cost = " << X->cost << endl;
    
	int n = set->get_set_cardinality ();

	ElementSubset * Y; 
    Y = new ElementSubset("", set);

    bool Neighbors[n];
    int cntneighbors = 0;
    for (int j = 0; j < n; j++) {
        Neighbors[j] = false;
        Y = new ElementSubset ("", set);
        Y->copy(X);
        if (X->has_element(j)) 
            Y->remove_element(j);
        else 
            Y->add_element(j);


        Y->cost = cost_function->cost (Y);
//        cout << Y->cost << " " << X->cost << endl;
//        cout << "Y cost = " << Y->cost << endl;
        if (Y->cost < X->cost) {
            cntneighbors++;
            Neighbors[j] = true;
        }
        delete Y;
    }

//    cout << cntneighbors << endl;
    if (cntneighbors == 0) return 1;

    int aleatory_choose = rand()%cntneighbors;
    int vis = 0;
//    cout << "watching neighbors" << endl;
    for (int j = 0; j < n; j++) if (Neighbors[j]) {
//        cout << j << endl;
        if (aleatory_choose == vis) {
//            cout << j << " " << aleatory_choose << " " << vis << endl;
            if (X->has_element(j)) 
                X->remove_element(j);
            else 
                X->add_element(j);
            break;
            //X = Y
        }
        vis++;
    }
////    cout << "end dfs" << endl;
    return 1 + dfs(X);
}

// get_minima_list return the local minima, and we run this algorithm n times
// this returns the number of steps to get to a local minima
// same as dfs but it only chooses minimum cost, doesnt 
int UCES::dfs2 (ElementSubset * X) {
//    cout << "start dfs" << endl;
    X->cost = cost_function->cost (X);  
//    cout << "X cost = " << X->cost << endl;
    
	int n = set->get_set_cardinality ();

	ElementSubset * Y; 
    Y = new ElementSubset("", set);

    bool Neighbors[n];
    int cntneighbors = 0;
    double current = X->cost;
    int curNeighbor = -1;
    for (int j = 0; j < n; j++) {
        Neighbors[j] = false;
        Y = new ElementSubset ("", set);
        Y->copy(X);
        if (X->has_element(j)) 
            Y->remove_element(j);
        else 
            Y->add_element(j);


        Y->cost = cost_function->cost (Y);
//        cout << Y->cost << " " << X->cost << endl;
//        cout << "Y cost = " << Y->cost << endl;
        if (Y->cost <= current) {
            if (curNeighbor != -1) Neighbors[curNeighbor] = false;
            current = Y->cost;
            Neighbors[j] = true;
            cntneighbors = 1;
        }
        delete Y;
    }

//    cout << cntneighbors << endl;
    if (cntneighbors == 0) return 1;

    int aleatory_choose = rand()%cntneighbors;
    int vis = 0;
//    cout << "watching neighbors" << endl;
    for (int j = 0; j < n; j++) if (Neighbors[j]) {
//        cout << j << endl;
        if (aleatory_choose == vis) {
//            cout << j << " " << aleatory_choose << " " << vis << endl;
            if (X->has_element(j)) 
                X->remove_element(j);
            else 
                X->add_element(j);
            break;
            //X = Y
        }
        vis++;
    }
////    cout << "end dfs" << endl;
    return 1 + dfs2(X);
}

//returns the minimum element for that epsilon and delta
//
double UCES::get_minima_eps_delta (unsigned int max_size_of_minima_list, double epsilon, double delta)
{
	timeval begin_program, end_program;
	gettimeofday (& begin_program, NULL);

	ElementSubset * X, * Y;
	int n = set->get_set_cardinality ();

    int cnt_iterations = 0;

    //this part is modified to use a calculated m
    //

    int iterationsPhat = 13572; //phat error = 0.01 with 99% of confidence, p <= p
    int iterations = iterationsPhat + 5;

    double curCost = 1e120;
    //
    // end of the part

    Y = new ElementSubset ("", set);

    vector <double> Results;

    int c_phat = 0;
    do {

        if (cnt_iterations == iterationsPhat) { //calculate real number of iterations to satisfy epsilon and delta
            double pHat = c_phat*1.0/iterationsPhat + 0.01; //add the error
            if (pHat*iterationsPhat < 10 ||  (1-pHat)*iterationsPhat < 10) {
                cout << "error estimation of p hat, np < 10 or n(1 − p) < 10" << endl;
            }
            iterations = pHat/epsilon*log(1.0/delta);
            cout << "phat = " << pHat << " and iterations = " << iterations << endl;
        }

        X = new ElementSubset("X", set);
        X->set_empty_subset (); // X starts with empty set

        int aleatory_subset = rand()%(1<<n);
        for (int i = 0; i < n; i++) {
            if (aleatory_subset & (1<<i))
                X->add_element(i);
        }

//        X->cost = cost_function->cost (X);  
//        cout << "start = " << aleatory_subset << " " <<  X->cost << endl;

        int sz = dfs(X);
        if (sz == 1) {
            c_phat++;
        }

//        cout << "end = " << X->cost << endl;

        //temporal improve time (only works for best solution!!)
        if (X->cost < curCost) {
            Y->copy(X);
            curCost = X->cost;
        }
        Results.push_back(curCost);
        //temporal improve time

        //list_of_minima.push_back (X);

        cnt_iterations++;
    } while (cnt_iterations < iterations);

    list_of_minima.push_back (Y);
//    for (list<ElementSubset *>::iterator it = list_of_minima.begin(); it != list_of_minima.end(); it++) {
//        cout << (*it)->cost << endl;
//    }

//    sort(list_of_minima.begin(), list_of_minima.end());
//    list_of_minima.resize(1);

	clean_list_of_minima (max_size_of_minima_list);

	// Exhaustive search, if implemented keeping just an element of minimum cost,
	// needs to compute the cost function 2^|S| times.
	//
	number_of_visited_subsets =  cost_function->get_number_of_calls_of_cost_function ();


	gettimeofday (& end_program, NULL);
	elapsed_time_of_the_algorithm = diff_us (end_program, begin_program);

    return Results[iterations-1];

}

//returns the list of minimum element for that epsilon and delta//
//
vector <double> UCES::get_minima_list_eps_delta (unsigned int max_size_of_minima_list, double epsilon, double delta, double &pHat)
{
	timeval begin_program, end_program;
	gettimeofday (& begin_program, NULL);

	ElementSubset * X, * Y;
	int n = set->get_set_cardinality ();

    int cnt_iterations = 0;

    //this part is modified to use a calculated m
    //

    int iterationsPhat = 13572; //phat error = 0.01 with 99% of confidence, p <= p
    int iterations = iterationsPhat + 5;

    double curCost = 1e120;
    //
    // end of the part

    Y = new ElementSubset ("", set);

    vector <double> Results;

    int c_phat = 0;
    do {

        if (cnt_iterations == iterationsPhat) { //calculate real number of iterations to satisfy epsilon and delta
//            cout << "counter = " << c_phat << " iterationsPhat = " << iterationsPhat << endl;
            pHat = c_phat*1.0/iterationsPhat; //add the error
            if (pHat*iterationsPhat < 10 ||  (1-pHat)*iterationsPhat < 10) {
                cout << "error estimation of p hat, np < 10 or n(1 − p) < 10" << endl;
                pHat = 10.0/iterationsPhat;
            }
            iterations = pHat/epsilon*log(1.0/delta);
            cout << "phat = " << pHat << " and iterations = " << iterations << endl;
        }

        X = new ElementSubset("X", set);
        X->set_empty_subset (); // X starts with empty set

        int aleatory_subset = rand()%(1<<n);
        for (int i = 0; i < n; i++) {
            if (aleatory_subset & (1<<i))
                X->add_element(i);
        }

//        X->cost = cost_function->cost (X);  
//        cout << "start = " << aleatory_subset << " " <<  X->cost << endl;

        int sz = dfs(X);
        if (sz == 1) {
            c_phat++;
        }

//        cout << "end = " << X->cost << endl;

        //temporal improve time (only works for best solution!!)
        if (X->cost < curCost) {
            Y->copy(X);
            curCost = X->cost;
        }
        Results.push_back(curCost);
        //temporal improve time

        //list_of_minima.push_back (X);

        cnt_iterations++;
    } while (cnt_iterations < iterations);

    list_of_minima.push_back (Y);
//    for (list<ElementSubset *>::iterator it = list_of_minima.begin(); it != list_of_minima.end(); it++) {
//        cout << (*it)->cost << endl;
//    }

//    sort(list_of_minima.begin(), list_of_minima.end());
//    list_of_minima.resize(1);

	clean_list_of_minima (max_size_of_minima_list);

	// Exhaustive search, if implemented keeping just an element of minimum cost,
	// needs to compute the cost function 2^|S| times.
	//
	number_of_visited_subsets =  cost_function->get_number_of_calls_of_cost_function ();


	gettimeofday (& end_program, NULL);
	elapsed_time_of_the_algorithm = diff_us (end_program, begin_program);

    return Results;

}



//returns the number of local minima specified using lower bound for local minima number of elements
//
void UCES::get_minima_list (unsigned int max_size_of_minima_list)
{
	timeval begin_program, end_program;
	gettimeofday (& begin_program, NULL);

	ElementSubset * X, * Y;
	int n = set->get_set_cardinality ();

    int cnt_size_of_minima_list = 0;

    //this part is modified to use a calculated m
    //

    //epsilon = 1/2^n, delta = 99.9%
    int values[24] = { 7, 8, 10, 13, 18, 24, 32, 43, 59, 81, 111, 153, 212, 294, 407, 637, 788, 1097, 1527, 2128, 2966, 4137, 5770, 8055};
    int iterations = values[n-1];

    double curCost = 1e120;
    //
    // end of the part

    Y = new ElementSubset ("", set);

    do {

        X = new ElementSubset("X", set);
        X->set_empty_subset (); // X starts with empty set

        int aleatory_subset = rand()%(1<<n);
        for (int i = 0; i < n; i++) {
            if (aleatory_subset & (1<<i))
                X->add_element(i);
        }

//        X->cost = cost_function->cost (X);  
//        cout << "start = " << aleatory_subset << " " <<  X->cost << endl;

        dfs(X);

//        cout << "end = " << X->cost << endl;

        //temporal improve time (only works for best solution!!)
        if (X->cost < curCost) {
            Y->copy(X);
            curCost = X->cost;
        }
        //temporal improve time

        //list_of_minima.push_back (X);

        cnt_size_of_minima_list++;
    } while (cnt_size_of_minima_list < iterations);



    list_of_minima.push_back (Y);

//    for (list<ElementSubset *>::iterator it = list_of_minima.begin(); it != list_of_minima.end(); it++) {
//        cout << (*it)->cost << endl;
//    }

//    sort(list_of_minima.begin(), list_of_minima.end());
//    list_of_minima.resize(1);

	clean_list_of_minima (max_size_of_minima_list);

	// Exhaustive search, if implemented keeping just an element of minimum cost,
	// needs to compute the cost function 2^|S| times.
	//
	number_of_visited_subsets =  cost_function->get_number_of_calls_of_cost_function ();


	gettimeofday (& end_program, NULL);
	elapsed_time_of_the_algorithm = diff_us (end_program, begin_program);

}

//returns the median of steps needed to get to the local minima, with the positions given by the user
double UCES::get_steps (vector <int> &Costs, vector <int> &results, vector <int> &positions, unsigned int max_size_of_minima_list) {

	timeval begin_program, end_program;
	gettimeofday (& begin_program, NULL);

	int n = (int) set->get_set_cardinality ();

    sort(positions.begin(), positions.end());

//    for (i = 0; i < positions.size(); i++) {
//        cout << positions[i] << " " << Costs[i] << "\n";
//    }

    //initiate results
    int max_value = (1<<30); //2^30
    results = vector <int> (positions.size(), (max_value)); //high value

	ElementSubset * X;

    int cnt_size_of_minima_list = 0;

//    for (int k = 0; k < 100; k++) {
//        cout << rand()%10000 << endl;
//    }

    int aleatory_subset, i, sum_steps_needed = 0;

    for (cnt_size_of_minima_list = 1; cnt_size_of_minima_list <= max_size_of_minima_list; cnt_size_of_minima_list++) {

        X = new ElementSubset("X", set);
        X->set_empty_subset (); // X starts with empty set


        aleatory_subset = rand()%(1<<n);
//        cout << cnt_size_of_minima_list << " " << aleatory_subset << endl;

//        cout << "aleatory " << rand()%1024 << " " << (1<<n) << " " << aleatory_subset << "\n";
        for (i = 0; i < n; i++) {
            if (aleatory_subset & (1<<i))
                X->add_element(i);
        }

        sum_steps_needed += dfs(X);

        list_of_minima.push_back (X);
        
        for (i = 0; i < positions.size(); i++) {
//            cout << i << " " << X->cost << " " << Costs.size() << endl;
            if (X->cost <= Costs[positions[i]]) {
                results[i] = min(cnt_size_of_minima_list, results[i]);
            }
        }

//        cout << X->cost << endl;

        if (results[0] != max_value) break;

    } 

//    if (cnt_size_of_minima_list != max_size_of_minima_list)
//        cnt_size_of_minima_list++;


//    for (i = 0; i < results.size(); i++) {
//        cout << "position i = " << positions[i] << " result = " << results[i] << " cost = " << Costs[i] << endl;
//    }

	clean_list_of_minima (max_size_of_minima_list);

	// Exhaustive search, if implemented keeping just an element of minimum cost,
	// needs to compute the cost function 2^|S| times.
	//
	number_of_visited_subsets =  cost_function->get_number_of_calls_of_cost_function ();

	gettimeofday (& end_program, NULL);
	elapsed_time_of_the_algorithm = diff_us (end_program, begin_program);

    return sum_steps_needed * 1.0/cnt_size_of_minima_list;

}

//returns the number of local minima in the lattice and all the costs sorted
int UCES::nLocalMinima (vector <int> &Costs) {
//    cout << "start dfs" << endl;
    // find all the costs and then sort them
	ElementSubset X ("X", set);
	int i;
	int n = (int) set->get_set_cardinality ();

	X.set_empty_subset (); // X starts with empty set
        
    ElementSubset * Y; 

    int cntMinima = 0;
    int contador = 0;
	do  // Amortized time per iteration is O(1) + O(f(n))
	{
		i = 0;
        contador++;
		while ((i < n) && (X.has_element (i)))
		{
			X.remove_element (i);
			i++;
		}
		if (i < n)
			X.add_element (i);

		if (store_visited_subsets)
			list_of_visited_subsets->add_subset (&X);

        X.cost = cost_function->cost (&X);  
        Costs.push_back(cost_function->cost (&X));

        Y = new ElementSubset("", set);

        int cntneighbors = 0;
        for (int j = 0; j < n; j++) {
            Y = new ElementSubset ("", set);
            Y->copy(&X);
            if (X.has_element(j)) 
                Y->remove_element(j);
            else 
                Y->add_element(j);


            Y->cost = cost_function->cost (Y);
    //        cout << Y->cost << " " << X->cost << endl;
    //        cout << "Y cost = " << Y->cost << endl;
            if (Y->cost >= X.cost) {
                cntneighbors++;
            }
            delete Y;
        }
        if (cntneighbors == n) {
//            cout << "minima = " << X.cost << endl;
            cntMinima++;
        }

//        cout << X.print_subset () <<  "cost: " << X.cost << endl;
        if (__builtin_popcount(contador) == 1) cout << contador << endl;
//        cout << i << endl;
        //add cost

	}
	while ( (i < n) );
//    cout << cnt << endl;
//

    sort(Costs.begin(), Costs.end());
//    Costs.resize(20); //modified only 20 best costs (take care)
//    cout << "X cost = " << X->cost << endl;
    return cntMinima;
}

//returns the number of local minima in the lattice and all the costs sorted
int UCES::nLocalMinimaMore (vector < pair<int, int> > &Costs) { //second = 1, false, =0, true
//    cout << "start dfs" << endl;
    // find all the costs and then sort them
	ElementSubset X ("X", set);
	int i;
	int n = (int) set->get_set_cardinality ();

	X.set_empty_subset (); // X starts with empty set
        
    ElementSubset * Y; 

    int cntMinima = 0;
    int contador = 0;
	do  // Amortized time per iteration is O(1) + O(f(n))
	{
		i = 0;
        contador++;
		while ((i < n) && (X.has_element (i)))
		{
			X.remove_element (i);
			i++;
		}
		if (i < n)
			X.add_element (i);

		if (store_visited_subsets)
			list_of_visited_subsets->add_subset (&X);

        X.cost = cost_function->cost (&X);  
        Costs.push_back(make_pair(cost_function->cost (&X), 1));

        Y = new ElementSubset("", set);

        int cntneighbors = 0;
        for (int j = 0; j < n; j++) {
            Y = new ElementSubset ("", set);
            Y->copy(&X);
            if (X.has_element(j)) 
                Y->remove_element(j);
            else 
                Y->add_element(j);


            Y->cost = cost_function->cost (Y);
    //        cout << Y->cost << " " << X->cost << endl;
    //        cout << "Y cost = " << Y->cost << endl;
            if (Y->cost >= X.cost) {
                cntneighbors++;
            }
            delete Y;
        }
        if (cntneighbors == n) {
//            cout << "minima = " << X.cost << endl;
            cntMinima++;
            Costs[Costs.size()-1].second = 0;
        }

//        cout << X.print_subset () <<  "cost: " << X.cost << endl;
        if (__builtin_popcount(contador) == 1) cout << contador << endl;
//        cout << i << endl;
        //add cost

	}
	while ( (i < n) );
//    cout << cnt << endl;
//

    sort(Costs.begin(), Costs.end());
//    Costs.resize(20); //modified only 20 best costs (take care)
//    cout << "X cost = " << X->cost << endl;
    return cntMinima;
}

//returns the number of local minima in the lattice and all the costs sorted
bool UCES::nCosts (vector <int> &Costs) {
//    cout << "start dfs" << endl;
    // find all the costs and then sort them
	ElementSubset X ("X", set);
	int i;
	int n = (int) set->get_set_cardinality ();

	X.set_empty_subset (); // X starts with empty set

	do  // Amortized time per iteration is O(1) + O(f(n))
	{
		i = 0;
		while ((i < n) && (X.has_element (i)))
		{
			X.remove_element (i);
			i++;
		}
		if (i < n)
			X.add_element (i);

		if (store_visited_subsets)
			list_of_visited_subsets->add_subset (&X);

        X.cost = cost_function->cost (&X);  
        Costs.push_back(cost_function->cost (&X));

	}
	while ( (i < n) );
//    cout << cnt << endl;
//

    sort(Costs.begin(), Costs.end());
//    cout << "X cost = " << X->cost << endl;
    return true;
}

