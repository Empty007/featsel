//
// UcurveBranchandBoundTest.h -- definition of the namespace
// "UcurveBranchandBoundTest".
//
//    This file is part of the featsel program
//    Copyright (C) 2015  Marcelo S. Reis
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

#ifndef UCURVE_BRANCH_AND_BOUND_TEST_H_
#define UCURVE_BRANCH_AND_BOUND_TEST_H_

#include "../../src/algorithms/UcurveBranchandBound.h"
#include "../../src/functions/SubsetSum.h"
#include "../../src/functions/HammingDistance.h"
#include "../../src/functions/Explicit.h"

namespace UcurveBranchandBoundTest
{

	bool a_lattice_with_constant_cost_should_be_totally_exausted ();

	bool it_should_find_a_global_minimum ();

	bool it_should_converge_for_large_hamming_sets ();

	bool it_should_converge_for_large_subset_sum_sets ();

	bool it_should_store_all_the_visited_subsets ();

	bool it_should_give_the_number_of_the_visited_subsets ();

	bool it_should_always_give_the_correct_answer ();

}

#endif /* UCURVE_BRANCH_AND_BOUND_TEST_H_ */
