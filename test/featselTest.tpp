//============================================================================
//
//    featselTest.cpp -- Unit tests of the featsel classes.
//
//    This file is part of the featsel program
//    Copyright (C) 2016  Marcelo S. Reis
//
//
//   If you use featsel in your publication, we kindly ask you to acknowledge us
//   by citing the paper that describes this framework:
//
//   M.S. Reis, G. Estrela, C.E. Ferreira and J. Barrera
//   "featsel: A Framework for Benchmarking of
//   Feature Selection Algorithms and Cost Functions"
//   https://github.com/msreis/featsel
//
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
//============================================================================


// Special data structures
//
#include "ElementTest.h"
#include "ElementSetTest.h"
#include "ElementSubsetTest.h"
#include "CollectionTest.h"

// Cost (objective) functions
//
// <COST FUNCTION TEMPLATE 1>

// Solvers (algorithms)
//
// <ALGORITHM TEMPLATE 1>

// Number of passed and failed tests
//
unsigned int number_of_passed_tests = 0;
unsigned int number_of_failed_tests = 0;

string current_class;

void result (string test_name, bool passed)
{
  cout << current_class << "::" << test_name;
  cout.flush ();

  if (passed)
  {
    cout << " OK\n";
    number_of_passed_tests++;
  }
  else
  {
    cout << " FAILED\n";
    number_of_failed_tests++;
  }
}

int main (void)
{
  cout << endl << "Starting Unit Tests... " << endl << endl;

  // Testing Class "Element"
  //
  current_class = "ElementTest";
  result ("an_element_should_have_a_name",
    ElementTest::an_element_should_have_a_name ());
  result ("a_new_element_should_not_have_added_values",
    ElementTest::a_new_element_should_not_have_added_values ());
  result ("an_element_without_all_values_added_should_allow_adding",
    ElementTest::an_element_without_all_values_added_should_allow_adding ());
  result ("an_element_with_all_values_added_should_not_allow_adding",
    ElementTest::an_element_with_all_values_added_should_not_allow_adding ());
  result ("an_added_value_should_be_found",
    ElementTest::an_added_value_should_be_found ());
  cout << endl;

  // Testing Class "ElementSet"
  //
  current_class = "ElementSetTest";
  result ("a_set_should_have_a_name",
    ElementSetTest::a_set_should_have_a_name ());
  result ("an_empty_set_should_not_have_elements",
    ElementSetTest::an_empty_set_should_not_have_elements ());
  result ("it_should_load_data_of_a_set_from_file",
    ElementSetTest::it_should_load_data_of_a_set_from_file ());
  result ("it_should_load_data_of_a_set_from_a_DAT_file",
    ElementSetTest::it_should_load_data_of_a_set_from_a_DAT_file ());
  result ("it_should_create_a_set_with_random_values",
    ElementSetTest::it_should_create_a_set_with_random_values ());
  result ("it_should_return_the_set_cardinality",
    ElementSetTest::it_should_return_the_set_cardinality ());
  result ("it_should_return_an_element_that_belongs_to_the_set",
    ElementSetTest::it_should_return_an_element_that_belongs_to_the_set ());
  result ("it_should_not_return_an_element_that_not_belongs_to_the_set",
    ElementSetTest::
    it_should_not_return_an_element_that_not_belongs_to_the_set ());
  result ("values_from_a_random_set_should_be_within_the_given_range",
    ElementSetTest::
    values_from_a_random_set_should_be_within_the_given_range ());
  cout << endl;

  // Testing Class "ElementSubset"
  //
  current_class = "ElementSubsetTest";
  result ("a_new_subset_should_be_an_empty_set",
    ElementSubsetTest::a_new_subset_should_be_an_empty_set ());
  result ("an_element_not_in_subset_should_be_added",
    ElementSubsetTest::an_element_not_in_subset_should_be_added ());
  result ("an_element_in_subset_should_be_removed",
    ElementSubsetTest::an_element_in_subset_should_be_removed ());
  result ("it_should_give_the_set_that_belongs_the_subset",
    ElementSubsetTest::it_should_give_the_set_that_belongs_the_subset ());
  result ("a_set_should_contains_its_subset",
    ElementSubsetTest::a_set_should_contains_its_subset ());
  result ("a_subset_should_be_contained_by_its_set",
    ElementSubsetTest::a_subset_should_be_contained_by_its_set ());
  result ("a_subset_should_be_successfully_cloned",
    ElementSubsetTest::a_subset_should_be_successfully_cloned ());
  result ("a_random_element_should_be_removed",
    ElementSubsetTest::a_random_element_should_be_removed ());
  result ("it_should_give_the_complement_of_the_set",
    ElementSubsetTest::it_should_give_the_complement_of_the_set ());
  cout << endl;

  // Testing Class "Collection"
  //
  current_class = "CollectionTest";
  result ("a_lower_covered_subset_should_be_found",
    CollectionTest::a_lower_covered_subset_should_be_found ());
  result ("a_non_lower_covered_subset_should_not_be_found",
    CollectionTest::a_non_lower_covered_subset_should_not_be_found ());
  result ("it_should_remove_lower_covered_subsets_in_a_collection",
    CollectionTest::it_should_remove_lower_covered_subsets_in_a_collection ());
  result ("the_evaluated_subset_should_not_be_deleted",
    CollectionTest::the_evaluated_subset_should_not_be_deleted ());
  result ("an_upper_covered_subset_should_be_found",
    CollectionTest::an_upper_covered_subset_should_be_found ());
  result ("a_non_upper_covered_subset_should_not_be_found",
    CollectionTest::a_non_upper_covered_subset_should_not_be_found ());
  result ("it_should_remove_upper_covered_subsets_in_a_collection",
    CollectionTest::it_should_remove_upper_covered_subsets_in_a_collection ());
  result ("it_should_copy_a_collection",
    CollectionTest::it_should_copy_a_collection ());
  cout << endl;

  // <COST FUNCTION TEMPLATE 2>

  // <ALGORITHM TEMPLATE 2>

  // Summary of the executed tests
  //
  cout << "Total " << number_of_passed_tests + number_of_failed_tests <<
  " test(s), ";
  cout << number_of_passed_tests << " test(s) passed, " <<
  number_of_failed_tests;
  cout << " test(s) failed." << endl << endl;

  // End of tests
  //
  return EXIT_SUCCESS;
}
