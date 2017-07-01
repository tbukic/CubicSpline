/*
 * main.cpp
 *
 *  Created on: Nov, 2015
 *      Author: Tomislav
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <memory>
#include <utility>
#include <string>
#include <stdexcept>
#include <sstream>

#include "CubicSpline/CubicSpline.h"


using node = CubicSpline::node;
using function = CubicSpline::function;
using function_ptr = CubicSpline::function_ptr;
using second_derivatives = std::pair<double, double>;
using function_data = std::pair<
		function_ptr,
		std::unique_ptr<second_derivatives>
	>;
using spline_ptr = std::unique_ptr<CubicSpline>;


std::istream & operator>>(std::istream & input, function_data & data) {
	data.first->clear();

	std::string readed;
	while (std::getline(input, readed)) {
		if (readed.empty()) {		// reads until empty line
			break;
		}

		std::istringstream s_stream(readed);
		double point, value;
		if (!(s_stream >> point) || !(s_stream >> value) || !s_stream.eof()) {
			throw std::runtime_error("Bad input!");
		}
		data.first->push_back(node(point, value));
	}

	std::getline(input, readed);
	std::istringstream s_stream(readed);
	if (
			!(s_stream >> data.second->first) ||
			!(s_stream >> data.second->second) ||
			!s_stream.eof()
	) {
		throw std::runtime_error("Bad input!");
	}

	return input;
}

void run_interpolation(spline_ptr spline) {

	const std::string EXIT_INPUT = "exit", PROMPT = ">";
	std::cout << "Function is defined on interval [" <<
			spline->min_defined() << ", " << spline->max_defined()
			<< "]." << std::endl << "Nodes are:"
			<< std::endl << std::endl << spline->nodes_printer
			<< std::endl  << "Interpolation starts. To quit program, "
			<< "type: '" << EXIT_INPUT << "'. " << std::endl;
	std::string input;

	bool show_message = true;
	while(true) {
		if (show_message) {
			std::cout << "Enter value you want to interpolate:" << std::endl;
		} else {
			show_message = true;
		}
		std::cout << PROMPT;
		std::getline(std::cin, input);
		if (EXIT_INPUT == input) {
			return;
		}
		if ("nodes" == input) {
			std::cout << spline->nodes_printer;
			continue;
		}
		std::stringstream s_stream(input);
		double x;
		if (!(s_stream >> x) || !s_stream.eof()) {
			std::cout << "Input is not a number." << std::endl;
			continue;
		}
		if (!(*spline)[x]) {
			std::cout << "Value " << x << " is not in domain of the spline."
					<< std::endl << "Spline is defined on interval [" <<
					spline->min_defined() << ", " << spline->max_defined()
					<< "]." << std::endl;
			continue;
		}
		std::cout << "spline(" << x << ")=" << (*spline)(x) << std::endl;
		show_message = false;
	}

}

int main(int argc, char * argv[]) {

	function_data data = std::make_pair(
			function_ptr(new function()),
			std::unique_ptr<second_derivatives>(new second_derivatives())
	);

	if (argc == 1) {			// reading from console:
		try {
			std::cout << "Enter data for nodes and second derivatives."
					 << std::endl << "In every row, enter one node and it's"
					 << " functional value, separated by blanks." << std::endl
					 << "After entering data for all nodes, enter empty line,"
					 << " and after that," << std::endl
					 << "enter values of second derivatives "
					 << "in first and last node, separated by blanks."
					 << std::endl;

			std::cin >> data;
		} catch (std::exception & e) {
			std::cout << "Error in input." << std::endl;
			return EXIT_FAILURE;
		}

	} else if (argc == 2) {		// reading from file (speeds up testing!):
		std::ifstream input(argv[1]);
		if (!input.good()) {
			input.close();
			std::cout << "File does not exist." << std::endl;
			return EXIT_FAILURE;
		}
		try {
			input >> data;
		} catch (std::exception & e) {
			input.close();
			std::cout << "File is badly formated." << std::endl;
			return EXIT_FAILURE;
		}
		input.close();

	} else {					// bad arguments
		std::cout <<
				"Program doesn't work with more than 1 argument."
				<< std::endl;
		return EXIT_FAILURE;
	}

	try {
		run_interpolation(
				std::move(spline_ptr(new CubicSpline(
						std::move(data.first),
						data.second->first,
						data.second->second
				)))
		);
		data.first = nullptr;

	} catch (std::runtime_error & e){
		std::cout << e.what() << std::endl;
		return EXIT_FAILURE;
	} catch (...) {
		std::cout << "Error in program." << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
