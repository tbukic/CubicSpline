/*
 * CubicSpline.cpp
 *
 *  Created on: Nov, 2015
 *      Author: Tomislav
 */

#include "CubicSpline.h"

#include <algorithm>
#include <string>
#include <cstdio>
#include <sstream>
#include <iostream>

CubicSpline::CubicSpline(
		function_ptr nodes,
		double first_node_2nd_d,
		double last_node_2nd_d) :
		nodes_printer(this), v_nodes(std::move(nodes)), v_coefficients() {
	if (!v_nodes) {
		throw std::logic_error("nullptr in constructor.");
	}

	std::sort(v_nodes->begin(), v_nodes->end());
	check_is_interpolable();
	check_is_function();

	build_spline(first_node_2nd_d, last_node_2nd_d);
}

void CubicSpline::check_is_interpolable() const {
	const size_t NO_NODES = v_nodes->size();
	if (0 == NO_NODES) {
		throw std::runtime_error("Empty function, nothing to interpolate.");
	}
	if (2 == NO_NODES) {
		throw std::runtime_error(
				"Can not interpolate cubic spline though only 2 points."
		);
	}
}

void CubicSpline::check_is_function() const {

	const int last_index = v_nodes->size() - 1;
	for (int index = 0; index < last_index; ++index) {
		if ((*v_nodes)[index].first == (*v_nodes)[index + 1].first) {
			throw::std::runtime_error("Node " + std::to_string((*v_nodes)[index].first)
				+ " is given two times.");
		}
	}
}

void CubicSpline::build_spline(double first_node_2nd_d, double last_node_2nd_d) {
	const index_type size = v_nodes->size() - 1;
	if (size == 0) {
		return;
	}

	std::vector<double> h_i(size);
	for (index_type index = 0; index < size; ++index) {
		h_i[index] = (*v_nodes)[index + 1].first - (*v_nodes)[index].first;
	}

	const index_type MATRIX_DIM = size - 1;
	std::vector<double> d_i(MATRIX_DIM);
	double last, current = ((*v_nodes)[1].second - (*v_nodes)[0].second) / h_i[0];
	for (index_type index = 1; index < size; ++index) {
		last = current;
		current = ((*v_nodes)[index + 1].second - (*v_nodes)[index].second) / h_i[index];
		d_i[index - 1] = 6.0 * (current - last) / (h_i[index] + h_i[index - 1]);
	}

	std::vector<double> m(size + 1);
	m[0] = first_node_2nd_d;
	m[size] = last_node_2nd_d;

	const index_type ALPHA = 0, DIAGONAL = 1, BETA = 2;

	using matrix = std::vector<std::vector<double> >;
	matrix trid_matrix(MATRIX_DIM);
	for (index_type index = 0; index < MATRIX_DIM; ++index) {
		trid_matrix[index].reserve(3);
		trid_matrix[index][ALPHA] = h_i[index]/(h_i[index] + h_i[index + 1]);
		trid_matrix[index][DIAGONAL] = 2.0;
		trid_matrix[index][BETA] = 1.0 - trid_matrix[index][ALPHA];
	}

	d_i[0] = d_i[0] - trid_matrix[0][ALPHA]*m[0];
	d_i[MATRIX_DIM - 1] = d_i[MATRIX_DIM - 1] - trid_matrix[MATRIX_DIM - 1][BETA]*m[size];

	double coefficient;
	for (index_type index = 1; index < MATRIX_DIM; ++index) {
		coefficient = trid_matrix[index][ALPHA]/trid_matrix[index - 1][DIAGONAL];
		// trid_matrix[index][ALPHA] = 0.0; // value won't be used, so no need for spending cycles.
		trid_matrix[index][DIAGONAL] -= (coefficient * trid_matrix[index - 1][BETA]);
		d_i[index] -= (coefficient * d_i[index - 1]);
	}

	m[MATRIX_DIM] = d_i[MATRIX_DIM - 1]/trid_matrix[MATRIX_DIM - 1][DIAGONAL];
	for (
			index_type index = MATRIX_DIM - 1, index_1 = index - 1;
			index > 0;
			--index, --index_1
	) {
		coefficient = trid_matrix[index_1][BETA]/trid_matrix[index][DIAGONAL];
		// trid_matrix[index_1][BETA] = 0.0; // irrelevant and never used, like in last loop
		// trid_matrix[index_1][DIAGONAL] -= (coefficient * trid_matrix[index][ALPHA]); // doesn't change
		d_i[index_1] -= (coefficient * d_i[index]);
		m[index] = d_i[index_1]/trid_matrix[index_1][DIAGONAL];
	}

	const index_type no_polynoms = v_nodes->size() - 1;
	v_coefficients.resize(no_polynoms);
	double x_i, s_i, s_i_min1, first_part, second_part, third_part;
	for (index_type index = 0; index < no_polynoms; ++index) {
		polynome & pol = v_coefficients[index];

		x_i = (*v_nodes)[index].first;
		s_i = (*v_nodes)[index + 1].second;
		s_i_min1 = (*v_nodes)[index].second;

		first_part = (s_i - s_i_min1)/h_i[index] -
				(m[index + 1] + 2.0 * m[index]) * h_i[index] / 6.0;
		second_part = m[index] / 2.0;
		third_part = (m[index + 1] - m[index]) / (6.0*h_i[index]);


		pol[0] = third_part;
		pol[1] = -third_part * 3.0 * x_i + second_part;
		pol[2] = x_i * (3.0 * x_i * third_part - 2.0 * second_part)
				+ first_part;
		pol[3] = s_i_min1 + x_i * (x_i * (second_part - x_i * third_part) - first_part);
	}
}

CubicSpline::index_type CubicSpline::find_interval(double argument) const {
	if (!(*this)[argument]) { // eliminates values in R\[first_node, last_node]
		throw new std::runtime_error(
				"Spline is not defined at " + std::to_string(argument) + "."
		);
	}
	index_type first = 0;
	index_type last = v_nodes->size() - 1;
	index_type middle;

	if ((*v_nodes)[last].first == argument) { // eliminates last
		return last;
	}

	while (first < last - 1) { // searches for maximal index, so node[index]<=argument
		middle = first + (last - first)/2;
		if ((*v_nodes)[middle].first <= argument) {
			first = middle;
		} else {
			last = middle;
		}
	}
	return first;  // index, so that argument is in interval [node[index], node[index+1]>
}

double CubicSpline::operator()(double argument) const {
	index_type interval_index = find_interval(argument);
	if ((*v_nodes)[interval_index].first == argument) {
		return (*v_nodes)[interval_index].second;
	}
	double approximation = 0.0;
	for (double coefficient : v_coefficients[interval_index]) {
		approximation = approximation * argument + coefficient;
	}
	return approximation;
}

CubicSpline::NodesPrinter::NodesPrinter(CubicSpline const * const spline) : v_spline(spline), cached(false) {
}

CubicSpline::NodesPrinter::~NodesPrinter() {
	v_spline = nullptr;
}


std::string format_at_length(std::string original, size_t length) {
	size_t original_length = original.length();
	if (original_length == length) {
		return original;
	}
	std::string formated = "";
	formated.append((length - original.length()) / 2, ' ');
	formated.append(original);
	formated.append(length - formated.length(), ' ');
	return formated;
}

std::ostream & operator<<(std::ostream & out, CubicSpline::NodesPrinter const & n_printer) {
	if (!n_printer.cached) {
		std::stringstream node_builder;

		const size_t length = n_printer.v_spline->v_nodes->size();
		std::vector<std::string> x_axis, y_axis;

		x_axis.reserve(length);
		y_axis.reserve(length);
		std::vector<size_t> lengths;

		for (size_t index = 0; index < length; ++index) {
			std::ostringstream str_x;
			str_x << (*(n_printer.v_spline->v_nodes))[index].first;
			x_axis.push_back(str_x.str());
			std::ostringstream str_y;
			str_y << (*(n_printer.v_spline->v_nodes))[index].second;
			y_axis.push_back(str_y.str());
			lengths.push_back(
					std::max(
							x_axis[index].length(),
							y_axis[index].length()
					)
			);
		}

		node_builder << "x";
		for (size_t index = 0; index < length; ++index) {
			node_builder << " | " << format_at_length(x_axis[index], lengths[index]);
		}
		node_builder << std::endl << "-";
		for (size_t index = 0; index < length; ++index) {
			std::string separator;
			separator.append(lengths[index], '-');
			node_builder << "-+-" << separator;
		}
		node_builder << "-" << std::endl << "y";
		for (size_t index = 0; index < length; ++index) {
			node_builder << " | " << format_at_length(y_axis[index], lengths[index]);
		}

		n_printer.cached_nodes = node_builder.str();
		n_printer.cached = true;
	}
	out << n_printer.cached_nodes << std::endl;

	return out;
}

CubicSpline::~CubicSpline() {
}
