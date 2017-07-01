/*
 * CubicSpline.h
 *
 *  Created on: Nov, 2015
 *      Author: Tomislav
 */

#ifndef CUBICSPLINE_H_
#define CUBICSPLINE_H_

#include <vector>
#include <utility>
#include <array>
#include <memory>
#include <stdexcept>
#include <ostream>

class CubicSpline {
public:
	class NodesPrinter {
	public:
		NodesPrinter()=delete;
		~NodesPrinter();
		friend std::ostream & operator<<(std::ostream &, NodesPrinter const &);
		friend CubicSpline;
	private:
		CubicSpline const * v_spline;
		bool mutable cached;
		std::string mutable cached_nodes;
		NodesPrinter(CubicSpline const * const);
	};

	friend NodesPrinter;

	friend std::ostream & operator<<(std::ostream &, NodesPrinter const &);

	static const int DEGREE = 3;

	CubicSpline(CubicSpline const &) = delete;
	CubicSpline(CubicSpline &&) = delete;
	CubicSpline & operator=(CubicSpline const &) = delete;
	CubicSpline & operator=(CubicSpline &&) = delete;

	const NodesPrinter nodes_printer;

	using node = std::pair<double, double>;
	using function = std::vector<node>;
	using function_ptr = std::unique_ptr<function>;
	using polynome = std::array<double, CubicSpline::DEGREE + 1>;
	using spline = std::vector<polynome>;
	using interval = std::pair<double, double>;

	CubicSpline(function_ptr, double, double);
	double operator()(double) const;
	inline bool operator[](double) const;

	inline interval get_interval_defined() const;
	inline double min_defined() const;
	inline double max_defined() const;

	inline bool is_in_domain(double) const;

	~CubicSpline();

private:
	function_ptr v_nodes;
	spline v_coefficients;

	using index_type = function::size_type;

	void build_spline(double, double);
	void check_is_function() const;
	void check_is_interpolable() const;
	index_type find_interval(double) const;
};

CubicSpline::interval CubicSpline::get_interval_defined() const {
	return interval(min_defined(), max_defined());
}

bool CubicSpline::is_in_domain(double argument) const {
	return argument >= min_defined() && argument <= max_defined();
}

bool CubicSpline::operator[](double argument) const {
	return is_in_domain(argument);
}

double CubicSpline::min_defined() const {
	return v_nodes->front().first;
}

double CubicSpline::max_defined() const {
	return v_nodes->back().first;
}

std::ostream & operator<<(std::ostream &, CubicSpline::NodesPrinter const &);

#endif /* CUBICSPLINE_H_ */
