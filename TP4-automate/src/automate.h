#ifndef __AUTOMATE_H__
#define __AUTOMATE_H__

#include <iostream>
#include <map>
#include <vector>
#include <cstring>

#include "matrix.h"
#include "parser.h"


class edge_t {
public:
	edge_t(std::string source, std::string target, std::string trans_name);
	std::string get_source();
	std::string get_target();
	std::string get_trans_name();
private:
	std::string    _source;
	std::string    _target;
	std::string    _trans_name;
};

class substitute_t {
public:
	substitute_t();
	void add_constraint(std::string a, std::string b);
	std::string substitute(std::string tosubst);
	friend std::ostream& operator<<(std::ostream& s, substitute_t& subst);
	void sort_data();
private:
	std::vector<std::vector<std::string> > _data;
};

class automate_t {
public:
	automate_t(char *filename, std::string prefix="");
	std::string get_initial_state();
	int nadjacent(std::string state);
	edge_t get_edge(std::string state, int idx_adjacent);
	int nprimitives(std::string state);
	fmatrix get_primitive(std::string state, int idx);
	friend std::ostream& operator<<(std::ostream& s, automate_t& a);

	smatrix get_strans(std::string trans_name);
	fmatrix get_ftrans(std::string trans_name);
	bool is_determined();
	void fill_free_variables();

private:
	std::string _initial_state;
	std::map<std::string, std::vector<edge_t>  > _adjacency_map;
	std::map<std::string, std::vector<fmatrix> > _primitives;
	std::map<std::string, smatrix > _strans;
	substitute_t _substitute;

	smatrix _multiply(cstr string);
	void _parse_automate_buffer(cstr buffer, std::string prefix);
	void _parse_constraints_buffer(cstr line);
	void _include_automate(cstr line, std::string prefix);
};

#endif //__AUTOMATE_H__
