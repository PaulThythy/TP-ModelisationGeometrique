#include <iostream>
#include <sstream>
#include <cstring>
#include <stdexcept>

#include "automate.h"


edge_t::edge_t(std::string source, std::string target, std::string trans_name) : _source(source), _target(target), _trans_name(trans_name) {
}

std::string edge_t::get_source() {
	return _source;
}

std::string edge_t::get_target() {
	return _target;
}

std::string edge_t::get_trans_name() {
	return _trans_name;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

substitute_t::substitute_t() {
}

std::string substitute_t::substitute(std::string tosubst) {
	for (unsigned int i=0; i<_data.size(); i++) {
		for (unsigned int j=0; j<_data[i].size(); j++) {
			if (_data[i][j]==tosubst) {
				return _data[i][0];
			}
		}
	}
	return tosubst;
}

void substitute_t::sort_data() {
	// just for fun let us sort the constraints
	for (unsigned int i=0; i<_data.size(); i++) {
		for (unsigned int j=0; j<_data[i].size(); j++) {
			for (unsigned int k=j+1; k<_data[i].size(); k++) {
				if (_data[i][k] < _data[i][j]) {
					std::string z = _data[i][j];
					_data[i][j] = _data[i][k];
					_data[i][k] = z;
				}
			}
		}

		if (_data[i].size()>=2) {
			float a;
			float b;
			bool afloat = parseT(_data[i][0].c_str(), a);
			bool bfloat = parseT(_data[i][1].c_str(), b);
			if (afloat && bfloat) {
				if (a!=b) std::cerr << *this << std::endl;
				assert(a==b);
			}
		}
	}
}

void substitute_t::add_constraint(std::string a, std::string b) {
	if (a==b) return; // nothing to add
	float foo;

	int aidx = -1;
	int bidx = -1;
	for (unsigned int i=0; (-1==aidx||-1==bidx) && i<_data.size(); i++) {
		for (unsigned int j=0; (-1==aidx||-1==bidx) && j<_data[i].size(); j++) {
			if (_data[i][j]==a) aidx = i;
			if (_data[i][j]==b) bidx = i;
		}
	}

	if (aidx==bidx && -1!=aidx) return; // processed
	if (-1==aidx && -1==bidx) {         // both do not exist
		std::vector<std::string> z;
		z.push_back(a);
		if (parseT(b.c_str(), foo))
			z.insert(z.begin(), b);
		else
			z.push_back(b);
		_data.push_back(z);
	} else if (-1==aidx) {              // have found b, but not a
		bool flag = true;
		for (unsigned int k=0; k<_data[bidx].size(); k++)
			if (_data[bidx][k]==a) flag = false;
		if (flag)
//			if (parseT(a.c_str(), foo))
//				_data[bidx].insert(_data[bidx].begin(), a);
//			else
				_data[bidx].push_back(a);
	} else if (-1==bidx) {              // have found a, but not b
		bool flag = true;
		for (unsigned int k=0; k<_data[aidx].size(); k++)
			if (_data[aidx][k]==b) flag = false;
		if (flag)
//			if (parseT(b.c_str(), foo))
//				_data[aidx].insert(_data[aidx].begin(), b);
//			else
				_data[aidx].push_back(b);
	} else {                            // both are found in different columns, let us merge the columns
		for (unsigned int j=0; j<_data[bidx].size(); j++) {
			bool flag = true;
			for (unsigned int k=0; k<_data[aidx].size(); k++)
				if (_data[aidx][k]==_data[bidx][j]) flag = false;
			if (flag)
//				if (parseT(_data[bidx][j].c_str(), foo))
//					_data[aidx].insert(_data[aidx].begin(), _data[bidx][j]);
//				else
					_data[aidx].push_back(_data[bidx][j]);
		}
		_data.erase(_data.begin()+bidx);
	}

}

std::ostream& operator<<(std::ostream& s, substitute_t& subst) {
	for (unsigned int i=0; i<subst._data.size(); i++) {
		for (unsigned int j=0; j<subst._data[i].size(); j++) {
			s << subst._data[i][j];
			if (j < subst._data[i].size()-1) s << " = ";
		}
		s << std::endl;
	}
	return s;
}

smatrix automate_t::_multiply(cstr string) {
	std::string s;
	smatrix cm;
	bool first = true;
	while ((string = parseT(string, s))) {
		std::map<std::string, smatrix >::iterator iter = _strans.find(s);
		smatrix temp;
		if (iter==_strans.end()) {
			temp = smatrix(1,1);
			temp[0][0] = s;
		} else {
			temp = iter->second;
		}
		cm = first ? temp : cm*temp;
		first = false;
	}
	return cm;
}

void automate_t::_include_automate(cstr line, std::string prefix) {
//	if ("undef"==_initial_state) _initial_state = trans_name;
	const char *b = strchr(line, '<');
	const char *e = strchr(line, '>');
	assert (b && e && e>b);
	std::string state_name = prefix + std::string(line, b-line);
	std::string file_name(b+1, e-b-1);
//	std::cout << state_name << std::endl << file_name << std::endl;
	automate_t include = automate_t((char *)(file_name.c_str()));
	std::string inc_in_st = include.get_initial_state();
//	for (std::map<std::string, fmatrix >::iterator iter = include._primitives.begin(); iter!=include._primitives.end(); ++iter) {
	for (std::map<std::string, std::vector<fmatrix> >::iterator iter = include._primitives.begin(); iter!=include._primitives.end(); ++iter) {
		std::string t_state_name = (inc_in_st==iter->first) ? state_name : state_name+iter->first ;
		_primitives[t_state_name] = iter->second;
	}
	for (std::map<std::string, smatrix >::iterator iter = include._strans.begin(); iter!=include._strans.end(); ++iter) {
		std::string t_state_name = (inc_in_st==iter->first) ? state_name : state_name+iter->first ;
		smatrix c = iter->second;
		for (int i=0; i<c.nrows(); i++) {
			for (int j=0; j<c.ncols(); j++) {
				std::string z = c[i][j].value();
				if ('{'==z[0]) {
					//.TODO probably a bug with initial state
					z.insert(1, state_name);
					c[i][j] = string_element_t(z);
				}
			}
		}
		//.TODO: {prefix
		_strans[t_state_name] = c;
	}
	for (std::map<std::string, std::vector<edge_t> >::iterator iter = include._adjacency_map.begin(); iter!=include._adjacency_map.end(); ++iter){
		int size = iter->second.size();
		for (int i=0; i<size; i++) {
			edge_t edge = iter->second[i];
			std::string source = edge.get_source();
			source = (source==inc_in_st) ? state_name : state_name+source;
			std::string target = edge.get_target();
			target = (target==inc_in_st) ? state_name : state_name+target;
			std::string trans_name = state_name+edge.get_trans_name();
			_adjacency_map[source].push_back(edge_t(source, target, trans_name));
//			smatrix z = a.get_strans(edge.get_trans_name());
//			s << edge.get_trans_name() << ":" << iter->first << "-" << edge.get_target() << std::endl << z << std::endl;
		}
	}
//	std::cout << include << std::endl;
}

void automate_t::_parse_constraints_buffer(cstr line) {
	const char *equal_mark = strchr(line, '=');
	if (!equal_mark) return;
	char *left = new char[equal_mark-line+1];
	strncpy(left, line, equal_mark-line);
	left[equal_mark-line] = '\0';
	smatrix cm  = _multiply(left);
	smatrix cm2 = _multiply(equal_mark+1);
	assert(cm.nrows()==cm2.nrows() && cm.ncols()==cm2.ncols());

	std::cout << "===========================\n";
	std::cout << line << std::endl << std::endl << cm << std::endl << cm2;

	for (int i=0; i<cm.nrows(); i++) {
		for (int j=0; j<cm.ncols(); j++) {
			assert(CANT_PROCESS!=cm[i][j].value() && CANT_PROCESS!=cm2[i][j].value());
			_substitute.add_constraint(cm[i][j].value(), cm2[i][j].value());
		}
	}
	std::cout << "===========================\n";
	delete [] left;
}

void automate_t::_parse_automate_buffer(cstr buffer, std::string prefix) {
	int max_line_size = 1024;
	char *line = new char[max_line_size];
	while (buffer && strlen(buffer)>0) {
		memset(line, '\0', max_line_size);
		buffer = skip_spaces(buffer, true);
		buffer = get_line(buffer, line, max_line_size, false);

		if (strchr(line, '=')) {
			_parse_constraints_buffer(line);
			continue;
		}

		if (strchr(line, '<') && strchr(line, '>')) {
			_include_automate(line, prefix);
			continue;
		}
		char *comma = strchr(line, ':');
		if (!comma) continue; // not a valid transformation/primitives description
		std::string trans_name = prefix + std::string(line, comma-line);
		if ("undef"==_initial_state) _initial_state = trans_name;
		char *dash  = strchr(line, '-');
		if (!dash) { // just some primitives to draw
			fmatrix foo;
//			buffer = parseT(buffer, _primitives[trans_name], true);
			buffer = parseT(buffer, foo, true);
			_primitives[trans_name].push_back(foo);
			continue;
		}

		// ok. we have a transformation
		smatrix curmatrix;
		buffer = parseT(buffer, curmatrix, true);
		if (!buffer) continue; // if there is a matrix

		std::string source(comma+1, dash-comma-1);
		std::string target(dash+1);
		if (1==curmatrix.nrows() && 1==curmatrix.ncols() && curmatrix[0][0]=="?") { // auto-fill the smatrix
//			curmatrix = smatrix(_primitives[source].nrows(), _primitives[target].nrows());
			//.TODO check for existence [0]
			curmatrix = smatrix(_primitives[source][0].nrows(), _primitives[target][0].nrows());
			for (int i=0; i<curmatrix.nrows(); i++) {
				for (int j=0; j<curmatrix.ncols(); j++) {
					curmatrix[i][j] = string_element_t("{"+trans_name + "_" +(char)(j+'a')+(char)(i+'a')+"}");
				}
			}
		}
		_adjacency_map[source].push_back(edge_t(source, target, trans_name));
		_strans[trans_name] = curmatrix;
	}
	delete [] line;
}

automate_t::automate_t(char *filename, std::string prefix) : _initial_state("undef") {
	std::ifstream in;
	in.open (filename, std::ifstream::in);
	if (in.fail()) return;
	in.seekg (0, std::ios::end);
	int length = in.tellg();
	in.seekg (0, std::ios::beg);
	char *buffer = new char [length];
	in.read (buffer,length);
	in.close();
	_parse_automate_buffer(buffer, prefix);
	delete [] buffer;
//	std::cout << _substitute << std::endl;
	_substitute.sort_data();
	for (std::map<std::string, smatrix >::iterator iter = _strans.begin(); iter!=_strans.end(); ++iter) {
		for (int i=0; i<iter->second.nrows(); i++)
			for (int j=0; j<iter->second.ncols(); j++)
				iter->second[i][j] = _substitute.substitute(iter->second[i][j].value());
	}

//	std::cout << *this << std::endl;
}


std::string automate_t::get_initial_state() {
	return _initial_state;
}

int automate_t::nadjacent(std::string state) {
	return _adjacency_map[state].size();
}

int automate_t::nprimitives(std::string state) {
	return _primitives[state].size();
}

edge_t automate_t::get_edge(std::string state, int idx_adjacent) {
	return _adjacency_map[state][idx_adjacent];
}

fmatrix automate_t::get_primitive(std::string state, int idx) {
	return _primitives[state][idx];
}

smatrix automate_t::get_strans(std::string trans_name) {
	return _strans[trans_name];
}

fmatrix automate_t::get_ftrans(std::string trans_name) {
	return s2fmatrix(_strans[trans_name]);
//	return _ftrans[trans_name];
}

std::ostream& operator<<(std::ostream& s, automate_t& a) {
//	for (std::map<std::string, fmatrix >::iterator iter = a._primitives.begin(); iter!=a._primitives.end(); ++iter) {
//		s << iter->first << ":" << std::endl << iter->second << std::endl;
	for (std::map<std::string, std::vector<fmatrix> >::iterator iter = a._primitives.begin(); iter!=a._primitives.end(); ++iter) {
		int size = iter->second.size();
		for (int i=0; i<size; i++) {
			s << iter->first << ":" << std::endl << iter->second[i] << std::endl;
		}
	}
	for (std::map<std::string, std::vector<edge_t> >::iterator iter = a._adjacency_map.begin(); iter!=a._adjacency_map.end(); ++iter){
		int size = iter->second.size();
		for (int i=0; i<size; i++) {
			edge_t edge = iter->second[i];
			smatrix z = a.get_strans(edge.get_trans_name());
			s << edge.get_trans_name() << ":" << iter->first << "-" << edge.get_target() << std::endl << z << std::endl;
		}
	}
	s << a._substitute << std::endl;
	return s;
}

bool automate_t::is_determined() {
	float foo;
	for (std::map<std::string, std::vector<edge_t> >::iterator iter = _adjacency_map.begin(); iter!=_adjacency_map.end(); ++iter){
		int size = iter->second.size();
		for (int i=0; i<size; i++) {
			edge_t edge = iter->second[i];
			smatrix z = get_strans(edge.get_trans_name());
			for (int i=0; i<z.nrows(); i++) {
				for (int j=0; j<z.ncols(); j++) {
					if (!parseT(z[i][j].value().c_str(), foo))
						return false;
				}
			}
		}
	}
	return true;
}
class BadConversion : public std::runtime_error {
public:
	BadConversion(const std::string& s) : std::runtime_error(s) { }
};

inline std::string stringify(float x) {
	std::ostringstream o;
	if (!(o << x))
		throw BadConversion("stringify(double)");
	return o.str();
}
void automate_t::fill_free_variables() {
	float foo;
	for (std::map<std::string, std::vector<edge_t> >::iterator iter = _adjacency_map.begin(); iter!=_adjacency_map.end(); ++iter){
		int size = iter->second.size();
		for (int i=0; i<size; i++) {
			edge_t edge = iter->second[i];
			smatrix z = get_strans(edge.get_trans_name());
			for (int j=0; j<z.ncols(); j++) {
				int cnt = 0;
				for (int i=0; i<z.nrows(); i++) {
					if (!parseT(z[i][j].value().c_str(), foo))
						cnt++;
				}
				for (int i=0; i<z.nrows(); i++) {
					if (!parseT(z[i][j].value().c_str(), foo)) {
						_strans[edge.get_trans_name()][i][j] = string_element_t(stringify(1.0/cnt));
					}
				}
			}
		}
	}
}

