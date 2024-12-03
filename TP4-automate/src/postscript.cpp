#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>

#include "automate.h"
#include "matrix.h"
#include "parser.h"
//#include "image.h"

#define SIZE_X 640
#define SIZE_Y 640
#define STEPS  3
int PRIMITIVES_ARE_CLOSED = 0;
int NSTEPS = 0;
int VRML   = 0;

struct current_state_t {
	std::string address;
	std::string state;
	fmatrix trans_matrix;
	current_state_t(std::string arg_address, std::string arg_state, fmatrix arg_trans_matrix) : address(arg_address), state(arg_state), trans_matrix(arg_trans_matrix) {}
	current_state_t() : address(""), state("z"), trans_matrix() {}
};

void iterate(automate_t &automate, fmatrix f, std::string state, int depth, std::ofstream &out) {
	if (depth+1<=NSTEPS) {
		int nedges = automate.nadjacent(state);
		for (int i=0; i<nedges; i++) {
			edge_t edge = automate.get_edge(state, i);
			fmatrix foo = automate.get_ftrans(edge.get_trans_name());
			iterate(automate, f*foo, edge.get_target(), depth+1, out);
		}
	} else {
		int nprimitives = automate.nprimitives(state);
		for (int i=0; i<nprimitives; i++) {
			fmatrix primitives = automate.get_primitive(state, i);
			primitives = f*primitives;

			if (VRML) {
				if (primitives.ncols()<2) continue;
				out << "Separator { Coordinate3 {\npoint [";
				for (int j=0; j<primitives.ncols(); j++)
					out << primitives[0][j] << " " << primitives[1][j] << " " << (primitives.nrows()>2 ? primitives[2][j] : 1) << ", ";
				if (PRIMITIVES_ARE_CLOSED)
					out << primitives[0][0] << " " << primitives[1][0] << " " << (primitives.nrows()>2 ? primitives[2][0] : 1) << ", ";
				out << "]}\n";

				out << "ShapeHints { vertexOrdering CLOCKWISE shapeType SOLID}";
				out << "Material { diffuseColor  1 0.2 0.2 }";
				out << "FaceSet{numVertices["<< (PRIMITIVES_ARE_CLOSED ? primitives.ncols()+1 : primitives.ncols()) << "]}\n";
				out << "Material { diffuseColor  1 0.8 0 }";
				out << "ShapeHints { vertexOrdering COUNTERCLOCKWISE shapeType SOLID}";
				out << "FaceSet{numVertices["<< (PRIMITIVES_ARE_CLOSED ? primitives.ncols()+1 : primitives.ncols()) << "]}\n";
				out << "}\n";
			} else {
				for (int j=0; j<primitives.ncols(); j++)
					out << primitives[0][j]*SIZE_X << " " << primitives[1][j]*SIZE_Y << " " << ((0==j) ? "m" : "l") << " ";
				if (1 != primitives.ncols()) out << (PRIMITIVES_ARE_CLOSED ? " CL" : "" ) << " S\n";
			}
		}
	}
}

int main(int argc, char **argv) {
	if (2>argc) { std::cout << "Usage: " << argv[0] << " automate.txt" << std::endl; exit(1); }
	if (argc==2 || !parseT(argv[2], NSTEPS)) {
		NSTEPS = STEPS;
	}

	automate_t automate = automate_t(argv[1]);
	std::cout << automate << std::endl;

	if (!automate.is_determined()) {
		automate.fill_free_variables();
		std::cout << "the automate contains non-defined variables" << std::endl;
//		exit(1);
	}

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///   non-recursive + postscript  //////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::stringstream buf;
	buf << time(NULL) << ".ps";

	std::ofstream out;
	out.open(buf.str().c_str());
	out << "%!PS-Adobe-3.0\n%%BoundingBox: 100 0 " << (SIZE_X-100) << " " << SIZE_Y << "\n%Pages: (atend)\n/circ{8 0 360 arc closepath fill}def /l{lineto}def /m{moveto}def /rm{rmoveto}def /rl{rlineto}def /S{stroke}def /CL{closepath}def /RGB{setrgbcolor}def /BLACK{0.1 0.1 0 RGB}def /GREEN{0 0 0 RGB}def /BLUE{0.174 1 1 RGB}def /RED{1 0 0 RGB}def\n\n%Page: 1\n0.1 setlinewidth\n";

	std::vector<current_state_t> l1,l2;
	l1.push_back(current_state_t("", automate.get_initial_state(), fmatrix()));
	for (int step=0; step<NSTEPS; step++) {
		std::cerr << step << "\n";
		for (int istate=0; istate<(int)l1.size(); istate++) {
			int nedges = automate.nadjacent(l1[istate].state);
			for (int i=0; i<nedges; i++) {
				edge_t edge = automate.get_edge(l1[istate].state, i);
				fmatrix foo = automate.get_ftrans(edge.get_trans_name());
				l2.push_back(current_state_t(l1[istate].address + "-" + edge.get_trans_name(), edge.get_target(), (0==step) ? foo : l1[istate].trans_matrix*foo));
			}
		}
		l1 = l2;
		l2.clear();
	}

	out << "showpage\n%Trailer\n%Pages: 1\n";
	out.close();
	l1.clear();

	std::stringstream buf2;
	buf2 << "gv " << buf.str();
	system(buf2.str().c_str());

*/
/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///  vrml  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::stringstream buf;
	buf << time(NULL) << ".iv";

	std::ofstream out;
	out.open(buf.str().c_str());
	out << "#Inventor V2.1 ascii\nSeparator { DrawStyle {pointSize 4 lineWidth 1 }\n";

	out <<  "Separator { Coordinate3 {\npoint [";
	fmatrix primitives = automate.get_ftrans("P");
	std::cout << primitives << std::endl;
	for (int j=0; j<primitives.ncols(); j++) {
		out << primitives[0][j] << " " << primitives[1][j] << " " <<primitives[2][j] << ", ";
	}
//	out << primitives[0][0] << " " << primitives[1][0] << " " <<primitives[2][0] << ", ";
	out << "]} PointSet{}}\n";
	for (int i=0; i<(int)l1.size(); i++) {
		out << "Separator { Coordinate3 {\npoint [";
		std::cout << l1[i].state << ": "<< l1[i].address << std::endl;
		fmatrix primitives = automate.get_primitives(l1[i].state);
		primitives = l1[i].trans_matrix*primitives;
		for (int j=0; j<primitives.ncols(); j++) {
				out << primitives[0][j] << " " << primitives[1][j] << " " <<primitives[2][j] << ", ";
		}
		out << primitives[0][0] << " " << primitives[1][0] << " " <<primitives[2][0] << ", ";
		out << "]} LineSet{numVertices["<< (primitives.ncols()+1) << "]}}\n";
	}
	out << "}\n";
	out.close();
	l1.clear();

	std::stringstream buf2;
	buf2 << "ivview " << buf.str();
	system(buf2.str().c_str());
	*/


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///  recursive + postscript  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::stringstream buf;
	buf << time(NULL) << (VRML ? ".iv" : ".ps");

	std::ofstream out;
	out.open(buf.str().c_str());

	if (VRML) {
		out << "#Inventor V2.1 ascii\nSeparator { DrawStyle {pointSize 4 lineWidth 2 }\n";

		out << "Material { ambientColor 1 1 1 diffuseColor  1.7 1.7 1.7 }";

//		out << "PointLight { location 0.5 0.5 3 color 1 1 .1 }";
//		out << "PerspectiveCamera { \n";

//out <<    "position 0.01553308 -2.4726472 5.1244702 orientation 0.99997836 0.0065651615 -0.0011135775  0.49828756";
//out << "    position -0.0071311165 -3.5054817 3.530942 orientation 0.99999475 -0.0027790614 0.0016285594  0.85787272";




//out << "    position 2.0601101 1.5432454 -0.58293366 orientation 0.85281414 -0.45001197 0.26494759  4.0525208 nearDistance 1.3148628 farDistance 3.0204513 focalDistance 2.1668057";

//out <<"    position 2.2348623 2.2255433 0.50505668 orientation 0.28894052 0.67937124 0.67451316  2.551384 nearDistance 1.5669345 farDistance 3.1559346 focalDistance 2.360642";
//out << "    position 0.90718871 1.3196561 0.6429531 orientation 0.17981249 0.73576158 0.6529336  2.838383 nearDistance 0.5975644 farDistance 1.233297 focalDistance 0.91511345";




//	out << "}\n";
//	out <<  "Separator { Coordinate3 {\npoint [ .5 .3 1]} PointSet{}}\n";
		out << "ShapeHints { vertexOrdering COUNTERCLOCKWISE shapeType SOLID}";



//		out <<  "Separator { Coordinate3 {\npoint [";
//		fmatrix primitives = automate.get_ftrans("P");
//		std::cout << primitives << std::endl;
//		for (int j=0; j<primitives.ncols(); j++) {
//			out << primitives[0][j] << " " << primitives[1][j] << " " <<primitives[2][j] << ", ";
//		}
//		out << "]} PointSet{}}\n";
	} else {
		out << "%!PS-Adobe-3.0\n%%BoundingBox: 0 0 " << SIZE_X << " " << SIZE_Y << "\n%Pages: (atend)\n/circ{2 0 360 arc closepath fill}def /l{lineto}def /m{moveto}def /rm{rmoveto}def /rl{rlineto}def /S{stroke}def /CL{closepath}def /RGB{setrgbcolor}def /BLACK{0.1 0.1 0 RGB}def /GREEN{0 0 0 RGB}def /BLUE{0.174 1 1 RGB}def /RED{1 0 0 RGB}def\n\n%Page: 1\n0.1 setlinewidth\n";
		fmatrix primitives = automate.get_ftrans("P");
		for (int j=0; j<primitives.ncols(); j++) {
			out << primitives[0][j]*SIZE_X << " " << primitives[1][j]*SIZE_Y << " circ\n";
		}
// ajout CG on met le polygone de controle en rouge		
out << "RED\n" ;

		for (int j=0; j<primitives.ncols(); j++)
			out << primitives[0][j]*SIZE_X << " " << primitives[1][j]*SIZE_Y << " " << ((0==j) ? "m" : "l") << " ";
		if (1 != primitives.ncols()) out << " S\n";
// ajout CG on remet en noir
out << "BLACK\n" ;

	}

	fmatrix f;
	f.identity(automate.get_primitive(automate.get_initial_state(), 0).nrows());
	iterate(automate,f, automate.get_initial_state(), 0, out);

	if (VRML)
		out << "}\n";
	else
		out << "showpage\n%Trailer\n%Pages: 1\n";
	out.close();

	std::stringstream buf2;
	buf2 << (VRML ? "ivview " :"evince ") << buf.str();
	system(buf2.str().c_str());

	return 0;
}

