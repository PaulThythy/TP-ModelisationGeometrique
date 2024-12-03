#include <iostream>
#include <fstream>
#include <sstream>

#include "image.h"

//bbox noinit
#define NOINIT 65535

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Image::Image(unsigned short w, unsigned short h, unsigned char init) :  width(w), height(h) {
	buffer = new unsigned char[w*h];
	for (int i=0; i<w*h; i++)
		buffer[i] = init;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Image::Image() : buffer(NULL), width(0), height(0) {}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Image::Image(char *filename) : buffer(NULL), width(0), height(0) {
	std::ifstream in;
	in.open (filename, std::ios::binary);

	// get length of file:
	in.seekg (0, std::ios::end);
	long length = in.tellg();
	in.seekg (0, std::ios::beg);

	unsigned short head_cnt = 0;
	char line[256];
	while (head_cnt<3 && !in.eof()) {
		in.getline(line, 256);
		if ('#'==*line) continue;
		if (2==++head_cnt) {
			std::istringstream dim(line);
			dim >> width >> height;
		}
	}

	if (length - in.tellg() != width*height) {
		throw ("Bad file");
	}

	buffer = new unsigned char[width*height];
	in.read ((char *)buffer,width*height);
	in.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Image & Image::operator =(const Image& img) {
	if (this != &img) {
		width  = img.width;
		height = img.height;
		long size = width*height;
		if (buffer) delete [] buffer;
		buffer = new unsigned char[size];
		for (long i=0; i<size; i++)
			buffer[i] = img.buffer[i];
	}
	return *this;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Image::Image(const Image &img) : width(img.width), height(img.height) {
	long size = width*height;
	buffer = new unsigned char[size];

	for (long i=0; i<size; i++)
		buffer[i] = img.buffer[i];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Image::~Image() {
	if (buffer) delete [] buffer;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

unsigned short Image::get_width() {
	return width;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

unsigned short Image::get_height() {
	return height;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline void Image::set(unsigned short x, unsigned short y, unsigned char val) {
	if (x>=width || y>=height) {
//		std::cerr << x << ">" << width << " " << y << ">" << height << "\n";
	} else {
		buffer[x+y*width] = val;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline unsigned char Image::get(unsigned short x, unsigned short y) {
	return buffer[x+y*width];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Image::dump(char *filename) {
	std::ofstream out;
	out.open (filename, std::ios::binary);
	out << "P5\n" << width << " " << height << "\n" << 255 << "\n";
	out.write((char *)buffer, width*height);
	out.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Image::scale(unsigned short new_width, unsigned short new_height, unsigned char threshold) {
	// horizontal scale
	unsigned char *temp_buffer = new unsigned char[new_width*height];
	float ratio = (float)width/(float)new_width;
	for (unsigned short j=0; j<height; j++) {
		long offset     = j*width;
		long new_offset = j*new_width;
		for (unsigned short i=0; i<new_width; i++) {
			unsigned short start = (unsigned short)(i*ratio);
			unsigned short end   = (unsigned short)((i+1)*ratio);
			long sum = 0;
			unsigned short cnt = 0;
			while (start+cnt<=end && start+cnt<width) sum += buffer[offset + start+cnt++];
			// no halftones
			temp_buffer[i+new_offset] = (sum/cnt>threshold ? 255 : 0);
		}
	}

	if (buffer) delete [] buffer;
	buffer = new unsigned char[new_width*new_height];
	// vertical scale
	ratio = (float)height/(float)new_height;
	for (unsigned short i=0; i<new_width; i++) {
		for (unsigned short j=0; j<new_height; j++) {
			unsigned short start = (unsigned short)(j*ratio);
			unsigned short end   = (unsigned short)((j+1)*ratio);
			long sum = 0;
			unsigned short cnt = 0;
			while (start+cnt<=end && start+cnt<height) sum += temp_buffer[i+(start+cnt++)*new_width];
			// no halftones again
			buffer[i+j*new_width] = (sum/cnt>threshold ? 255 : 0);
		}
	}

	delete [] temp_buffer;
	width  = new_width;
	height = new_height;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Image::median(Image &img2, Image &median) {
	if (median.width!=width || median.height!=height) {
		if (median.buffer) delete [] median.buffer;
		median.buffer = new unsigned char[width*height];
		median.width  = width;
		median.height = height;
	}
	long size = width*height;
	for (long i=0; i<size; i++) median.buffer[i] = 255;

	unsigned short x1, y1, w1, h1, x2, y2, w2, h2, mx, my, mw, mh;

	bbox(x1, y1, w1, h1);
	img2.bbox(x2, y2, w2, h2);
	if (x1==NOINIT||y1==NOINIT) {x1=x2; y1=y2; w1=w2; h1=h2;}

	Image sub1, sub2;
	sub_rect(x1, y1, w1, h1, sub1);
	img2.sub_rect(x2, y2, w2, h2, sub2);

	mw = w1 + w2 >> 1;
	mh = h1 + h2 >> 1;
	mx = x1 + x2 >> 1;
	my = y1 + y2 >> 1;
	sub1.scale(mw, mh);
	sub2.scale(mw, mh);

	size = mw*mh;

	unsigned short *nearx1 = new unsigned short[size];
	unsigned short *neary1 = new unsigned short[size];
	unsigned short *nearx2 = new unsigned short[size];
	unsigned short *neary2 = new unsigned short[size];

	Image sm(mw, mh);
	for (long i=0; i<size; i++) sm.buffer[i] = 255;
	sub1.feature_transform(nearx1, neary1);
	sub2.feature_transform(nearx2, neary2);
	for (unsigned short y=0; y<mh; y++) {
		for (unsigned short x=0; x<mw; x++) {
			if (0==sub2.get(x, y)) {
				unsigned short xx = (unsigned short)(x + 0.8*(nearx1[y*mw+x]-x));
				unsigned short yy = (unsigned short)(y + 0.8*(neary1[y*mw+x]-y));
				sm.set(xx, yy, 0);
			}
			if (0==sub1.get(x, y)) {
				unsigned short xx = (unsigned short)(x + 0.8*(nearx2[y*mw+x]-x));
				unsigned short yy = (unsigned short)(y + 0.8*(neary2[y*mw+x]-y));
				sm.set(xx, yy, 0);
			}
		}
	}

//	sm.dump("out/sm.pgm");
	median.redraw_rect(sm, mx, my);

	delete [] nearx1;
	delete [] neary1;
	delete [] nearx2;
	delete [] neary2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Image::bbox(unsigned short &x, unsigned short &y, unsigned short &w, unsigned short &h) {
	x = y = w = h = NOINIT;
	long offset = 0;
	for (unsigned short j=0; j<height; j++) {
		for (unsigned short i=0; i<width; i++) {
			if (255!=buffer[i+offset]) {
				if (x>i || NOINIT==x) x = i;
				if (y>j || NOINIT==y) y = j;
				if (w<i || NOINIT==w) w = i;
				if (h<j || NOINIT==h) h = j;
			}
		}
		offset += width;
	}
	w -= x-1;
	h -= y-1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Image::sub_rect(unsigned short x, unsigned short y, unsigned short w, unsigned short h, Image &result) {
	if (result.width!=w || result.height!=h) {
		if (result.buffer) delete [] result.buffer;
		result.buffer = new unsigned char[width*height];
		result.width  = w;
		result.height = h;
	}

	for (unsigned short j=0; j<h; j++)
		for (unsigned short i=0; i<w; i++)
			result.buffer[i+j*w] = get(i+x, j+y);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Image::redraw_rect(Image &src, unsigned short x, unsigned short y) {
	for (unsigned short j=0; j<src.height; j++)
		for (unsigned short i=0; i<src.width; i++)
			buffer[i+x+(j+y)*width] = src.get(i, j);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Image::clear_rect(unsigned short x, unsigned short y, unsigned short w, unsigned short h, unsigned char val) {
	for (unsigned short j=0; j<h; j++)
		for (unsigned short i=0; i<w; i++)
			buffer[i+x+(j+y)*width] = val;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Image::unite(Image &img) {
	long size = height*width;
	for (long i=0; i<size; i++)
		if (0==img.buffer[i])
			buffer[i] = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define A 10
#define B 14
#define INFTY 0xFFFF

void Image::feature_transform(unsigned short *nearx, unsigned short *neary) {
	long size = width*height;
	long *dist = new long[size];
	unsigned short x, y;
	x = y = 0;
	for (long i=0; i<size; i++) {
		dist[i]  = INFTY;
		nearx[i] = neary[i] = 0;
		if (0==buffer[i]) {
			dist[i]  = 0;
			nearx[i] = x;
			neary[i] = y;
		}
		if (++x==width) x = 0, y++;
	}

	x = y = 0;
	long min, minidx;
	for (long i=0; i<size; i++) {
		min = INFTY;
		if (y>0              && min>dist[i-width]  +A) {min=dist[i-width]+A;  minidx=i-width;}
		if (y>0 && x>0       && min>dist[i-width-1]+B) {min=dist[i-width-1]+B,minidx=i-width-1;}
		if (y>0 && x<width-1 && min>dist[i-width+1]+B) {min=dist[i-width+1]+B;minidx=i-width+1;}
		if (       x>0       && min>dist[i-1]      +A) {min=dist[i-1]    +A;  minidx=i-1;}
		if (dist[i]>min) {
			dist[i] = min;
			nearx[i] = nearx[minidx];
			neary[i] = neary[minidx];
		}
		if (++x==width) {x=0; y++;}
	}

	x = width-1;
	y = height-1;
	for (long i=size-1; i>=0; i--) {
		min = INFTY;
		if (y<height-1              && min>dist[i+width]  +A) {min=dist[i+width]  +A; minidx=i+width;}
		if (y<height-1 && x>0       && min>dist[i+width-1]+B) {min=dist[i+width-1]+B; minidx=i+width-1;}
		if (y<height-1 && x<width-1 && min>dist[i+width+1]+B) {min=dist[i+width+1]+B; minidx=i+width+1;}
		if (              x<width-1 && min>dist[i+1]      +A) {min=dist[i+1]      +A; minidx=i+1;}
		if (dist[i]>min) {
			dist[i]  = min;
			nearx[i] = nearx[minidx];
			neary[i] = neary[minidx];
		}
		if (x--==0) {x = width-1; y--;}
	}

	delete [] dist;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Image::inverse() {
	long size = height*width;
	for (long i=0; i<size; i++)
		buffer[i] = 255-buffer[i];
}

