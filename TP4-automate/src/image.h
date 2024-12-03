class Image {
private:
	unsigned char *buffer;
	unsigned short width, height;

public:
	Image();
	Image(char *filename);
	Image(unsigned short w, unsigned short h, unsigned char init=0);
	~Image();

	Image(const Image &i);
	Image &operator =(const Image& img);

	void scale(unsigned short width, unsigned short height, unsigned char threshold=128);
	unsigned short get_width();
	unsigned short get_height();
	unsigned char get(unsigned short x, unsigned short y);
	void set(unsigned short x, unsigned short y, unsigned char val);
	void inverse();
	void dump(char *filename);

	void unite(Image &img);
	void bbox(unsigned short &x, unsigned short &y, unsigned short &w, unsigned short &h);
	void sub_rect(unsigned short x, unsigned short y, unsigned short w, unsigned short h, Image &result);
	void redraw_rect(Image &src, unsigned short x, unsigned short y);
	void clear_rect(unsigned short x, unsigned short y, unsigned short w, unsigned short h, unsigned char val=255);

	void median(Image &img2, Image &median);
	void feature_transform(unsigned short *nearx, unsigned short *neary);

};

typedef Image * imgPtr;
