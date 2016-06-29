#ifndef __LF_DEM__MatrixBlocks__
#define __LF_DEM__MatrixBlocks__
#include <array>

struct ODBlock {
	std::array<double, 5> col0;
	std::array<double, 3> col1;
	std::array<double, 1> col2;
	std::array<double, 5> col3;
	std::array<double, 3> col4;
	std::array<double, 1> col5;
	int bla;
};

inline void resetODBlock(struct ODBlock &b)
{
	b.col0.fill(0);
	b.col1.fill(0);
	b.col2.fill(0);
	b.col3.fill(0);
	b.col4.fill(0);
	b.col5.fill(0);
}

struct DBlock {
	std::array<double, 5> col0;
	std::array<double, 3> col1;
	std::array<double, 1> col2;
	std::array<double, 3> col3;
	std::array<double, 2> col4;
	std::array<double, 1> col5;

	inline DBlock& operator += (const DBlock &b)
	{
		for (unsigned int i=0; i<col0.size(); i++) {
			col0[i] += b.col0[i];
		}
		for (unsigned int i=0; i<col1.size(); i++) {
			col1[i] += b.col1[i];
		}
		col2[0] += b.col2[0];
		for (unsigned int i=0; i<col3.size(); i++) {
			col3[i] += b.col3[i];
		}
		for (unsigned int i=0; i<col4.size(); i++) {
			col4[i] += b.col4[i];
		}
		col5[0] += b.col5[0];
		return *this;
	}
};

inline void resetDBlock(struct DBlock& b)
{
	b.col0.fill(0);
	b.col1.fill(0);
	b.col2.fill(0);
	b.col3.fill(0);
	b.col4.fill(0);
	b.col5.fill(0);
}

inline void addDBlock(struct DBlock& b)
{
	b.col0.fill(0);
	b.col1.fill(0);
	b.col2.fill(0);
	b.col3.fill(0);
	b.col4.fill(0);
	b.col5.fill(0);
}

#endif
