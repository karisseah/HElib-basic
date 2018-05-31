//
// Created by junjie on 5/14/18.
//

#include <gtest/gtest.h>
#include <iostream>
#include <random>
#include <FHE.h>
#include <EncryptedArray.h>
#include <iomanip>
#include "powerful.h"

using namespace std;
using namespace NTL;
using namespace chrono;

typedef high_resolution_clock Clock;

TEST(BootstrappingTest, BootstrappingTest_init)
{
	cout << "Hello, World!" << endl;

	cout << "-------------------- Initialization --------------------" << endl;
	auto begin_init = Clock::now();

	static long mValues[][14] = {
	  //{ p, phi(m),  m,    d, m1,  m2, m3,   g1,    g2,    g3,ord1,ord2,ord3, c_m}
		{  2,    48,   105, 12,  3,  35,  0,    71,    76,    0,  2,  2,   0, 100},
		{  2,   600,  1023, 10, 11,  93,  0,   838,   584,    0, 10,  6,   0, 100}, // m=(3)*11*{31} m/phim(m)=1.7    C=24  D=2 E=1
		{  2,  1200,  1705, 20, 11, 155,  0,   156,   936,    0, 10,  6,   0, 100}, // m=(5)*11*{31} m/phim(m)=1.42   C=34  D=2 E=2
		{  2,  1728,  4095, 12,  7,  5, 117,  2341,  3277, 3641,  6,  4,   6, 100}, // m=(3^2)*5*7*{13} m/phim(m)=2.36 C=26 D=3 E=2
		{  2,  2304,  4641, 24,  7,  3, 221,  3979,  3095, 3760,  6,  2,  -8, 100}, // m=3*7*(13)*{17} :-( m/phim(m)=2.01 C=45 D=4 E=3
		{  2,  4096,  4369, 16, 17, 257,  0,   258,  4115,    0, 16,-16,   0, 100}, // m=17*(257) :-( m/phim(m)=1.06 C=61 D=3 E=4
		{  2, 12800, 17425, 40, 41, 425,  0,  5951,  8078,    0, 40, -8,   0, 100}, // m=(5^2)*{17}*41 m/phim(m)=1.36 C=93  D=3 E=3
		{  2, 15004, 15709, 22, 23, 683,  0,  4099, 13663,    0, 22, 31,   0, 100}, // m=23*(683) m/phim(m)=1.04      C=73  D=2 E=1
		{  2, 16384, 21845, 16, 17,   5,257,  8996, 17477, 21591, 16, 4, -16,1600}, // m=5*17*(257) :-( m/phim(m)=1.33 C=65 D=4 E=4
		{  2, 18000, 18631, 25, 31, 601,  0, 15627,  1334,    0, 30, 24,   0,  50}, // m=31*(601) m/phim(m)=1.03      C=77  D=2 E=0
		{  2, 18816, 24295, 28, 43, 565,  0, 16386, 16427,    0, 42, 16,   0, 100}, // m=(5)*43*{113} m/phim(m)=1.29  C=84  D=2 E=2
		{  2, 21168, 27305, 28, 43, 635,  0, 10796, 26059,    0, 42, 18,   0, 100}, // m=(5)*43*{127} m/phim(m)=1.28  C=86  D=2 E=2
		{  2, 23040, 28679, 24, 17,  7, 241, 15184,  4098,28204, 16,  6, -10,1000}, // m=7*17*(241) m/phim(m)=1.24    C=63  D=4 E=3
		{  2, 24000, 31775, 20, 41, 775,  0,  6976, 24806,    0, 40, 30,   0, 100}, // m=(5^2)*{31}*41 m/phim(m)=1.32 C=88  D=2 E=2
		{  2, 26400, 27311, 55, 31, 881,  0, 21145,  1830,    0, 30, 16,   0, 100}, // m=31*(881) m/phim(m)=1.03      C=99  D=2 E=0
		{  2, 27000, 32767, 15, 31,  7, 151, 11628, 28087,25824, 30,  6, -10, 150},
		{  2, 31104, 35113, 36, 37, 949,  0, 16134,  8548,    0, 36, 24,   0, 400}, // m=(13)*37*{73} m/phim(m)=1.12  C=94  D=2 E=2
		{  2, 34848, 45655, 44, 23, 1985, 0, 33746, 27831,    0, 22, 36,   0, 100}, // m=(5)*23*{397} m/phim(m)=1.31  C=100 D=2 E=2
		{  2, 42336, 42799, 21, 127, 337, 0, 25276, 40133,    0,126, 16,   0,  20}, // m=127*(337) m/phim(m)=1.01     C=161 D=2 E=0
		{  2, 45360, 46063, 45, 73, 631,  0, 35337, 20222,    0, 72, 14,   0, 100}, // m=73*(631) m/phim(m)=1.01      C=129 D=2 E=0
		{  2, 46080, 53261, 24, 17, 13, 241, 43863, 28680,15913, 16, 12, -10, 100}, // m=13*17*(241) m/phim(m)=1.15   C=69  D=4 E=3
		{  2, 49500, 49981, 30, 151, 331, 0,  6952, 28540,    0,150, 11,   0, 100}, // m=151*(331) m/phim(m)=1        C=189 D=2 E=1
		{  2, 54000, 55831, 25, 31, 1801, 0, 19812, 50593,    0, 30, 72,   0, 100}, // m=31*(1801) m/phim(m)=1.03     C=125 D=2 E=0
		{  2, 60016, 60787, 22, 89, 683,  0,  2050, 58741,    0, 88, 31,   0, 200}, // m=89*(683) m/phim(m)=1.01      C=139 D=2 E=1
		// no. 24
		{  7,    36,    57,  3,  3,  19,  0,    20,    40,    0,  2, -6,   0, 100}, // m=3*(19) :-( m/phim(m)=1.58 C=14 D=3 E=0
		// no. 25
		{ 17,    48,   105, 12,  3,  35,  0,    71,    76,    0,  2,  2,   0, 100}, // m=3*(5)*{7} m/phim(m)=2.18 C=14 D=2 E=2
		{ 17,   576,  1365, 12,  7,   3, 65,   976,   911,  463,  6,  2,   4, 100}, // m=3*(5)*7*{13} m/phim(m)=2.36  C=22  D=3
		{ 17, 18000, 21917, 30, 101, 217, 0,  5860,  5455,    0, 100, 6,   0, 100}, // m=(7)*{31}*101 m/phim(m)=1.21  C=134 D=2
		{ 17, 30000, 34441, 30, 101, 341, 0,  2729, 31715,    0, 100, 10,  0, 100}, // m=(11)*{31}*101 m/phim(m)=1.14 C=138 D=2
		{ 17, 40000, 45551, 40, 101, 451, 0, 19394,  7677,    0, 100, 10,  0,2000}, // m=(11)*{41}*101 m/phim(m)=1.13 C=148 D=2
		{ 17, 46656, 52429, 36, 109, 481, 0, 46658,  5778,    0, 108, 12,  0, 100}, // m=(13)*{37}*109 m/phim(m)=1.12 C=154 D=2
		{ 17, 54208, 59363, 44, 23, 2581, 0, 25811,  5199,    0, 22, 56,   0, 100}, // m=23*(29)*{89} m/phim(m)=1.09  C=120 D=2
		{ 17, 70000, 78881, 10, 101, 781, 0, 67167, 58581,    0, 100, 70,  0, 100}, // m=(11)*{71}*101 m/phim(m)=1.12 C=178 D=2
		// no. 33
		{127,   576,  1365, 12,  7,   3, 65,   976,   911,  463,  6,  2,   4, 100}, // m=3*(5)*7*{13} m/phim(m)=2.36   C=22  D=3
		{127,  1200,  1925, 20,  11, 175, 0,  1751,   199,    0, 10, 6,    0, 100}, //  m=(5^2)*{7}*11 m/phim(m)=1.6   C=34 D=2
		{127,  2160,  2821, 30,  13, 217, 0,   652,   222,    0, 12, 6,    0, 100}, // m=(7)*13*{31} m/phim(m)=1.3     C=46 D=2
		{127, 18816, 24295, 28, 43, 565,  0, 16386, 16427,    0, 42, 16,   0, 100}, // m=(5)*43*{113} m/phim(m)=1.29   C=84  D=2
		{127, 26112, 30277, 24, 17, 1781, 0, 14249, 10694,    0, 16, 68,   0, 100}, // m=(13)*17*{137} m/phim(m)=1.15  C=106 D=2
		{127, 31752, 32551, 14, 43,  757, 0,  7571, 28768,    0, 42, 54,   0, 100}, // m=43*(757) :-( m/phim(m)=1.02   C=161 D=3
		{127, 46656, 51319, 36, 37, 1387, 0, 48546, 24976,    0, 36, -36,  0, 200}, //m=(19)*37*{73}:-( m/phim(m)=1.09 C=141 D=3
		{127, 49392, 61103, 28, 43, 1421, 0,  1422, 14234,    0, 42, 42,   0,4000}, // m=(7^2)*{29}*43 m/phim(m)=1.23  C=110 D=2
		{127, 54400, 61787, 40, 41, 1507, 0, 30141, 46782,    0, 40, 34,   0, 100}, // m=(11)*41*{137} m/phim(m)=1.13  C=112 D=2
		{127, 72000, 77531, 30, 61, 1271, 0,  7627, 34344,    0, 60, 40,   0, 100}  // m=(31)*{41}*61 m/phim(m)=1.07   C=128 D=2
	};

	long m = 105;                                       // Specific modulus
	long p = 2;                                         // Plaintext base
	long r = 1;                                         // Lifting
	long L = 25;                                        // Number of levels in the modulus chain
	long c = 1;                                         // Number of columns in key-switching matrix
	long w = 2;                                         // Hamming weight of secret key
	long d = 1;                                         // Degree of the field extension
	long k = 80;                                        // Security parameter
	long s = 0;                                         // Minimum number of slots

	Vec<long> mvec;
	vector<long> gens;
	vector<long> ords;

	int idx = 0;
	append(mvec, mValues[idx][4]);
	if (mValues[idx][5]>1) append(mvec, mValues[idx][5]);
	if (mValues[idx][6]>1) append(mvec, mValues[idx][6]);
	gens.push_back(mValues[idx][7]);
	if (mValues[idx][8]>1) gens.push_back(mValues[idx][8]);
	if (mValues[idx][9]>1) gens.push_back(mValues[idx][9]);
	ords.push_back(mValues[idx][10]);
	if (abs(mValues[idx][11])>1) ords.push_back(mValues[idx][11]);
	if (abs(mValues[idx][12])>1) ords.push_back(mValues[idx][12]);
//	long one = 1;
//	long three = 3;
//	long five = 5;
//	long seven = 7;
//	long threefive = 35;
//	append(mvec, three);
//	append(mvec, threefive);

	// Finding m
//	m = FindM(k, L, c, p, d, s, 0);                     // Find a value for m given the params
	cout << "m = " << m << endl;
	cout << "p = " << p << endl;

	// Initializing context
	FHEcontext context(m, p, r, gens, ords);                        // Initialize context
	buildModChain(context, L, c);                       // Modify the context, adding primes to the modulus chain
	context.makeBootstrappable(mvec);

	// Creating polynomial
	ZZX G =  context.alMod.getFactorsOverZZ()[0];       // Creates the polynomial for a plaintext slot
	cout << "degree: " << deg(G) << endl;
	cout << "G: " << G << endl;

	// Generating keys
	FHESecKey secretKey(context);                       // Construct a secret key structure
	FHEPubKey& publicKey = secretKey;             // An "upcast": FHESecKey is a subclass of FHEPubKey - Creates publicKey from secretKey
	secretKey.GenSecKey(w);                             // Generate a secret key with Hamming weight w
	addSome1DMatrices(secretKey);                       // apparently for key switching
	addFrbMatrices(secretKey);                          // for Ctxt rotate

	// Bootstrapping
	secretKey.genRecryptData();

	// Helper Class for encryption and decryption
	EncryptedArray ea(context, G);                      // Construct EncryptedArray object ea, associated with context and G
	long nslots = ea.size();                            // Number of slots in the plaintext encoding
	cout << "slots: " << nslots << endl;

	auto end_init = Clock::now();
	cout << "FHE Ready!" << endl;
	cout << "It took: " << duration_cast<seconds>(end_init - begin_init).count() << " seconds." << endl;


	cout << "-------------------- Encryption --------------------" << endl;
	auto begin_encrypt = Clock::now();
	ZZX ptxt_poly;
	SetCoeff(ptxt_poly,1);

	cout << "ptxt_poly: " << ptxt_poly << endl;

	Ctxt enc(publicKey);
	publicKey.Encrypt(enc,ptxt_poly);
	cout << "stats: " << enc.findBaseLevel() << ", " << enc.getNoiseVar() <<endl;

	auto end_encrypt = Clock::now();
	cout << "Encryption Over!" << endl;
	cout << "It took: " << duration_cast<seconds>(end_encrypt - begin_encrypt).count() << " seconds." << endl;


	cout << "-------------------- Operation --------------------" << endl;
	Ctxt temp = enc;
	int counter = 0;

	for (int i=0 ; i<25 ; i++)
	{
		enc.multiplyBy(temp);
		counter++;
		cout << i << ": " << enc.findBaseLevel() << ", " << enc.getNoiseVar() <<endl;
	}

	publicKey.reCrypt(enc);
	cout << "next : " << enc.findBaseLevel() << ", " << enc.getNoiseVar() <<endl;

	for (int i=0 ; i<15 ; i++)
	{
		enc.multiplyBy(temp);
		counter++;
		cout << i << ": " << enc.findBaseLevel() << ", " << enc.getNoiseVar() <<endl;
	}

	publicKey.reCrypt(enc);
	cout << "next : " << enc.findBaseLevel() << ", " << enc.getNoiseVar() <<endl;

	for (int i=0 ; i<6 ; i++)
	{
		enc.multiplyBy(temp);
		counter++;
		cout << i << ": " << enc.findBaseLevel() << ", " << enc.getNoiseVar() <<endl;
	}


	cout << "-------------------- Decryption --------------------" << endl;
	ZZX ans;
	secretKey.Decrypt(ans,enc);
	cout << "Plaintext:  " << ans << endl;
	cout << "Multiplication Count: " << counter << endl;
}