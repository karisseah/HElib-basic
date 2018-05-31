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

int main()
{
	cout << "Hello, World!" << endl;

	cout << "-------------------- Initialization --------------------" << endl;
	auto begin_init = Clock::now();

	long m = 8;                                         // Specific modulus
	long p = 17;                                        // Plaintext base
	long r = 1;                                         // Lifting
	long L = 5;                                         // Number of levels in the modulus chain
	long c = 1;                                         // Number of columns in key-switching matrix
	long w = 2;                                         // Hamming weight of secret key
	long d = 1;                                         // Degree of the field extension
	long k = 80;                                        // Security parameter
	long s = 0;                                         // Minimum number of slots

	// Finding m
//	m = FindM(k, L, c, p, d, s, 0);                     // Find a value for m given the params
	cout << "m = " << m << endl;

	// Initializing context
	FHEcontext context(m, p, r);                        // Initialize context
	buildModChain(context, L, c);                       // Modify the context, adding primes to the modulus chain

	// Creating polynomial
	ZZX G = context.alMod.getFactorsOverZZ()[0];       // Creates the polynomial for a plaintext slot
	cout << "degree: " << deg(G) << endl;
	cout << "G: " << G << endl;

	// Generating keys
	FHESecKey secretKey(context);                       // Construct a secret key structure
	FHEPubKey &publicKey = secretKey;                   // An "upcast": FHESecKey is a subclass of FHEPubKey - Creates publicKey from secretKey
	secretKey.GenSecKey(w);                             // Generate a secret key with Hamming weight w
	addSome1DMatrices(secretKey);                       // apparently for key switching
	addFrbMatrices(secretKey);                          // for Ctxt rotate

	// Helper Class for encryption and decryption
	EncryptedArray ea(context, G);                      // Construct EncryptedArray object ea, associated with context and G
	long nslots = ea.size();                            // Number of slots in the plaintext encoding
	cout << "slots: " << nslots << endl;

	auto end_init = Clock::now();
	cout << "FHE Ready!" << endl;
	cout << "It took: " << duration_cast<seconds>(end_init - begin_init).count() << " seconds." << endl;


	cout << "-------------------- Encryption --------------------" << endl;
	auto begin_encrypt = Clock::now();
	zz_pX poly_p;
	SetCoeff(poly_p, 1);
	ZZX ptxt_poly = conv<ZZX>(poly_p);
	PolyRed(ptxt_poly, 17, true); // reduce to the symmetric interval

	cout << "ptxt_poly: " << ptxt_poly << endl;

	Ctxt enc(publicKey);
	publicKey.Encrypt(enc, ptxt_poly);

	auto end_encrypt = Clock::now();
	cout << "Encryption Over!" << endl;
	cout << "It took: " << duration_cast<seconds>(end_encrypt - begin_encrypt).count() << " seconds." << endl;


	cout << "-------------------- Operation --------------------" << endl;
	Ctxt temp = enc;
	for (int i=0 ; i<4 ; i++)
	{
		enc.multiplyBy(temp);
	}


	cout << "-------------------- Decryption --------------------" << endl;
	ZZX ans;
	secretKey.Decrypt(ans, enc);
	cout << "Plaintext:  " << ans << endl;

	return 0;
}