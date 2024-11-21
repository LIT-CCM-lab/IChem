/* ****************************************************************************
 *
 * Copyright 2013 informaticien77
 *
 * This file is part of EasyCrypt - the AES and RSA implementation in C++.
 *
 * EasyCrypt is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * EasyCrypt is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rsa.  If not, see <http://www.gnu.org/licenses/>.
 *
 * 				EasyCrypt.h
 *
 * Author: informaticien77
 * Release date: 22th of May 2013
 *
 * This file contains the declaration for the EasyCrypt class.
 *
 * ****************************************************************************
 */

#ifndef EASYCRYPT_H
#define EASYCRYPT_H

#include <string>
#include <fstream>

namespace ec
{
	class AES
	{
	public:

		static std::string GenerateKey(int keybits = 128);

		static std::string Encrypt(std::string plaintext, std::string key, int keybits = 128);

		static std::string Decrypt(std::string ciphertext, std::string key, int keybits = 128);

	};
	class RSA
	{
	public:

		static bool GenerateKeys(unsigned long int digitCount, std::string *moduluskeys, std::string *privatexponent, std::string *publicexponent, unsigned long int precisionOfPrimality = 3);

		static std::string Encrypt(std::string plaintext, std::string keymod, std::string keyexp);

		static std::string Decrypt(std::string ciphertext, std::string keymod, std::string keyexp);

	};
}

//#define KEYLENGTH(keybits) ((keybits)/8)
//#define RKLENGTH(keybits)  ((keybits)/8+28)
//#define NROUNDS(keybits)   ((keybits)/32+6)

#endif
