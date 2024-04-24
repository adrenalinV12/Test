
Caesar cipher

#include <iostream>
#include <cctype>
using namespace std;

void encrypt(string& plaintext, const string& key) {
    int key_index = 0;
    for (char& c : plaintext) {
        if (isalpha(c)) {
            int shift = toupper(key[key_index++ % key.length()]) - 'A';
            c = ((((c - 'A') + shift) % 26) + 'A');
        }
    }
}

void decrypt(string& ciphertext, const string& key) {
    int key_index = 0;
    for (char& c : ciphertext) {
        if (isalpha(c)) {
            int shift = toupper(key[key_index++ % key.length()]) - 'A';
            c = ((((c - 'A') - shift + 26) % 26) + 'A');
        }
    }
}

int main() {
    string plaintext = "SHARUQ IS SLEEPING";
    string key = "SCOPE";
    cout << "Plaintext: " << plaintext << endl;
    encrypt(plaintext, key);
    cout << "Ciphertext: " << plaintext << endl;
    decrypt(plaintext, key);
    cout << "Decrypted text: " << plaintext << endl;
    return 0;
}



Hill cipher

#include <iostream>
#include <vector>
#include <random>
using namespace std;

vector<int> matrix_multiply(const vector<vector<int>>& keyMatrix, const vector<int>& textVector) {
    vector<int> result(2, 0);
    for (int i = 0; i < 2; i++) {
        for (int k = 0; k < 2; k++) {
            result[i] += keyMatrix[i][k] * textVector[k];
        }
        result[i] %= 26;
    }
    return result;
}

int main() {
    vector<vector<int>> keyMatrix = {{3, 17}, {6, 5}};
    string plaintext = "GOOD";
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(1, 25);
    int randomOffset = dis(gen);
    vector<vector<int>> plaintextMatrix;

    // Split plaintext into pairs of characters and convert to numbers
    for (int i = 0; i < plaintext.length(); i += 2) {
        vector<int> pair;
        for (int j = 0; j < 2 && i + j < plaintext.length(); j++) {
            pair.push_back(plaintext[i + j] - 'A');
        }
        plaintextMatrix.push_back(pair);
    }

    // Encrypt each pair of plaintext
    string ciphertext = "";
    for (const vector<int>& pair : plaintextMatrix) {
        vector<int> encryptedPair = matrix_multiply(keyMatrix, pair);
        for (int& num : encryptedPair) {
            num = (num + randomOffset) % 26;
            ciphertext += char('A' + num);
        }
    }
    cout << "Ciphertext: " << ciphertext << endl;
    return 0;
}



Vernam cipher

#include <iostream>
#include <string>
using namespace std;

string encrypt(const string& plaintext, const string& key) {
    string ciphertext = "";
    for (int i = 0; i < plaintext.length(); i++) {
        int p = plaintext[i] - 'A';
        int k = key[i % key.length()] - 'A';
        int c = (p + k) % 26;
        ciphertext += char('A' + c);
    }
    return ciphertext;
}

string decrypt(const string& ciphertext, const string& key) {
    string plaintext = "";
    for (int i = 0; i < ciphertext.length(); i++) {
        int c = ciphertext[i] - 'A';
        int k = key[i % key.length()] - 'A';
        int p = (c - k + 26) % 26;
        plaintext += char('A' + p);
    }
    return plaintext;
}

int main() {
    string plaintext = "SHANTI";
    string key = "VIT";
    string ciphertext = encrypt(plaintext, key);
    cout << "Ciphertext: " << ciphertext << endl;
    string decryptedText = decrypt(ciphertext, key);
    cout << "Decrypted text: " << decryptedText << endl;
    return 0;
}


RAIL CIPHER 

#include <iostream>
#include <string>

using namespace std;

string railEncrypt(string plaintext, string key) {
    int keyLength = key.length();
    string ciphertext = "";
    string rails[keyLength];

    int rail = 0;
    bool down = false;

    for (char& c : plaintext) {
        rails[rail] += c;

        if (rail == 0 || rail == keyLength - 1)
            down = !down;

        if (down)
            rail++;
        else
            rail--;
    }

    for (int i = 0; i < keyLength; i++) {
        ciphertext += rails[key.find(key[i])];
    }

    return ciphertext;
}

string railDecrypt(string ciphertext, string key) {
    int keyLength = key.length();
    string plaintext = "";
    string rails[keyLength];

    int rail = 0;
    bool down = false;

    // Calculate the number of characters in each rail
    int charCount[keyLength] = {0};
    for (int i = 0; i < ciphertext.length(); i++) {
        charCount[rail]++;
        if (rail == 0 || rail == keyLength - 1)
            down = !down;
        if (down)
            rail++;
        else
            rail--;
    }

    // Fill in the rails with ciphertext characters
    int index = 0;
    for (int i = 0; i < keyLength; i++) {
        for (int j = 0; j < charCount[key.find(key[i])]; j++) {
            rails[key.find(key[i])] += ciphertext[index++];
        }
    }

    // Reconstruct the plaintext
    rail = 0;
    down = false;
    for (int i = 0; i < ciphertext.length(); i++) {
        plaintext += rails[rail][0];
        rails[rail].erase(0, 1);
        if (rail == 0 || rail == keyLength - 1)
            down = !down;
        if (down)
            rail++;
        else
            rail--;
    }

    return plaintext;
}

int main() {
    string plaintext = "HELLO WORLD";
    string key = "SCOPE";

    string encrypted = railEncrypt(plaintext, key);
    cout << "Encrypted: " << encrypted << endl;

    string decrypted = railDecrypt(encrypted, key);
    cout << "Decrypted: " << decrypted << endl;

    return 0;
}


DSS 

#include <iostream>
using namespace std;

// Function to calculate modular exponentiation
int modExp(int base, int exponent, int modulus) {
    int result = 1;
    base = base % modulus;
    while (exponent > 0) {
        if (exponent % 2 == 1)
            result = (result * base) % modulus;
        exponent = exponent >> 1;
        base = (base * base) % modulus;
    }
    return result;
}

// Function to calculate modular inverse
int mod_inverse(int a, int m) {
    int m0, x0, x1;
    m0 = m;
    x0 = 0;
    x1 = 1;
    while (a > 1) {
        int q = a / m;
        int temp = m;
        m = a % m;
        a = temp;
        int temp_x = x0;
        x0 = x1 - q * x0;
        x1 = temp_x;
    }
    return (x1 < 0) ? x1 + m0 : x1;
}

int main() {
    int p, q, H, G, X, k;
    cout << "Give the value of p (prime 1): ";
    cin >> p;
    cout << "Give the value of q (prime 2): ";
    cin >> q;
    cout << "Give the value of H(M) (message, less than N): ";
    cin >> H;
    cout << "Give the value of G: ";
    cin >> G;
    cout << "Give the value of X (Private Key): ";
    cin >> X;

    // Calculate public key y
    int y = modExp(G, X, p);
    cout << "The value of y (public key) is: " << y << endl;

    cout << "Give the value of k chosen: ";
    cin >> k;

    // Calculate r
    int r = modExp(G, k, p) % q;
    cout << "The value of r is: " << r << endl;

    // Calculate s
    int w = 0;
    while ((w * k) % q != 1)
        w++;
    int s = (w * (H + X * r)) % q;
    cout << "The value of s is: " << s << endl;

    // Verification
    cout << "For verification:" << endl;
    int u1 = (H * w) % q;
    cout << "The value of w is: " << w << endl;
    cout << "The value of u1 is: " << u1 << endl;
    int u2 = (r * w) % q;
    cout << "The value of u2 is: " << u2 << endl;

    int v = ((modExp(G, u1, p) * modExp(y, u2, p)) % p) % q;

    cout << "Applying to checking formula:" << endl;
    cout << r << "==" << v << endl;

    if (r == v)
        cout << "Signature verified successfully." << endl;
    else
        cout << "Signature verification failed." << endl;

    return 0;
}


ELGAMAL 

#include<iostream>
#include<cmath>
using namespace std;


int power(int base, unsigned int exp, int mod) {
    int result = 1;
    base = base % mod;
    while (exp > 0) {
        if (exp % 2 == 1)
            result = (result * base) % mod;
        exp = exp >> 1;
        base = (base * base) % mod;
    }
    return result;
}

int main() {
    
    int q = 19;
    int alpha = 10;
    int Xa = 5;
    int k = 6;
    int M = 17;

    // Calculate Public Key Ya
    int Ya = power(alpha, Xa, q);

    // Encryption
    int C1 = power(alpha, k, q);
    int C2 = (M * power(Ya, k, q)) % q;

    // Output
    cout << "Encryption\n";
    cout << "C1 = " << C1 << endl;
    cout << "C2 = " << C2 << endl;

    // Decryption
    int K = power(C1, Xa, q);
    int M_decrypted = (C2 * power(K, q - 2, q)) % q;

    // Output
    cout << "\nDecryption\n";
    cout << "K = " << K << endl;
    cout << "M = " << M_decrypted << endl;

    return 0;
}	 	  	 	 	  	    	 	     	  	        	 	


DIFFIE.C

#include <stdio.h>
#include <math.h>

// Function to perform modular exponentiation
int mod_exp(int base, int exponent, int modulus) {
    int result = 1;
    base = base % modulus;
    while (exponent > 0) {
        if (exponent % 2 == 1)
            result = (result * base) % modulus;
        exponent = exponent >> 1;
        base = (base * base) % modulus;
    }
    return result;
}

// Function to perform Diffie-Hellman key exchange
void diffie_hellman(int q, int alpha, int xa, int xb, int *ka, int *kb) {
    // Calculate public keys
    int ya = mod_exp(alpha, xa, q);
    int yb = mod_exp(alpha, xb, q);
    
    // Calculate shared secret keys
    *ka = mod_exp(yb, xa, q);
    *kb = mod_exp(ya, xb, q);
}

int main() {
    int q, alpha, xa, xb, ka, kb;
    
    // Sample Input
    printf("Enter a Prime Number 'q': ");
    scanf("%d", &q);
    printf("Enter alpha: ");
    scanf("%d", &alpha);
    printf("Enter a No 'xa' which is less than value of q: ");
    scanf("%d", &xa);
    printf("Enter a No 'xb' which is less than value of q: ");
    scanf("%d", &xb);
    
    // Calculate keys
    diffie_hellman(q, alpha, xa, xb, &ka, &kb);
    
    // Sample Output
    printf("ka = %d\n", ka);
    printf("kb = %d\n", kb);
    
    return 0;
}	 	  	 	 	  	    	 	     	  	        	 	


SHA512

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>

// Function to perform SHA-512 hashing
std::string sha512(const std::string& message) {
    // Initialize SHA-512 constants
    const uint64_t k[80] = {
        0x428a2f98d728ae22ULL, 0x7137449123ef65cdULL, 0xb5c0fbcfec4d3b2fULL, 0xe9b5dba58189dbbcULL,
        0x3956c25bf348b538ULL, 0x59f111f1b605d019ULL, 0x923f82a4af194f9bULL, 0xab1c5ed5da6d8118ULL,
        0xd807aa98a3030242ULL, 0x12835b0145706fbeULL, 0x243185be4ee4b28cULL, 0x550c7dc3d5ffb4e2ULL,
        0x72be5d74f27b896fULL, 0x80deb1fe3b1696b1ULL, 0x9bdc06a725c71235ULL, 0xc19bf174cf692694ULL,
        0xe49b69c19ef14ad2ULL, 0xefbe4786384f25e3ULL, 0x0fc19dc68b8cd5b5ULL, 0x240ca1cc77ac9c65ULL,
        0x2de92c6f592b0275ULL, 0x4a7484aa6ea6e483ULL, 0x5cb0a9dcbd41fbd4ULL, 0x76f988da831153b5ULL,
        0x983e5152ee66dfabULL, 0xa831c66d2db43210ULL, 0xb00327c898fb213fULL, 0xbf597fc7beef0ee4ULL,
        0xc6e00bf33da88fc2ULL, 0xd5a79147930aa725ULL, 0x06ca6351e003826fULL, 0x142929670a0e6e70ULL,
        0x27b70a8546d22ffcULL, 0x2e1b21385c26c926ULL, 0x4d2c6dfc5ac42aedULL, 0x53380d139d95b3dfULL,
        0x650a73548baf63deULL, 0x766a0abb3c77b2a8ULL, 0x81c2c92e47edaee6ULL, 0x92722c851482353bULL,
        0xa2bfe8a14cf10364ULL, 0xa81a664bbc423001ULL, 0xc24b8b70d0f89791ULL, 0xc76c51a30654be30ULL,
        0xd192e819d6ef5218ULL, 0xd69906245565a910ULL, 0xf40e35855771202aULL, 0x106aa07032bbd1b8ULL,
        0x19a4c116b8d2d0c8ULL, 0x1e376c085141ab53ULL, 0x2748774cdf8eeb99ULL, 0x34b0bcb5e19b48a8ULL,
        0x391c0cb3c5c95a63ULL, 0x4ed8aa4ae3418acbULL, 0x5b9cca4f7763e373ULL, 0x682e6ff3d6b2b8a3ULL,
        0x748f82ee5defb2fcULL, 0x78a5636f43172f60ULL, 0x84c87814a1f0ab72ULL, 0x8cc702081a6439ecULL,
        0x90befffa23631e28ULL, 0xa4506cebde82bde9ULL, 0xbef9a3f7b2c67915ULL, 0xc67178f2e372532bULL,
        0xca273eceea26619cULL, 0xd186b8c721c0c207ULL, 0xeada7dd6cde0eb1eULL, 0xf57d4f7fee6ed178ULL,
        0x06f067aa72176fbaULL, 0x0a637dc5a2c898a6ULL, 0x113f9804bef90daeULL, 0x1b710b35131c471bULL,
        0x28db77f523047d84ULL, 0x32caab7b40c72493ULL, 0x3c9ebe0a15c9bebcULL, 0x431d67c49c100d4cULL,
        0x4cc5d4becb3e42b6ULL, 0x597f299cfc657e2aULL, 0x5fcb6fab3ad6faecULL, 0x6c44198c4a475817ULL
    };

    // Initialize hash values
    uint64_t h[8] = {
        0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL, 0x3c6ef372fe94f82bULL, 0xa54ff53a5f1d36f1ULL,
        0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL, 0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL
    };

    // Pre-processing: padding the message
    std::string paddedMessage = message;
    paddedMessage += '\x80'; // append single '1' bit
    while ((paddedMessage.size() % 128) != 112) // append '0' bits until message length ≡ 896 (mod 1024)
        paddedMessage += '\x00';
    uint64_t bitLength = message.size() * 8; // append length of message in bits as 64-bit big-endian integer
    for (int i = 7; i >= 0; --i)
        paddedMessage += static_cast<char>((bitLength >> (i * 8)) & 0xFF);

    // Process message in 1024-bit blocks
    for (size_t chunk = 0; chunk < paddedMessage.size(); chunk += 128) {
        // Initialize message schedule (W)
        uint64_t W[80];
        for (int t = 0; t < 16; ++t) {
            W[t] = 0;
            for (int j = 0; j < 8; ++j)
                W[t] |= static_cast<uint64_t>(static_cast<uint8_t>(paddedMessage[chunk + t * 8 + j])) << ((7 - j) * 8);
        }
        for (int t = 16; t < 80; ++t) {
            uint64_t s0 = ((W[t - 15] >> 1) | (W[t - 15] << 63)) ^ ((W[t - 15] >> 8) | (W[t - 15] << 56)) ^ (W[t - 15] >> 7);
            uint64_t s1 = ((W[t - 2] >> 19) | (W[t - 2] << 45)) ^ ((W[t - 2] >> 61) | (W[t - 2] << 3)) ^ (W[t - 2] >> 6);
            W[t] = W[t - 16] + s0 + W[t - 7] + s1;
        }

        // Initialize working variables
        uint64_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], h0 = h[7];

        // Compression function main loop
        for (int t = 0; t < 80; ++t) {
            uint64_t S1 = ((e >> 14) | (e << 50)) ^ ((e >> 18) | (e << 46)) ^ ((e >> 41) | (e << 23));
            uint64_t ch = (e & f) ^ (~e & g);
            uint64_t temp1 = h0 + S1 + ch + k[t] + W[t];
            uint64_t S0 = ((a >> 28) | (a << 36)) ^ ((a >> 34) | (a << 30)) ^ ((a >> 39) | (a << 25));
            uint64_t maj = (a & b) ^ (a & c) ^ (b & c);
            uint64_t temp2 = S0 + maj;

            h0 = g;
            g = f;
            f = e;
            e = d + temp1;
            d = c;
            c = b;
            b = a;
            a = temp1 + temp2;
        }

        // Add the compressed chunk to the current hash value
        h[0] += a;
        h[1] += b;
        h[2] += c;
        h[3] += d;
        h[4] += e;
        h[5] += f;
        h[6] += g;
        h[7] += h0;
    }

    // Produce the final hash value (big-endian)
    std::stringstream ss;
    for (int i = 0; i < 8; ++i)
        ss << std::hex << std::setw(16) << std::setfill('0') << h[i];
    return ss.str();
}

int main() {
    std::string message;
    std::cout << "Enter the message: ";
    std::getline(std::cin, message);

    // Calculate the SHA-512 hash of the message
    std::string hash = sha512(message);

    std::cout << "The SHA-512 hash of \"" << message << "\" is:\n" << hash << std::endl;

    // Extract the first 64 characters as the MAC
    std::string mac = hash.substr(0, 64);

    std::cout << "The MAC for the message is: " << mac << std::endl;

    return 0;
}
