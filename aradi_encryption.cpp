/*=========aradi_block_cipher_encryption.cpp=========*/
#include <iostream>
#include <vector>
#include <tuple>
#include <cstdint>
#include <random>
bool last_linear_layer = true;

using namespace std;
string toBinaryString(unsigned int value) {
    string binary = "";
    for (int i = 31; i >= 0; --i) {
        binary += ((value >> i) & 1) ? '1' : '0';
    }
    return binary;
}

//==============Left circular shift=================
uint32_t left_circular_shift(uint32_t x, int shift) {
    return (x << shift) | (x >> (32 - shift));
}

uint16_t left_circular_shift16(uint16_t x, int shift) {
    return (x << shift) | (x >> (16 - shift));
}

//=========Function L=============================
uint32_t L(int a, int b, int c, uint32_t value) {
    uint16_t x = (uint16_t)(value >> 16);
    uint16_t y = (uint16_t)(value & 0xFFFF);

    uint32_t newX = ((x ^ left_circular_shift16(x, a) ^ left_circular_shift16(y, c)) << 16);
    uint32_t newY = (y ^ left_circular_shift16(y, a) ^ left_circular_shift16(x, b));

    return newX | newY;
}

//===========Function M=============================
pair<uint32_t, uint32_t> M(int i, int j, uint32_t X, uint32_t Y) {
    uint32_t S_i_Y = left_circular_shift(Y, i);
    uint32_t S_j_X = left_circular_shift(X, j);
    return make_pair(S_i_Y ^ S_j_X ^ X, S_i_Y ^ X);
}

//===========Key expansion==========================
void key_expansion(vector<vector<uint32_t>>& k, vector<uint32_t>& K, int nround) {
    for (int i = 0; i < nround; i++) {
        int j = i % 2;
        k[i][0] = K[4 * j + 0];
        k[i][1] = K[4 * j + 1];
        k[i][2] = K[4 * j + 2];
        k[i][3] = K[4 * j + 3];

        tie(K[1], K[0]) = M(1, 3, K[1], K[0]);
        tie(K[3], K[2]) = M(9, 28, K[3], K[2]);
        tie(K[5], K[4]) = M(1, 3, K[5], K[4]);
        tie(K[7], K[6]) = M(9, 28, K[7], K[6]);
        K[7] ^= i;
        if (j == 0) {
            swap(K[1], K[2]);

            swap(K[5], K[6]);
        }
        else {
            swap(K[1], K[4]);
            swap(K[3], K[6]);
        }
    }
    k[nround][0] = K[0];
    k[nround][1] = K[1];
    k[nround][2] = K[2];
    k[nround][3] = K[3];
}

//==========Encryption function=====================
void aradi_encrypt(uint32_t& w, uint32_t& x, uint32_t& y, uint32_t& z, vector<vector<uint32_t>>& k, int nround) {
    int a[4] = {11, 10, 9, 8};
    int b[4] = {8, 9, 4, 9};
    int c[4] = {14, 11, 14, 7};

    for (int i = 0; i < nround; i++) {
        z ^= k[i][3];
        y ^= k[i][2];
        x ^= k[i][1];
        w ^= k[i][0];

        x ^= (w & y);
        z ^= (x & y);
        y ^= (w & z);
        w ^= (x & z);

        int j = i % 4;
        z = L(a[j], b[j], c[j], z);
        y = L(a[j], b[j], c[j], y);
        x = L(a[j], b[j], c[j], x);
        w = L(a[j], b[j], c[j], w);

    }
    if(nround == 16){
        z ^= k[nround][3];
        y ^= k[nround][2];
        x ^= k[nround][1];
        w ^= k[nround][0];
    }
}

int main() {
    cout<<"======== Aradi paper's test vector =============="<<endl;
    uint32_t w = 0x0, x = 0x0, y = 0x0, z = 0x0;
    cout << "Plaintext: "<< hex << w << " " << x << " " << y << " " << z << endl;
    vector<uint32_t> K = {0x03020100, 0x07060504, 0x0b0a0908, 0x0f0e0d0c,
                          0x13121110, 0x17161514, 0x1b1a1918, 0x1f1e1d1c};
    cout << "Master key: "<< hex << K[0] << " " << K[1] << " " << K[2] << " " << K[3] << K[4] << " " << K[5] << " " << K[6] << " " << K[7] << endl;
    int nround=16;
    vector<vector<uint32_t>> k(nround + 1, vector<uint32_t>(4));
    key_expansion(k, K, nround);
    aradi_encrypt(w, x, y, z, k, nround);
    cout << "Ciphertext: "<< hex << w << " " << x << " " << y << " " << z << endl;

    cout<<"\n======== Distingusher for 5 rounds Aradi =============="<<endl;
    // Initialize random number generator
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint32_t> dis(0, 0xFFFFFFFF);
    nround=5;
    int nb_random_keys = 16;
    for (int j = 0; j < nb_random_keys; ++j) { // j iterates over number of master keys
        // generate a random key K
        vector<uint32_t> K(8);
        for (int i = 0; i < 8; ++i) {
            K[i] = dis(gen);
        }

        // key expansion
        vector<vector<uint32_t>> k(nround + 1, vector<uint32_t>(4));
        key_expansion(k, K, nround);

        uint32_t final_w = 0, final_x = 0, final_y = 0, final_z = 0;
        for (uint32_t i = 0; i < (1 << 13); ++i){
            uint32_t w = i<<7;
            uint32_t x = 0x0;
            uint32_t y = 0x0;
            uint32_t z = 0x0;

            aradi_encrypt(w, x, y, z, k, nround);

            final_w ^= w;
            final_x ^= x;
            final_y ^= y;
            final_z ^= z;
        }
        cout<<"=========================================="<<endl;
        cout << "key "<<dec<<j<<": " << hex << final_w << " " << final_x << " " <<final_y << " " << final_z << endl;
    }
    cout<<endl;
    cout<<"For 5 rounds, we can observe that for " << dec << nb_random_keys <<
    " random keys, the consecutive bytes are equal in the words x and z, for the selected cube-index sets." <<endl;
    cout<<endl;

    return 0;
}