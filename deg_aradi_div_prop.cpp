#include "gurobi_c++.h"
#include <iostream>
#include <vector>
#include <tuple>
#include <cstdint>
#include <iomanip>
#include <numeric>
#include<bitset>
#include<algorithm>
#include<map>
#include<iomanip>
#include<cmath>
bool last_linear_layer = true;

using namespace std;
vector<GRBVar> circular_shift(const vector<GRBVar> &, int);
vector<GRBVar> add_xor_constraints(GRBModel &, const vector<GRBVar> &, const vector<GRBVar> &);
void add_copy_constraints(GRBModel &, vector<GRBVar> &, vector<vector<GRBVar>> &);
void apply_constraints_and_optimize(GRBModel &, vector<vector<GRBVar>> &, vector<vector<GRBVar>> &, vector<vector<GRBVar>> &, vector<vector<GRBVar>> &, int);

vector<GRBVar> add_and_constraints(GRBModel &model, const vector<GRBVar> &a, const vector<GRBVar> &b) {
    int size = a.size();
    vector<GRBVar> a_and_b(size);
    for (int i = 0; i < size; ++i) {
        a_and_b[i] = model.addVar(0, 1, 0, GRB_BINARY);
        model.addConstr(a_and_b[i] == a[i] );
        model.addConstr(a_and_b[i] == b[i]);
    }
    return a_and_b;
}

tuple<vector<GRBVar>, vector<GRBVar>, vector<GRBVar>, vector<GRBVar>>
aradi_nonlinear_layer(GRBModel& model, vector<GRBVar>& w, vector<GRBVar>& x, vector<GRBVar>& y, vector<GRBVar>& z) {

    vector<vector<GRBVar>> w_copy(3, vector<GRBVar>(32));
    add_copy_constraints(model, w, w_copy);
    vector<vector<GRBVar>> y_copy(3, vector<GRBVar>(32));
    add_copy_constraints(model, y, y_copy);

    // modelling of the first equation: x = x + w * y
    auto w_and_y = add_and_constraints(model, w_copy[0], y_copy[0]);
    auto x_new  = add_xor_constraints(model, x, w_and_y);

    vector<vector<GRBVar>> x_new_copy(3, vector<GRBVar>(32));
    add_copy_constraints(model, x_new, x_new_copy);

    //modelling of the second equation: z  = z + x * y

    auto x_new_and_y = add_and_constraints(model, x_new_copy[0], y_copy[1]);
    auto z_new  = add_xor_constraints(model, z, x_new_and_y);

    vector<vector<GRBVar>> z_new_copy(3, vector<GRBVar>(32));
    add_copy_constraints(model, z_new, z_new_copy);

    // modeling of third equation: y = y + z * w

    auto z_new_and_w = add_and_constraints(model, z_new_copy[0], w_copy[1]);
    auto y_new  = add_xor_constraints(model, y_copy[2], z_new_and_w);

    //modelling of fourth equation: w = w + x *z

    auto x_new_and_z_new = add_and_constraints(model, x_new_copy[1], z_new_copy[1]);
    auto w_new  = add_xor_constraints(model, w_copy[2], x_new_and_z_new);

    return make_tuple(w_new, x_new_copy[2], y_new, z_new_copy[2]);
}

vector<GRBVar> L(GRBModel &model, int a, int b, int c, vector<GRBVar> &value) {


    string varName = value[0].get(GRB_StringAttr_VarName);
    string v=string(1, varName[0]);

    vector<vector<GRBVar>> value_copy(3, vector<GRBVar>(32));
    add_copy_constraints(model, value, value_copy);

    vector<GRBVar> x(value_copy[0].begin(), value_copy[0].begin() + 16);
    vector<GRBVar> y(value_copy[0].begin() + 16, value_copy[0].end());

    vector<GRBVar> x_a(value_copy[1].begin(), value_copy[1].begin() + 16);
    vector<GRBVar> y_a(value_copy[1].begin() + 16, value_copy[1].end());

    vector<GRBVar> x_b(value_copy[2].begin(), value_copy[2].begin() + 16);
    vector<GRBVar> y_c(value_copy[2].begin() + 16, value_copy[2].end());
    model.update();

    auto shiftX_a = circular_shift(x_a, a);
    auto shiftY_c = circular_shift(y_c, c);
    model.update();

    auto Lx = add_xor_constraints(model, x, shiftX_a);
    auto Lxy = add_xor_constraints(model, Lx, shiftY_c);

    auto shiftY_a = circular_shift(y_a, a);
    auto shiftX_b = circular_shift(x_b, b);
    auto Ly = add_xor_constraints(model, y, shiftY_a);
    auto Lyx = add_xor_constraints(model, Ly, shiftX_b);

    vector<GRBVar> newX(Lxy.begin(), Lxy.end());
    vector<GRBVar> newY(Lyx.begin(), Lyx.end());

    newX.insert(newX.end(), newY.begin(), newY.end());
     model.update();
    return newX;
}

void add_copy_constraints(GRBModel &model, vector<GRBVar> &orig, vector<vector<GRBVar>> &copies) {
    int num_bits = orig.size();
    int num_copies = copies.size();
    for (int i = 0; i < num_bits; i++) {
        GRBLinExpr sum_of_copies = 0;
        for (int c = 0; c < num_copies; c++) {
            copies[c][i]=model.addVar(0, 1, 0, GRB_BINARY);
            model.addConstr(orig[i] >= copies[c][i]);
            sum_of_copies += copies[c][i];  // sum of copies
        }
        model.addConstr(sum_of_copies >= orig[i]);
    }
}

vector<GRBVar> add_xor_constraints(GRBModel &model, const vector<GRBVar> &a, const vector<GRBVar> &b) {
    int size = a.size();
    vector<GRBVar> a_xor_b(size);
    for (int i = 0; i < size; ++i) {
        a_xor_b[i] = model.addVar(0, 1, 0, GRB_BINARY);
        model.addConstr(a_xor_b[i] == a[i] + b[i]);
    }
    return a_xor_b;
}

void add_xor_constant_constraints(GRBModel &model, vector<GRBVar> &K7, const int i){

    for(int it=0; it<4; it++){
        if((i>>it) & 1){
         GRBVar t;
         t=model.addVar(0, 1, 0, GRB_BINARY);
         model.addConstr(t >= K7[31-it]);
         K7[31-it] = t;
        }
    }

}

vector<GRBVar> circular_shift(const vector<GRBVar> &input, int shift) {

    int size = input.size();
    vector<GRBVar> output(size);
    for (int i = 0; i < size; ++i) {
        output[i] = input[(i + shift) % size];
    }
    return output;
}

pair<vector<GRBVar>, vector<GRBVar>> M(GRBModel &model, int i, int j,  vector<GRBVar> &X_copy1,  vector<GRBVar> &X_copy2, vector<GRBVar> &Y) {
    auto S_i_Y = circular_shift(Y, i);
    model.update();
    auto second = add_xor_constraints(model, S_i_Y, X_copy1);

    auto S_j_X = circular_shift(X_copy2, j);
    model.update();
    vector<vector<GRBVar>> s_copy(2, vector<GRBVar>(32));
    add_copy_constraints(model, second, s_copy);
     model.update();
    auto first = add_xor_constraints(model, s_copy[0], S_j_X);
    model.update();
    return make_pair(first, s_copy[1]);
}

void key_expansion(GRBModel &model, vector<vector<vector<GRBVar>>> &k, vector<vector<GRBVar>>& K, int rnd) {
    if(rnd == 1){
        k[0][0] = K[0];
        k[0][1] = K[1];
        k[0][2] = K[2];
        k[0][3] = K[3];

        for(int i=0; i<32; i++){
           model.addConstr(K[4][i] == 0);
           model.addConstr(K[5][i] == 0);
           model.addConstr(K[6][i] == 0);
           model.addConstr(K[7][i] == 0);
        }
    }
    else{
        for(int i = 0; i < rnd-1; i++) {
            int j = i % 2;
            if (j == 0){
                vector<vector<GRBVar>> K0_copy(2, vector<GRBVar>(32));
                add_copy_constraints(model, K[0], K0_copy);

                vector<vector<GRBVar>> K1_copy(3, vector<GRBVar>(32));
                add_copy_constraints(model, K[1], K1_copy);

                vector<vector<GRBVar>> K2_copy(2, vector<GRBVar>(32));
                add_copy_constraints(model, K[2], K2_copy);

                vector<vector<GRBVar>> K3_copy(3, vector<GRBVar>(32));
                add_copy_constraints(model, K[3], K3_copy);
                k[i][0] = K0_copy[0];
                k[i][1] = K1_copy[0];
                k[i][2] = K2_copy[0];
                k[i][3] = K3_copy[0];

                model.update();
                pair<vector<GRBVar>, vector<GRBVar>> result1 = M(model, 1, 3, K1_copy[1], K1_copy[2], K0_copy[1]);
                K[1] = result1.first;
                K[0] = result1.second;


                tie(K[3], K[2]) = M(model, 9, 28, K3_copy[1], K3_copy[2], K2_copy[1]);
                vector<vector<GRBVar>> K5_copy(2, vector<GRBVar>(32));
                add_copy_constraints(model, K[5], K5_copy);

                tie(K[5], K[4]) = M(model, 1, 3, K5_copy[0],K5_copy[1], K[4]);


                vector<vector<GRBVar>> K7_copy(2, vector<GRBVar>(32));
                add_copy_constraints(model, K[7], K7_copy);

                tie(K[7], K[6]) = M(model, 9, 28, K7_copy[0], K7_copy[1], K[6]);

                model.update();
                add_xor_constant_constraints(model, K[7], i);

                swap(K[1], K[2]);
                swap(K[5], K[6]);
                model.update();
            }
            else{
                vector<vector<GRBVar>> K4_copy(2, vector<GRBVar>(32));
                add_copy_constraints(model, K[4], K4_copy);

                vector<vector<GRBVar>> K5_copy(3, vector<GRBVar>(32));
                add_copy_constraints(model, K[5], K5_copy);

                vector<vector<GRBVar>> K6_copy(2, vector<GRBVar>(32));
                add_copy_constraints(model, K[6], K6_copy);

                vector<vector<GRBVar>> K7_copy(3, vector<GRBVar>(32));
                add_copy_constraints(model, K[7], K7_copy);

                k[i][0] = K4_copy[0];
                k[i][1] = K5_copy[0];
                k[i][2] = K6_copy[0];
                k[i][3] = K7_copy[0];
                model.update();

                vector<vector<GRBVar>> K1_copy(2, vector<GRBVar>(32));
                add_copy_constraints(model, K[1], K1_copy);

                tie(K[1], K[0]) = M(model, 1, 3, K1_copy[0], K1_copy[1], K[0]);

                vector<vector<GRBVar>> K3_copy(2, vector<GRBVar>(32));
                add_copy_constraints(model, K[3], K3_copy);

                tie(K[3], K[2]) = M(model, 9, 28, K3_copy[0],K3_copy[1], K[2]);
                tie(K[5], K[4]) = M(model, 1, 3, K5_copy[1], K5_copy[2], K4_copy[1]);
                tie(K[7], K[6]) = M(model, 9, 28, K7_copy[1], K7_copy[2], K6_copy[1]);

                add_xor_constant_constraints(model, K[7], i);

                swap(K[1], K[4]);
                swap(K[3], K[6]);
                model.update();
           }
        }
    }
    if((rnd > 1) && (rnd % 2 == 0)){
        k[rnd-1][0] = K[4];
        k[rnd-1][1] = K[5];
        k[rnd-1][2] = K[6];
        k[rnd-1][3] = K[7];

        for(int i=0; i<32; i++){
            model.addConstr(K[0][i] == 0);
            model.addConstr(K[1][i] == 0);
            model.addConstr(K[2][i] == 0);
            model.addConstr(K[3][i] == 0);
        }
    }
    if ((rnd > 1) && (rnd % 2 == 1)){
        k[rnd-1][0] = K[0];
        k[rnd-1][1] = K[1];
        k[rnd-1][2] = K[2];
        k[rnd-1][3] = K[3];

        for (int i = 0; i < 32; i++){
            model.addConstr(K[4][i] == 0);
            model.addConstr(K[5][i] == 0);
            model.addConstr(K[6][i] == 0);
            model.addConstr(K[7][i] == 0);
        }
    }
}

void model_cipher(GRBModel &model, int rnd) {

   const int num_bits = 32;

    vector<vector<GRBVar>> w(rnd+1 , vector<GRBVar>(32));
    vector<vector<GRBVar>> x(rnd+1 , vector<GRBVar>(32));
    vector<vector<GRBVar>> y(rnd+1 , vector<GRBVar>(32));
    vector<vector<GRBVar>> z(rnd+1 , vector<GRBVar>(32));

    for (int r = 0; r < 1; r++) { //only first time vars need to be added
        for (int i = 0; i < num_bits; i++)
            w[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "w_" + to_string(r) + "_" + to_string(i));
    }
    for (int r = 0; r < 1; r++) {
        for (int i = 0; i < num_bits; i++)
            x[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "x_" + to_string(r) + "_" + to_string(i));
    }
    for (int r = 0; r < 1; r++) {
        for (int i = 0; i < num_bits; i++)
            y[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "y_" + to_string(r) + "_" + to_string(i));
    }
    for (int r = 0; r < 1; r++) {
        for (int i = 0; i < num_bits; i++)
            z[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "z_" + to_string(r) + "_" + to_string(i));
    }
    model.update();
    //   master key
    vector<vector<GRBVar>> K(8, vector<GRBVar>(32));

    for (int i = 0; i < 8; i++) {
        for(int j = 0; j < 32; j++)
        K[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "K_" + to_string(i) + "_" + to_string(j));
    }

    model.update();
    // round keys
    vector<vector<vector<GRBVar>>> k(rnd, vector<vector<GRBVar>>(4, vector<GRBVar>(32))); //round keys
    key_expansion(model, k, K, rnd);
    model.update();

     int a[4] = {11, 10, 9, 8};
     int b[4] = {8, 9, 4, 9};
     int c[4] = {14, 11, 14, 7};
     int j = 0;

    for (int round_it = 0; round_it < rnd; round_it++){

        auto w_tmp= add_xor_constraints(model, w[round_it], k[round_it][0]);
        auto x_tmp= add_xor_constraints(model, x[round_it], k[round_it][1]);
        auto y_tmp= add_xor_constraints(model, y[round_it], k[round_it][2]);
        auto z_tmp= add_xor_constraints(model, z[round_it], k[round_it][3]);

        auto result = aradi_nonlinear_layer(model, w_tmp, x_tmp, y_tmp, z_tmp);

        if ((round_it == rnd -1) && !last_linear_layer){
            w[round_it+1] = get<0>(result);
            x[round_it+1] = get<1>(result);
            y[round_it+1] = get<2>(result);
            z[round_it+1] = get<3>(result);
            model.update();
        }
        else{
            auto w_new = get<0>(result);
            auto x_new = get<1>(result);
            auto y_new = get<2>(result);
            auto z_new = get<3>(result);
            model.update();

            j = round_it % 4;
            z[round_it+1] = L(model, a[j], b[j], c[j], z_new);
            y[round_it+1] = L(model, a[j], b[j], c[j], y_new);
            x[round_it+1] = L(model, a[j], b[j], c[j], x_new);
            w[round_it+1] = L(model, a[j], b[j], c[j], w_new);
            model.update();
        }
    }
    /* adding final round key */
   /*
    auto wo= add_xor_constraints(model, w[rnd], k[rnd][0]);
    auto xo= add_xor_constraints(model, x[rnd], k[rnd][1]);
    auto yo= add_xor_constraints(model, y[rnd], k[rnd][2]);
    auto zo= add_xor_constraints(model, z[rnd], k[rnd][3]);
    model.update();
    */

   apply_constraints_and_optimize(model, w, x, y, z, rnd);

}

void apply_constraints_and_optimize(GRBModel &model,
                                    vector<vector<GRBVar>> &w,
                                    vector<vector<GRBVar>> &x,
                                    vector<vector<GRBVar>> &y,
                                    vector<vector<GRBVar>> &z,
                                    int rnd) {

    GRBLinExpr ll = 0;
    for(int i = 0 ; i < 32; i++){
	    ll += w[rnd][i];
	    ll += x[rnd][i];
	    ll += y[rnd][i];
	    ll += z[rnd][i];
    }
    model.addConstr(ll == 1) ;

    // setting constraints corresponding to cube-index sets for 4 rounds
    for(int i = 0 ; i < 32; i++){
        if(rnd == 4){
            if( i < 16 ){
                model.addConstr(w[0][i] >= 0); // 16 free vars
                model.addConstr(y[0][i] == 0);
                model.addConstr(x[0][i] == 0);
                model.addConstr(z[0][i] == 0);
            }
            else{
                model.addConstr(w[0][i] == 0);
                model.addConstr(x[0][i] == 0);
                model.addConstr(y[0][i] == 0);
                model.addConstr(z[0][i] == 0);
            }
        }
    }
    GRBLinExpr objective = 0;
    for(int i = 0 ; i<32; i++){
	    objective +=w[0][i];
	    objective +=x[0][i];
	    objective +=y[0][i];
	    objective +=z[0][i];
    }
    model.update();

    cout<<"model specification for "<< rnd<< " rounds "<<endl;
    cout<<"number of variables: "<< model.get(GRB_IntAttr_NumVars)<<endl;
    cout<<"number of constraint: "<<model.get(GRB_IntAttr_NumConstrs)<<endl;

    int deg;
    cout<<"========================================"<<endl;
    for (int target=0; target <32; target++){

        GRBModel model_copy = GRBModel(model);
        model_copy.update();
        model_copy.setObjective(objective, GRB_MAXIMIZE);
        model_copy.addConstr(w[rnd][target]  ==1);
        model_copy.optimize();
        deg = round(model_copy.getObjective().getValue()) ;
        cout << "Target bit: " <<"w["<< target <<"], deg: " << deg << endl;
    }
    cout<<"========================================"<<endl;
    for (int target=0; target <32; target++){
        GRBModel model_copy = GRBModel(model);
        model_copy.update();
        model_copy.setObjective(objective, GRB_MAXIMIZE);
        model_copy.addConstr(x[rnd][target]  ==1);
        model_copy.optimize();
        deg = round(model_copy.getObjective().getValue()) ;
        cout << "Target bit: " <<"x["<< target <<"], deg: " << deg << endl;
    }
    cout<<"========================================"<<endl;
    for (int target=0; target <32; target++){
        GRBModel model_copy = GRBModel(model);
        model_copy.update();
        model_copy.setObjective(objective, GRB_MAXIMIZE);
        model_copy.addConstr(y[rnd][target]  ==1);
        model_copy.optimize();
        deg = round(model_copy.getObjective().getValue()) ;
        cout << "Target bit: " <<"y["<< target <<"], deg: " << deg << endl;
    }
    cout<<"========================================"<<endl;
    for (int target=0; target <32; target++){
        GRBModel model_copy = GRBModel(model);
        model_copy.update();
        model_copy.setObjective(objective, GRB_MAXIMIZE);
        model_copy.addConstr(z[rnd][target]  ==1);
        model_copy.optimize();
        deg = round(model_copy.getObjective().getValue()) ;
        cout << "Target bit: " <<"z["<< target <<"], deg: " << deg << endl;
    }
}


int main() {

    try {
        GRBEnv env = GRBEnv(true);
        env.start();
        //env.set(GRB_IntParam_Threads, 32);
        env.set(GRB_IntParam_LogToConsole, 0);
        env.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_BESTBOUND);
        GRBModel model = GRBModel(env);

        int rounds = 4;
        model_cipher(model, rounds);
    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode( ) << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during optimization" << endl;
    }
    return 0;
}

