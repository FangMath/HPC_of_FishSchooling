/* Define the class Wing 
 * Copyright Fang Fang, 04/27/2015, Courant Inst.
 */
#ifndef FISH_WING_H
#define FISH_WING_H

#include <Eigen/Dense>
using namespace Eigen;

namespace Fish {
    class Wing{
        public:
        int K;
        VectorXd zetaf_r, zetaf_i;
        VectorXd gammaf;

        VectorXd zetabody_r, zetabody_i;
        VectorXd nu;

        VectorXd source_r, source_i;
        VectorXd distr, delta;

        VectorXd pot_r, pot_i;

        /*
        int pK;
        VectorXd pointVX_r, pointVX_i;
        VectorXd CirlVX_r, Cirl_VX_i;

        VectorXd nu;
        VectorXd zetabody_r, zetabody_i;
        */

        Wing (){
            K = 0;
        }

    };
/*
        void set_values (int,int);
        int area() {return width*height;}
    void Wing::set_values (int x, int y) {
        width = x;
        height = y;
    }
    */

}// end of namespace Fish
#endif // FISH_WING_H 
