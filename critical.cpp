#pragma omp parallel
{
    VectorXd private_r = VectorXd::Zero(nt,1), private_i = VectorXd::Zero(nt,1);

#pragma omp for nowait 
    for (int ib = 0; ib < Nw; ++ib){
        VectorXd temp_r = VectorXd::Zero(nt,1), temp_i = VectorXd::Zero(nt,1);
        byfree(zztreal, zztimag, W[ib].zetaf_r, W[ib].zetaf_i, W[ib].gammaf, temp_r, temp_i);
        private_r += temp_r;
        private_i += temp_i;
    }

#pragma omp critical
    {
        pot_r += private_r;
        pot_i += private_i;
    }

}

