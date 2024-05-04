void NuOsc::initialize()
{
/*
    Initialize the components of \rho and \bar{\rho}
    here.
*/
    /*double signu =  10000.0;  //0.6;
    double sigbnu = 10000.0; //0.5;*/
    double alpha = 0.7; //0.9;
    double epsn;
    double epsnq;
    /*double epsp;
    double epspq;*/
    std::ofstream g_file(ID + "_G0.bin",std::ofstream::out | std::ofstream::binary);
    if(!g_file)
    {
        std::cout << "Unable to open " << ID+"_G0.bin" << " file from NuOsc::initialise." 
        << "Will not be storing initial angular profiles.\n";
    }
    /*for (int i = 0; i < nvz; i++)
    {
        for (int j = 0; j < nz; j++)
        {*/
    /*
    epsn=eps(0.0, 0.0, perturbation_size);
    epsnq=sqrt(1.0 - epsn * epsn);
    epsp=epsn;
    epspq=epsnq;
    epsp=eps(0.0, 0.0, perturbation_size);
    epspq=sqrt(1.0 - epsn * epsn);
    */
    for (int j = 0; j < nz; j++)
    {
        /*
        //epsn=0;
        //epsnq=1;
        // epsn=eps(0.0, 0.0, perturbation_size);
        // epsnq=sqrt(1.0 - epsn * epsn); 
        */
        for (int i = 0; i < nvz; i++)
        {
            //epsn=eps(Z[j], 0.0, perturbation_size);
            //epsq=sqrt(1.0 - epsn * epsn);
            epsn=eps(0.0, 0.0, perturbation_size);
            epsnq=sqrt(1.0 - epsn * epsn);  
            //G0->G[idx(i, j)] = g(vz[i], 1.0, signu);
            //G0->bG[idx(i, j)] = alpha * g(vz[i], 1.0, sigbnu);
            G0->G[idx(i, j)] = 1.0 / (float)nvz; /* //g(vz[i], 1.0, signu);*/
            G0->bG[idx(i, j)] = alpha * 1.0 / (float)nvz; /* //* g(vz[i], 1.0, sigbnu);*/
            
            v_stat->ee[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0 + epsnq); 
            v_stat->xx[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0 - epsnq);
            //v_stat->ee[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0); 
            //v_stat->xx[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0);
            v_stat->ex_re[idx(i, j)] = 0.5 * G0->G[idx(i, j)] * (0.0 + epsn);
            v_stat->ex_im[idx(i, j)] = -0.0;
            //epsn=eps(Z[j], 0.0, perturbation_size);
            //epsnq=sqrt(1.0 - epsn * epsn);
            v_stat->bee[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0 + epsnq); 
            v_stat->bxx[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0 - epsnq);
            //v_stat->bee[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0); 
            //v_stat->bxx[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0);
            v_stat->bex_re[idx(i, j)] = 0.5 * G0->bG[idx(i, j)] * (0.0 + epsn);
            v_stat->bex_im[idx(i, j)] = 0.0;

            g_file.write((char *)&G0->G [idx(i, j)], sizeof(double)); 
            g_file.write((char *)&G0->bG[idx(i, j)], sizeof(double));

            
        }
    }
    /*
    for (int j = 0; j < nz; j++)
    {
        //epsn=0;
        //epsq=1;
        for (int i = 0; i < nvz/2; i++)
        {
            //epsn=eps(Z[j], 0.0, perturbation_size);
            //epsq=sqrt(1.0 - epsn * epsn);
            //G0->G[idx(i, j)] = g(vz[i], 1.0, signu);
            //G0->bG[idx(i, j)] = alpha * g(vz[i], 1.0, sigbnu);
            G0->G[idx(i, j)] = 1.0 / (float)nvz; //g(vz[i], 1.0, signu);
            G0->bG[idx(i, j)] = alpha * 1.0 / (float)nvz; //* g(vz[i], 1.0, sigbnu);
            
            v_stat->ee[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0 + epsnq); 
            v_stat->xx[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0 - epsnq);
            //v_stat->ee[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0); 
            //v_stat->xx[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0);
            v_stat->ex_re[idx(i, j)] = 0.5 * G0->G[idx(i, j)] * (0.0 + epsn);
            v_stat->ex_im[idx(i, j)] = -0.0;
            //epsn=eps(Z[j], 0.0, perturbation_size);
            //epsnq=sqrt(1.0 - epsn * epsn);
            v_stat->bee[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0 + epsnq); 
            v_stat->bxx[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0 - epsnq);
            //v_stat->bee[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0); 
            //v_stat->bxx[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0);
            v_stat->bex_re[idx(i, j)] = 0.5 * G0->bG[idx(i, j)] * (0.0 + epsn);
            v_stat->bex_im[idx(i, j)] = 0.0;

            g_file.write((char *)&G0->G [idx(i, j)], sizeof(double)); 
            g_file.write((char *)&G0->bG[idx(i, j)], sizeof(double));

            
        }
        for (int i = nvz/2; i < nvz; i++)
        {
            //epsn=eps(Z[j], 0.0, perturbation_size);
            //epsnq=sqrt(1.0 - epsn * epsn);
            //G0->G[idx(i, j)] = g(vz[i], 1.0, signu);
            //G0->bG[idx(i, j)] = alpha * g(vz[i], 1.0, sigbnu);
            G0->G[idx(i, j)] = 1.0 / (float)nvz; //g(vz[i], 1.0, signu);
            G0->bG[idx(i, j)] = alpha * 1.0 / (float)nvz; //* g(vz[i], 1.0, sigbnu);
            
            v_stat->ee[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0 + epspq); 
            v_stat->xx[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0 - epspq);
            //v_stat->ee[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0); 
            //v_stat->xx[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0);
            v_stat->ex_re[idx(i, j)] = 0.5 * G0->G[idx(i, j)] * (0.0 + epsp);
            v_stat->ex_im[idx(i, j)] = -0.0;
            //epsn=eps(Z[j], 0.0, perturbation_size);
            //epspq=sqrt(1.0 - epsn * epsn);
            v_stat->bee[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0 + epspq); 
            v_stat->bxx[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0 - epspq);
            //v_stat->bee[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0); 
            //v_stat->bxx[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0);
            v_stat->bex_re[idx(i, j)] = 0.5 * G0->bG[idx(i, j)] * (0.0 + epsp);
            v_stat->bex_im[idx(i, j)] = 0.0;

            g_file.write((char *)&G0->G [idx(i, j)], sizeof(double)); 
            g_file.write((char *)&G0->bG[idx(i, j)], sizeof(double));

            
        }
    }*/
    updateBufferZone(v_stat);
    std::cout << "Simulation state initialized." << std::endl;
    g_file.close();
}
