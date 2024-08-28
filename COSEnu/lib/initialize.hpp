void NuOsc::initialize()
{
/*
    Initialize the components of \rho and \bar{\rho}
    here.
*/
    // double signu =  10000.0;  //0.6;
    // double sigbnu = 10000.0; //0.5;
    double alpha = 0.7; //0.9;
    double epsn, epsnr, epsni, bepsnr, bepsni ,bepsnq;
    double epsnq;
    double epsp;
    double epspq;
    double rdph[20000];
    std::ofstream g_file(ID + "_G0.bin",std::ofstream::out | std::ofstream::binary);
    if(!g_file)
    {
        std::cout << "Unable to open " << ID+"_G0.bin" << " file from NuOsc::initialise." 
        << "Will not be storing initial angular profiles.\n";
    }



    //  JORDAN ADD JORDAN ADD ///////////////////////////////////////////////////////////////////////////
    // std::printf("從文件中讀取的浮點數組 12345：\n");
    // std::ifstream inputFile("data1.txt");
    // // if (!inputFile) {
    // //     std::cerr << "Unable to open, JW" << std::endl;
    // //     // return 1;
    // // }
    // std::vector<float> numbers;
    // float number;
    // while (inputFile >> number) {
    //     numbers.push_back(number);
    // }

    // inputFile.close();

    // std::printf("從文件中讀取的浮點數組：\n");
    // for (size_t i = 0; i < numbers.size(); ++i) {
    //     std::printf("%.5f ", numbers[i]);
    // }
    // std::printf("\n");

    /////////////////////////////////////////////////////////////////////////////////////////////





    /*for (int i = 0; i < nvz; i++)
    {874
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
    // FILE* f;
    /*---------------------------------------------*/
    unsigned int seed = 77777;
    srand(seed);
    for(int i = 0; i<20000;i++){
        rdph[i] = (double)rand() / RAND_MAX * 2 * M_PI;
    }
    /*--------------------------------------------------*/
    for (int j = 0; j < nz; j++)
    {
        epsn=0.0;
        //epsnq=1;
        // epsn=eps(0.0, 0.0, perturbation_size);
        // epsnq=sqrt(1.0 - epsn * epsn); 
        for (int i = 0; i < nvz; i++)
        {
            // rdph = (double)rand() / RAND_MAX * 2 * M_PI;
            // epsn=eps(Z[j], 0.0, perturbation_size);
            epsnr = epsr(Z[j], 0.0, perturbation_size, rdph);
            epsni = epsi(Z[j], 0.0, perturbation_size, rdph);
            // epsq=sqrt(1.0 - epsn * epsn);
            // epsn=eps(0.0, 0.0, perturbation_size);
            // epsn=eps(Z[j], 0.0, perturbation_size);
            
            // epsn=0.6;
            //////////////////////////////////////////////////////////////////////////////////////////////
            // G0->G[idx(i, j)] = 1.0/2.0; //1.0 / (float)nvz; g(vz[i], 1.0, signu); f(Z[i])
            // G0->bG[idx(i, j)] = alpha * 1.0/2.0; 
            G0->G[idx(i, j)] = 1.0 / (float)nvz;// g(vz[i], 1.0, signu); f(Z[i])
            G0->bG[idx(i, j)] = alpha * 1.0/(float)nvz;
            epsnq=sqrt(1.0 - epsnr * epsnr - epsni * epsni);

            /*--------------------------------------------------------------------------------*/  
            // G0->G[idx(i, j)] = g(vz[i], 1.0, signu);
            //G0->bG[idx(i, j)] = alpha * g(vz[i], 1.0, sigbnu);
            // G0->G[idx(i, j)] = 1.0 / (float)nvz; // g(vz[i], 1.0, signu);
            // G0->bG[idx(i, j)] = alpha * 1.0 / (float)nvz;
            // if(j > 5000 && i==0){
            //     printf("%.25f\n", eig(j));

            // }
            // if ((i == 0) ||( i == nvz  -1))
            // {
            //     G0->G[idx(i, j)] = nvz/4.0; //1.0 / (float)nvz; // g(vz[i], 1.0, signu);
            //     G0->bG[idx(i, j)] = alpha * nvz /4.0; //1.0 / (float)nvz; //* g(vz[i], 1.0, sigbnu);
            // }
            // else
            // {
            //     G0->G[idx(i, j)] = 0.0;
            //     G0->bG[idx(i, j)] = 0.0;
            // }
            // printf("TEST %f\n", numbers[j]);
            // printf("TEST\n");
            /*-------------------------------------------------------------------------------------*/
            v_stat->ee[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0 + epsnq); 
            v_stat->xx[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0 - epsnq);
            //v_stat->ee[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0); 
            //v_stat->xx[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0);
            // v_stat->ex_re[idx(i, j)] = 0.5 * G0->G[idx(i, j)] * (0.0 + epsn);
            v_stat->ex_re[idx(i, j)] = 0.5 * G0->G[idx(i, j)] * (0.0 + epsnr);
            // v_stat->ex_im[idx(i, j)] = -0.0;
            v_stat->ex_im[idx(i, j)] = 0.5 * G0->G[idx(i, j)] * (0.0 + epsni);
            //epsn=eps(Z[j], 0.0, perturbation_size);
            //epsnq=sqrt(1.0 - epsn * epsn);
            v_stat->bee[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0 + epsnq); 
            v_stat->bxx[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0 - epsnq);
            //v_stat->bee[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0); 
            //v_stat->bxx[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0);
            // v_stat->bex_re[idx(i, j)] = 0.5 * G0->bG[idx(i, j)] * (0.0 + epsn);
            v_stat->bex_re[idx(i, j)] = 0.5 * G0->bG[idx(i, j)] * (0.0 + epsnr);
            v_stat->bex_im[idx(i, j)] = 0.5 * G0->bG[idx(i, j)] * (0.0 + epsni);
            // v_stat->bex_im[idx(i, j)] = 0.0;

            //////////// for manual input initialize file///////////////////////////////////////////////////////////
            // if(j < 20 && i==0){
            //     printf("%.20f\n", eig(j));
            //     epsn = perturbation_size * eig(j);
            //     printf("%.20f\n",epsn);

            // }
            // epsnr = perturbation_size * eigr(j,i);
            // epsni = perturbation_size * eigr(j+nz,i);
            // bepsnr = perturbation_size * eigr(j+2 * nz,i);
            // bepsni = perturbation_size * eigr(j+ 3 * nz,i);
            // // printf("%f\n",epsn);
            // epsnq=sqrt(1.0 - epsnr * epsnr - epsni * epsni);
            // bepsnq=sqrt(1.0 - bepsnr * bepsnr - bepsni * bepsni);  
            // v_stat->ee[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0 + epsnq); 
            // v_stat->xx[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0 - epsnq);
            // //v_stat->ee[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0); 
            // //v_stat->xx[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0);
            // v_stat->ex_re[idx(i, j)] = 0.5 * G0->G[idx(i, j)] * (0.0 + epsnr);
            // v_stat->ex_im[idx(i, j)] = 0.5 * G0->G[idx(i, j)] * (0.0 + epsni);
            // //epsn=eps(Z[j], 0.0, perturbation_size);
            // //epsnq=sqrt(1.0 - epsn * epsn);
            // v_stat->bee[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0 + bepsnq); 
            // v_stat->bxx[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0 - bepsnq);
            // //v_stat->bee[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0); 
            // //v_stat->bxx[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0);
            // v_stat->bex_re[idx(i, j)] = 0.5 * G0->bG[idx(i, j)] * (0.0 + bepsnr);
            // v_stat->bex_im[idx(i, j)] = 0.5 * G0->bG[idx(i, j)] * (0.0 + bepsni);
            

            // g_file.write((char *)&G0->G [idx(i, j)], sizeof(double)); 
            // g_file.write((char *)&G0->bG[idx(i, j)], sizeof(double));
            //////////end of manual innitialize ////////////////////////////////////////////////////////////////
            // printf("rho00 %f\n",0.5 * G0->G[idx(i, j)] * (1.0 + epsnq));
            
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
