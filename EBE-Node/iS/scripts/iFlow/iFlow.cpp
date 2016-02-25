//===============================================================================
//  iIntegrateSpectraAndFlow
//===============================================================================
//
//  Chun Shen & Zhi Qiu
//    shen.201@asc.ohio-state.edu
//    qiu.24@asc.ohio-state.edu
//
//        Date: 03/2012
// Update version info in the main function.

// dy/deta=pT/mT; pT*sh(eta)=mT*sh(y); p=pT*ch(eta); E=mT*ch(y)

#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include<vector>

#include "Table.h"
#include "arsenal.h"

#define GAUSS_N 30 // used in the integrated vn calculation with pT cut

using namespace std;


Table pT_tab("tables/pT_gauss_table.dat");
int pT_tab_length = pT_tab.getNumberOfRows();
Table phi_tab("tables/phi_gauss_table.dat");
int phi_tab_length = phi_tab.getNumberOfRows();
Table eta_table("tables/particle_eta_uni_0.5_table.dat");
int eta_table_length = eta_table.getNumberOfRows();


void calculate_dN_deta_from_dN_dy_boostInv(Table* dN_deta, Table& dN_dy, double eta, double mass)
// Assuming boost invariance, calculate dN/(d(eta) pT dpT dphi) from dN/(dy pT dpT dphi)(y)=dN/(dy pT dpT dphi)(0).
// Both dN/d(eta) and dN/dy are arrays of pT and phi.
{
    for (int i=1; i<=pT_tab_length; i++) // pT loop
    {
        double pT = pT_tab.get(1,i)+1e-100; // division-friendly
        double mT = sqrt(mass*mass + pT*pT);
        double sinh_eta = sinh(eta);
        double dy_deta = cosh(eta)/sqrt(mT*mT/pT/pT + sinh_eta*sinh_eta);

        for (int j=1; j<=phi_tab_length; j++) // phi loop
        {
            dN_deta->set(i, j, dN_dy.get(i,j)*dy_deta);

            // if (i==pT_tab_length)
            // {
            //     cout << dN_dy.get(i,j) << " vs " << dN_deta->get(i,j) << endl;
            // }

        }

    }
}





//-------------------------------------------------------------------------
void calculate_flows(Table& dN_ptdptdphi, ostream& os_diff, ostream& os_inte, double mass, int to_order=9, double pT_cut_min=-1, double pT_cut_max=-1, int interp_mode=2)
// Calculate flow from order from_order to to_order and store them to files.
// If pT_cut>0, then the actual pT integral is performed from pT_cut to the
// maximum value in pT_tab.
{
    int from_order = 1;

    int number_of_flows = to_order-from_order+1;

    Table vn_diff(3+number_of_flows*3, pT_tab_length); // line format: pT, mT, dN/(pT dpT), flow_1_real, flow_1_imag, flow_1_norm, ...
    Table vn_inte(6, to_order+1); // line format: order# (starting from 0), numerator_real, numerator_imag, flow_real, flow_imag, flow_norm

    //---------------------
    // differential flow
    //---------------------
    double normalization[pT_tab_length]; // normalization factor
    for (int i=0; i<pT_tab_length; i++) normalization[i] = 0.0;

    double vn[pT_tab_length][number_of_flows][2]; // diff_flow numerators; 2: 0,1.real,imag
    for (int i=0; i<pT_tab_length; i++)
    for (int t=0; t<number_of_flows; t++)
    {vn[i][t][0]=0; vn[i][t][1]=0;}

    for (int i=0; i<pT_tab_length; i++)
    //for (int i=0; i<1; i++) // for debug
    {
    double pT = pT_tab.get(1,i+1); // pT_weight = pT_tab.get(2,i+1); // unused
    double mT = sqrt(mass*mass + pT*pT);

    // phi integration
    for(int j=0; j<phi_tab_length; j++)
    //for(int j=0; j<1; j++) // for debug
    {
      double phi = phi_tab.get(1,j+1), phi_weight = phi_tab.get(2,j+1);
      double dN = dN_ptdptdphi.get(i+1,j+1);

      normalization[i] += dN*phi_weight;
      for (int order=from_order; order<=to_order; order++)
      {
        vn[i][order-from_order][0] += dN*phi_weight*cos(order*phi);
        vn[i][order-from_order][1] += dN*phi_weight*sin(order*phi);
      }
    }

    normalization[i] = normalization[i] + 1e-30;
    // store values
    vn_diff.set(1, i+1, pT);
    vn_diff.set(2, i+1, mT-mass);
    vn_diff.set(3, i+1, normalization[i]/(2.0*M_PI)); // 2*pi: azimuthal angle averaged
    for (int t=0; t<number_of_flows; t++)
    {
      vn_diff.set(4+t*3, i+1, vn[i][t][0]/normalization[i]);
      vn_diff.set(5+t*3, i+1, vn[i][t][1]/normalization[i]);
      vn_diff.set(6+t*3, i+1, sqrt(vn[i][t][0]*vn[i][t][0]+vn[i][t][1]*vn[i][t][1])/normalization[i]);
    }

    }

    //---------------------
    // integrated flow
    //---------------------
    // initialization
    double normalizationi = 0;
    double vni[number_of_flows][2]; // integrated_flow numerators; 2: 0,1.real,imag
    for (int t=0; t<number_of_flows; t++) {vni[t][0]=0; vni[t][1]=0;}

    if (pT_cut_min < 0) // Calculate flow without pT_cut
    {

        for (int i=0; i<pT_tab_length; i++)
        {
            double pT = pT_tab.get(1,i+1), pT_weight = pT_tab.get(2,i+1);

            normalizationi += normalization[i]*pT*pT_weight;

            for (int order=from_order; order<=to_order; order++)
            {
              vni[order-from_order][0] += vn[i][order-from_order][0]*pT*pT_weight;
              vni[order-from_order][1] += vn[i][order-from_order][1]*pT*pT_weight;
            }
        }

    }
    else // Use a pT_cut
    {
        double max_pT = pT_tab.getLast(1);
        if (pT_cut_min>pT_cut_max)
        {
            cout << "calculate_flows error: pT_cut_min=" << pT_cut_min << " is bigger than the pT_cut_max=" << pT_cut_max << " in the table." << endl;
            exit(-1);
        }
        if (pT_cut_max>max_pT)
        {
            cout << "calculate_flows error: pT_cut_max=" << pT_cut_max << " is bigger than the maximum pT=" << max_pT << " in the table." << endl;
            exit(-1);
        }
        double gauss_x[GAUSS_N], gauss_w[GAUSS_N];
        GaussLegendre_getWeight(GAUSS_N, gauss_x,gauss_w,pT_cut_min,pT_cut_max);
        Table pT_dN_vn(2+(to_order-from_order+1)*2, pT_tab_length, 0); // 1st col: pT, 2nd col: spectra, 3rd-nth: differential vn (real & imag)
        for (int i=0; i<pT_tab_length; i++)
        {
            int col_id = 1;
            pT_dN_vn.set(col_id, i+1, pT_tab.get(1, i+1)); col_id++;
            if(normalization[i] < 0) normalization[i] = 1e-100; // possible negativity issue
            pT_dN_vn.set(col_id, i+1, log(normalization[i])); col_id++;
            for (int order=from_order; order<=to_order; order++)
            {
                pT_dN_vn.set(col_id, i+1, vn[i][order-from_order][0]); col_id++;
                pT_dN_vn.set(col_id, i+1, vn[i][order-from_order][1]); col_id++;
            }
        }
        Table pT_dN_vn_table(pT_dN_vn);

        for (int i=0; i<GAUSS_N; i++)
        {
            double pT = gauss_x[i], pT_weight = gauss_w[i];

            double normalization_pT = exp(pT_dN_vn_table.interp(1, 2, pT));
            normalizationi += normalization_pT*pT*pT_weight;

            for (int order=from_order; order<=to_order; order++)
            {
                double vn_pT_real = pT_dN_vn_table.interp(1, 2+(order-from_order)*2+1, pT);
                double vn_pT_imag = pT_dN_vn_table.interp(1, 2+(order-from_order)*2+2, pT);
                vni[order-from_order][0] += vn_pT_real*pT*pT_weight;
                vni[order-from_order][1] += vn_pT_imag*pT*pT_weight;
            }
        }


    }


    // store values
    // To mimic:
    // WRITE(30,941) " ", N, " ",X(N)," ",Y(N),
    //   &  " ",X(N)/X(0)," ",Y(N)/X(0)," ",
    //   &  sqrt(X(N)*X(N)+Y(N)*Y(N))/X(0)
    vn_inte.set(1, 1, 0);
    vn_inte.set(2, 1, normalizationi);
    vn_inte.set(3, 1, 0);
    vn_inte.set(4, 1, 1);
    vn_inte.set(5, 1, 0);
    vn_inte.set(6, 1, 1);

    for (int t=0; t<number_of_flows; t++)
    {
    vn_inte.set(1, t+2, from_order+t);
    vn_inte.set(2, t+2, vni[from_order+t-1][0]);
    vn_inte.set(3, t+2, vni[from_order+t-1][1]);
    vn_inte.set(4, t+2, vni[from_order+t-1][0]/normalizationi);
    vn_inte.set(5, t+2, vni[from_order+t-1][1]/normalizationi);
    vn_inte.set(6, t+2, sqrt(vni[from_order+t-1][0]*vni[from_order+t-1][0]+vni[from_order+t-1][1]*vni[from_order+t-1][1])/normalizationi);
    }

    vn_diff.printTable(os_diff);
    vn_inte.printTable(os_inte);

}




//-----------------------------------------------------------------------
void perform_eta_integration(Table* total_N, Table& dN_dy, double mass)
// Calculate N=integrate(dN/dy(eta), eta in eta_table)
{
    Table dN_deta(pT_tab_length, phi_tab_length, 0);
    total_N->setAll(0);
    for (int k=1; k<=eta_table_length; k++)
    {
        double eta = eta_table.get(1, k);
        double eta_weight = eta_table.get(2, k);
        calculate_dN_deta_from_dN_dy_boostInv(&dN_deta, dN_dy, eta, mass);
        for (int i=1; i<=pT_tab_length; i++)
        for (int j=1; j<=phi_tab_length; j++)
        {
            total_N->set(i, j, total_N->get(i,j) + dN_deta.get(i,j)*eta_weight);
        }

    }
}



//-----------------------------------------------------------------------
void calculate_and_output_spectra_and_vn(Table& dN_dy, string particle_name, double mass, double pT_cut_min=-1, double pT_cut_max=-1, int to_order=9, int interp_mode=2)
{
    ofstream of1,of2;
    string vndata_filename = "%s_vndata.dat";
    string vndata_inte_filename = "%s_integrated_vndata.dat";

    char buffer[300];
    sprintf(buffer, vndata_filename.c_str(), particle_name.c_str());
    of1.open(buffer);
    sprintf(buffer, vndata_inte_filename.c_str(), particle_name.c_str());
    of2.open(buffer);

    calculate_flows(dN_dy,of1,of2,mass,to_order,pT_cut_min,pT_cut_max,interp_mode);

    of1.close(); of2.close();
}



//-----------------------------------------------------------------------
int main(int argc, char* argv[])
{

    string particle_name, dN_file;
    if (argc<3)
    {
        printf("Usage: iFlow.e dN_dy_maxtrix_file flow_filename (no ext) [pT_min] [pT_max]\n");
        exit(-1);
    }

    double mass=0, pT_cut_min=-1, pT_cut_max=-1;
    particle_name = argv[2];
    Table pion_dN_dy(argv[1]);
    Table pion_total_N(pT_tab_length, phi_tab_length);
    if (argc>=5)
    {
        pT_cut_min = stringToDouble(argv[3]);
        pT_cut_max = stringToDouble(argv[4]);
    }
    calculate_and_output_spectra_and_vn(pion_dN_dy, particle_name, mass, pT_cut_min, pT_cut_max);

}

