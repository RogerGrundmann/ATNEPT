/*
 * Neptune General Circulation Modell(ATNEPT)applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and nh3 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to write sequel, transfer and paraview files
*/

#include "cNeptuneModel.h"

using namespace std;
//using namespace AtomUtils;

namespace ParaViewNeptune{
    void dump_array(const string &name, Array &a, double multiplier, ofstream &f) {
        f <<  "    <DataArray type=\"Float32\" Name=\"" << name << "\" format=\"ascii\">\n";
        for (int k = 0; k < a.km; k++){
            for (int j = 0; j < a.jm; j++){
                for (int i = 0; i < a.im; i++){
                    f << (a.x[i][j][k] * multiplier) << endl;
                }
                f << "\n";
            }
            f << "\n";
        }
        f << "\n";
        f << "    </DataArray>\n";
    }
/*
 * 
*/
    void dump_radial(const string &desc, Array &a, double multiplier, int i, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;
        for (int j = 0; j < a.jm; j++){
            for (int k = 0; k < a.km; k++){
                f << (a.x[i][j][k] * multiplier) << endl;
            }
        }
    }
/*
 * 
*/
    void dump_radial_2d(const string &desc, Array_2D &a, double multiplier, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;
        for (int j = 0; j < a.jm; j++){
            for (int k = 0; k < a.km; k++){
                f << (a.y[j][k] * multiplier) << endl;
            }
        }
    }
/*
 * 
*/
    void dump_zonal(const string &desc, Array &a, double multiplier, int k, ofstream &f){
        f <<  "SCALARS " << desc << " float " << 1 << endl;
        f <<  "LOOKUP_TABLE default" << endl;
        for(int i = 0; i < a.im; i++){
            for(int j = 0; j < a.jm; j++){
                f << (a.x[i][j][k] * multiplier) << endl;
            }
        }
    }
/*
 * 
*/
    void dump_longal(const string &desc, Array &a, double multiplier, int j, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;
        for (int i = 0; i < a.im; i++){
            for (int k = 0; k < a.km; k++){
                f << (a.x[i][j][k] * multiplier) << endl;
            }
        }
    }
}
/*
 * 
*/
void cNeptuneModel::paraview_panorama_vts(int n){
    using namespace ParaViewNeptune;
    double x, y, z, dx, dy, dz;
    double r_mix_plus = r_mix * 1e6;
    string Neptune_panorama_vts_File_Name = output_path + "/Neptune_panorama_" 
        + std::to_string(n) + ".vts";
    ofstream Neptune_panorama_vts_File;
    Neptune_panorama_vts_File.precision(4);
    Neptune_panorama_vts_File.setf(ios::fixed);
    Neptune_panorama_vts_File.open(Neptune_panorama_vts_File_Name);
    if(!Neptune_panorama_vts_File.is_open()){
        cerr << "ERROR: could not open shpere_vts file " << __FILE__ 
            << " at line " << __LINE__ << "\n";
        abort();
    }
    Neptune_panorama_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
    Neptune_panorama_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
    Neptune_panorama_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Neptune_panorama_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Neptune_panorama_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Temperature PressureDynamic PressureStatic NH3 NH3Cloud NH3Ice H2O H2OCloud H2OIce Q_Latent Q_Sensible BuoyancyForce \">\n"  << endl;

    Neptune_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_panorama_vts_File << u.x[i][j][k] << " " << v.x[i][j][k] << " " << w.x[i][j][k] << endl;
            }
            Neptune_panorama_vts_File <<  "\n"  << endl;
        }
        Neptune_panorama_vts_File <<  "\n"  << endl;
    }
    Neptune_panorama_vts_File <<  "\n"  << endl;
    Neptune_panorama_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_panorama_vts_File << t.x[i][j][k] * t_ref - 273.15 << endl;
            }
            Neptune_panorama_vts_File <<  "\n"  << endl;
        }
        Neptune_panorama_vts_File <<  "\n"  << endl;
    }
    Neptune_panorama_vts_File <<  "\n"  << endl;
    Neptune_panorama_vts_File <<  "    </DataArray>\n" << endl;
    dump_array("u-component", u, u_0, Neptune_panorama_vts_File);
    dump_array("v-component", v, u_0, Neptune_panorama_vts_File);
    dump_array("w-component", w, u_0, Neptune_panorama_vts_File);

    dump_array("PressureDyn", p_dyn, 1e3, Neptune_panorama_vts_File);
    dump_array("rho_mix", rho_mix, 1.0, Neptune_panorama_vts_File);

    dump_array("CoriolisForce", CoriolisForce, 1.0, Neptune_panorama_vts_File);
//    dump_array("CentrifugalForce", CentrifugalForce, 1e3, Neptune_panorama_vts_File);
    dump_array("BuoyancyForce", BuoyancyForce, 1.0, Neptune_panorama_vts_File);


    dump_array("H2O", h2o, r_mix, Neptune_panorama_vts_File);
    dump_array("H2OCloud", h2o_cloud, r_mix, Neptune_panorama_vts_File);
    dump_array("H2OIce", h2o_ice, r_mix, Neptune_panorama_vts_File);

    dump_array("CH4", ch4, r_mix, Neptune_panorama_vts_File);
    dump_array("CH4Cloud", ch4_cloud, r_mix, Neptune_panorama_vts_File);
    dump_array("CH4Ice", ch4_ice, r_mix, Neptune_panorama_vts_File);

    dump_array("H2S", h2s, r_mix, Neptune_panorama_vts_File);
    dump_array("H2SCloud", h2s_cloud, r_mix, Neptune_panorama_vts_File);
    dump_array("H2SIce", h2s_ice, r_mix, Neptune_panorama_vts_File);
    dump_array("w_h2s", w_h2s, 1.0, Neptune_panorama_vts_File);
//    dump_array("j_h2s", j_h2s, 1.0, Neptune_panorama_vts_File);
//    dump_array("jT_h2s", jT_h2s, 1.0, Neptune_panorama_vts_File);

    dump_array("NH3", nh3, 1.0, Neptune_panorama_vts_File);
    dump_array("NH3Cloud", nh3_cloud, 1.0, Neptune_panorama_vts_File);
    dump_array("NH3Ice", nh3_ice, 1.0, Neptune_panorama_vts_File);
    dump_array("w_nh3", w_nh3, 1.0, Neptune_panorama_vts_File);
//    dump_array("j_nh3", j_nh3, 1.0, Neptune_panorama_vts_File);
//    dump_array("jT_nh3", jT_nh3, 1.0, Neptune_panorama_vts_File);

    dump_array("NH4SH", nh4sh, r_mix_plus, Neptune_panorama_vts_File);
    dump_array("w_nh4sh", w_nh4sh, r_mix_plus, Neptune_panorama_vts_File);
//    dump_array("j_nh4sh", j_nh4sh, 1.0, Neptune_panorama_vts_File);
//    dump_array("jT_nh4sh", jT_nh4sh, 1.0, Neptune_panorama_vts_File);

//    dump_array("Q_Latent", Q_Latent, 1.0, Neptune_panorama_vts_File);
//    dump_array("Q_Sensible", Q_Sensible, 1.0, Neptune_panorama_vts_File);

    Neptune_panorama_vts_File <<  "   </PointData>\n" << endl;
    Neptune_panorama_vts_File <<  "   <Points>\n"  << endl;
    Neptune_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"  << endl;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    dx = 0.1;
    dy = 0.1;
    dz = 0.1;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(k == 0 || j == 0) x = 0.0;
                else x = x + dx;
                Neptune_panorama_vts_File << x << " " << y << " " << z  << endl;
            }
            x = 0;
            y = y + dy;
            Neptune_panorama_vts_File <<  "\n"  << endl;
        }
        y = 0.0;
        z = z + dz;
        Neptune_panorama_vts_File <<  "\n"  << endl;
    }
    Neptune_panorama_vts_File <<  "    </DataArray>\n"  << endl;
    Neptune_panorama_vts_File <<  "   </Points>\n"  << endl;
    Neptune_panorama_vts_File <<  "  </Piece>\n"  << endl;
    Neptune_panorama_vts_File <<  " </StructuredGrid>\n"  << endl;
    Neptune_panorama_vts_File <<  "</VTKFile>\n"  << endl;
    Neptune_panorama_vts_File.close();
    cout << "   File:  " << "Nept_panorama_" 
        << n << ".vts" << "  has been written to Directory:  " 
        << output_path << endl;
    return;
}
/*
 * 
*/
void cNeptuneModel::paraview_vtk_radial(int n, int i_radial){
    using namespace ParaViewNeptune;
    double x, y, z, dx, dy;
    double r_mix_plus = r_mix * 1e6;
    string Neptune_radial_File_Name = output_path + "/Neptune_radial_" 
        + std::to_string(i_radial) + "_" + std::to_string(n) + ".vtk";
    ofstream Neptune_vtk_radial_File;
    Neptune_vtk_radial_File.precision (4);
    Neptune_vtk_radial_File.setf(ios::fixed);
    Neptune_vtk_radial_File.open(Neptune_radial_File_Name);
    if(!Neptune_vtk_radial_File.is_open()){
        cerr << "ERROR: could not open paraview_vtk file " << __FILE__ 
            << " at line " << __LINE__ << "\n";
        abort();
    }
    Neptune_vtk_radial_File <<  "# vtk DataFile Version 3.0" << endl;
    Neptune_vtk_radial_File <<  "Radial_Data_Nept_Circulation\n";
    Neptune_vtk_radial_File <<  "ASCII" << endl;
    Neptune_vtk_radial_File <<  "DATASET STRUCTURED_GRID" << endl;
    Neptune_vtk_radial_File <<  "DIMENSIONS " << km << " "<< jm << " " << 1 << endl;
    Neptune_vtk_radial_File <<  "POINTS " << jm * km << " float" << endl;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    dx = 0.1;
    dy = 0.1;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(k == 0) y = 0.0;
            else y = y + dy;
            Neptune_vtk_radial_File << x << " " << y << " "<< z << endl;
        }
        y = 0.0;
        x = x + dx;
    }
    Neptune_vtk_radial_File <<  "POINT_DATA " << jm * km << endl;
    dump_radial("u-Component", u, u_0, i_radial, Neptune_vtk_radial_File);
    dump_radial("v-Component", v, u_0, i_radial, Neptune_vtk_radial_File);
    dump_radial("w-Component", w, u_0, i_radial, Neptune_vtk_radial_File);
    Neptune_vtk_radial_File <<  "SCALARS Temperature float " << 1 << endl;
    Neptune_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Neptune_vtk_radial_File << t.x[i_radial][j][k] * t_ref - 273.15 << endl;
        }
    }

    dump_radial("thermalmassflux", thermalmassflux, 1.0, i_radial, Neptune_vtk_radial_File);

    dump_radial("H2O", h2o, r_mix, i_radial, Neptune_vtk_radial_File);
    dump_radial("H2OCloud", h2o_cloud, r_mix, i_radial, Neptune_vtk_radial_File);
    dump_radial("H2OIce", h2o_ice, r_mix, i_radial, Neptune_vtk_radial_File);

    dump_radial("CH4", ch4, r_mix, i_radial, Neptune_vtk_radial_File);
    dump_radial("CH4Cloud", ch4_cloud, r_mix, i_radial, Neptune_vtk_radial_File);
    dump_radial("CH4Ice", ch4_ice, r_mix, i_radial, Neptune_vtk_radial_File);

    dump_radial("H2S", h2s, r_mix, i_radial, Neptune_vtk_radial_File);
    dump_radial("H2SCloud", h2s_cloud, r_mix, i_radial, Neptune_vtk_radial_File);
    dump_radial("H2SIce", h2s_ice, r_mix, i_radial, Neptune_vtk_radial_File);
    dump_radial("w_h2s", w_h2s, 1.0, i_radial, Neptune_vtk_radial_File);
    dump_radial("j_h2s", j_h2s, 1.0, i_radial, Neptune_vtk_radial_File);
    dump_radial("jT_h2s", jT_h2s, 1.0, i_radial, Neptune_vtk_radial_File);
    dump_radial("massflux_h2s", massflux_h2s, 1.0, i_radial, Neptune_vtk_radial_File);
    dump_radial("difflux_h2s", difflux_h2s, 1.0, i_radial, Neptune_vtk_radial_File);

    dump_radial("NH3", nh3, r_mix, i_radial, Neptune_vtk_radial_File);
    dump_radial("NH3Cloud", nh3_cloud, r_mix, i_radial, Neptune_vtk_radial_File);
    dump_radial("NH3Ice", nh3_ice, r_mix, i_radial, Neptune_vtk_radial_File);
    dump_radial("w_nh3", w_nh3, 1.0, i_radial, Neptune_vtk_radial_File);
    dump_radial("j_nh3", j_nh3, 1.0, i_radial, Neptune_vtk_radial_File);
    dump_radial("jT_nh3", jT_nh3, 1.0, i_radial, Neptune_vtk_radial_File);
    dump_radial("massflux_nh3", massflux_nh3, 1.0, i_radial, Neptune_vtk_radial_File);
    dump_radial("difflux_nh3", difflux_nh3, 1.0, i_radial, Neptune_vtk_radial_File);

    dump_radial("NH4SH", nh4sh, r_mix_plus, i_radial, Neptune_vtk_radial_File);
    dump_radial("w_nh4sh", w_nh4sh, r_mix_plus, i_radial, Neptune_vtk_radial_File);
    dump_radial("j_nh4sh", j_nh4sh, r_mix_plus, i_radial, Neptune_vtk_radial_File);
    dump_radial("jT_nh4sh", jT_nh4sh, r_mix_plus, i_radial, Neptune_vtk_radial_File);
   dump_radial("massflux_nh4sh", massflux_nh4sh, r_mix_plus, i_radial, Neptune_vtk_radial_File);
    dump_radial("difflux_nh3", difflux_nh4sh, r_mix_plus, i_radial, Neptune_vtk_radial_File);


    dump_radial("PressureDyn", p_dyn, 1e3, i_radial, Neptune_vtk_radial_File);
    dump_radial("PressureStat", p_stat, 1.0, i_radial, Neptune_vtk_radial_File);
    dump_radial("rho_mix", rho_mix, 1.0, i_radial, Neptune_vtk_radial_File);

    dump_radial("CoriolisForce", CoriolisForce, 1.0, i_radial, Neptune_vtk_radial_File);
    dump_radial("CentrifugalForce", CentrifugalForce, 1e3, i_radial, Neptune_vtk_radial_File);
    dump_radial("BuoyancyForce", BuoyancyForce, 1.0, i_radial, Neptune_vtk_radial_File);
    dump_radial("PresGradForce", PresGradForce, 1.0, i_radial, Neptune_vtk_radial_File);

    dump_radial("Q_Latent", Q_Latent, 1.0, i_radial, Neptune_vtk_radial_File);
    dump_radial("Q_Sensible", Q_Sensible, 1.0, i_radial, Neptune_vtk_radial_File);

    Neptune_vtk_radial_File <<  "VECTORS v-w-Cell float " << endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Neptune_vtk_radial_File << v.x[i_radial][j][k] << " " << w.x[i_radial][j][k] << " " << z << endl;
        }
    }
    Neptune_vtk_radial_File.close();
    cout << "   File:  " << "Nept_radial_" 
        << i_radial << "_" << n << ".vtk" 
        << "  has been written to Directory:  " << output_path << endl;
    return;
}
/*
 * 
*/
void cNeptuneModel::paraview_vtk_zonal(int n, int k_zonal){
    using namespace ParaViewNeptune;
    double x, y, z, dx, dy;
    double r_mix_plus = r_mix * 1e6;
    string Neptune_zonal_File_Name = output_path + "/Neptune_zonal_" 
        + std::to_string(k_zonal) + "_" + std::to_string(n) + ".vtk";
    ofstream Neptune_vtk_zonal_File;
    Neptune_vtk_zonal_File.precision(4);
    Neptune_vtk_zonal_File.setf(ios::fixed);
    Neptune_vtk_zonal_File.open(Neptune_zonal_File_Name);
    if(!Neptune_vtk_zonal_File.is_open()){
        cerr << "ERROR: could not open vtk_zonal file " << __FILE__ 
            << " at line " << __LINE__ << "\n";
        abort();
    }
    Neptune_vtk_zonal_File <<  "# vtk DataFile Version 3.0" << endl;
    Neptune_vtk_zonal_File <<  "Zonal_Data_Nept_Circulation\n";
    Neptune_vtk_zonal_File <<  "ASCII" << endl;
    Neptune_vtk_zonal_File <<  "DATASET STRUCTURED_GRID" << endl;
    Neptune_vtk_zonal_File <<  "DIMENSIONS " << jm << " "<< im << " " << 1 << endl;
    Neptune_vtk_zonal_File <<  "POINTS " << im * jm << " float" << endl;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    dx = 0.1;
    dy = 0.05;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            if(j == 0) y = 0.0;
            else y = y + dy;
            Neptune_vtk_zonal_File << x << " " << y << " "<< z << endl;
        }
        y = 0.0;
        x = x + dx;
    }
    Neptune_vtk_zonal_File <<  "POINT_DATA " << im * jm << endl;
    dump_zonal("u-Component", u, u_0, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("v-Component", v, u_0, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("w-Component", w, u_0, k_zonal, Neptune_vtk_zonal_File);
    Neptune_vtk_zonal_File <<  "SCALARS Temperature float " << 1 << endl;
    Neptune_vtk_zonal_File <<  "LOOKUP_TABLE default"  <<endl;

    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Neptune_vtk_zonal_File << t.x[i][j][k_zonal] * t_ref - 273.15 << endl;
            aux.x[i][j][k_zonal] = get_layer_height(i);
        }
    }
    dump_zonal("thermalmassflux", thermalmassflux, 1.0, k_zonal, Neptune_vtk_zonal_File);

    dump_zonal("height", aux, 1.0, k_zonal, Neptune_vtk_zonal_File);

    dump_zonal("H2O", h2o, r_mix, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("H2OCloud", h2o_cloud, r_mix, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("H2OIce", h2o_ice, r_mix, k_zonal, Neptune_vtk_zonal_File);

    dump_zonal("CH4", ch4, r_mix, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("CH4Cloud", ch4_cloud, r_mix, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("CH4Ice", ch4_ice, r_mix, k_zonal, Neptune_vtk_zonal_File);

    dump_zonal("H2S", h2s, r_mix, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("H2SCloud", h2s_cloud, r_mix, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("H2SIce", h2s_ice, r_mix, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("w_h2s", w_h2s, 1.0, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("j_h2s", j_h2s, 1.0, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("jT_h2s", jT_h2s, 1.0, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("massflux_h2s", massflux_h2s, 1.0, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("difflux_h2s", difflux_h2s, 1.0, k_zonal, Neptune_vtk_zonal_File);

    dump_zonal("NH3", nh3, r_mix, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("NH3Cloud", nh3_cloud, r_mix, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("NH3Ice", nh3_ice, r_mix, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("w_nh3", w_nh3, 1.0, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("j_nh3", j_nh3, 1.0, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("jT_nh3", jT_nh3, 1.0, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("massflux_nh3", massflux_nh3, 1.0, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("difflux_nh3", difflux_nh3, 1.0, k_zonal, Neptune_vtk_zonal_File);

    dump_zonal("NH4SH", nh4sh, r_mix_plus, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("w_nh4sh", w_nh4sh, r_mix_plus, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("j_nh4sh", j_nh4sh, r_mix_plus, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("jT_nh4sh", jT_nh4sh, r_mix_plus, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("massflux_nh4sh", massflux_nh4sh, r_mix_plus, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("difflux_nh3", difflux_nh4sh, r_mix_plus, k_zonal, Neptune_vtk_zonal_File);

    dump_zonal("PressureDyn", p_dyn, 1e3, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("PressureStat", p_stat, 1.0, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("rho_mix", rho_mix, 1.0, k_zonal, Neptune_vtk_zonal_File);

    dump_zonal("CoriolisForce", CoriolisForce, 1.0, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("CentrifugalForce", CentrifugalForce, 1e3, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("BuoyancyForce", BuoyancyForce, 1.0, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("PresGradForce", PresGradForce, 1.0, k_zonal, Neptune_vtk_zonal_File);

    dump_zonal("Q_Latent", Q_Latent, 1.0, k_zonal, Neptune_vtk_zonal_File);
    dump_zonal("Q_Sensible", Q_Sensible, 1.0, k_zonal, Neptune_vtk_zonal_File);

    Neptune_vtk_zonal_File <<  "VECTORS u-v-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Neptune_vtk_zonal_File << u.x[i][j][k_zonal] << " " << v.x[i][j][k_zonal] << " " << z << endl;
        }
    }
    Neptune_vtk_zonal_File.close();
    cout << "   File:  " << "Nept_zonal_" 
        << k_zonal << "_" << n << ".vtk" 
        << "  has been written to Directory:  " << output_path << endl;
    return;
}
/*
 * 
*/
void cNeptuneModel::paraview_vtk_longal(int n, int j_longal){
    using namespace ParaViewNeptune;
    double x, y, z, dx, dz;
    double r_mix_plus = r_mix * 1e6;
    string Neptune_longal_File_Name = output_path + "/Neptune_longal_" 
        + std::to_string(j_longal) + "_" + std::to_string(n) + ".vtk";
    ofstream Neptune_vtk_longal_File;
    Neptune_vtk_longal_File.precision(4);
    Neptune_vtk_longal_File.setf(ios::fixed);
    Neptune_vtk_longal_File.open(Neptune_longal_File_Name);
    if(!Neptune_vtk_longal_File.is_open()){
        cerr << "ERROR: could not open vtk_longal file " 
            << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Neptune_vtk_longal_File <<  "# vtk DataFile Version 3.0" << endl;
    Neptune_vtk_longal_File <<  "Longitudinal_Data_Nept_Circulation\n";
    Neptune_vtk_longal_File <<  "ASCII" << endl;
    Neptune_vtk_longal_File <<  "DATASET STRUCTURED_GRID" << endl;
    Neptune_vtk_longal_File <<  "DIMENSIONS " << km << " "<< im << " " << 1 << endl;
    Neptune_vtk_longal_File <<  "POINTS " << im * km << " float" << endl;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    dx = 0.1;
    dz = 0.025;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            if(k == 0){
                z = 0.0;
            }else{
                z = z + dz;
            }
            Neptune_vtk_longal_File << x << " " << y << " "<< z << endl;
        }
        z = 0.0;
        x = x + dx;
    }
    Neptune_vtk_longal_File <<  "POINT_DATA " << im * km << endl;
    dump_longal("u-Component", u, u_0, j_longal, Neptune_vtk_longal_File);
    dump_longal("v-Component", v, u_0, j_longal, Neptune_vtk_longal_File);
    dump_longal("w-Component", w, u_0, j_longal, Neptune_vtk_longal_File);
    Neptune_vtk_longal_File <<  "SCALARS Temperature float " << 1 << endl;
    Neptune_vtk_longal_File <<  "LOOKUP_TABLE default"  <<endl;

    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Neptune_vtk_longal_File << t.x[i][j_longal][k] * t_ref - 273.15 << endl;
            aux.x[i][j_longal][k] = get_layer_height(i);
        }
    }

    dump_longal("thermalmassflux", thermalmassflux, 1.0, j_longal, Neptune_vtk_longal_File);

    dump_longal("height", aux, 1.0, j_longal, Neptune_vtk_longal_File);

    dump_longal("H2O", h2o, r_mix, j_longal, Neptune_vtk_longal_File);
    dump_longal("H2OCloud", h2o_cloud, r_mix, j_longal, Neptune_vtk_longal_File);
    dump_longal("H2OIce", h2o_ice, r_mix, j_longal, Neptune_vtk_longal_File);

    dump_longal("CH4", ch4, r_mix, j_longal, Neptune_vtk_longal_File);
    dump_longal("CH4Cloud", ch4_cloud, r_mix, j_longal, Neptune_vtk_longal_File);
    dump_longal("CH4Ice", ch4_ice, r_mix, j_longal, Neptune_vtk_longal_File);

    dump_longal("H2S", h2s, r_mix, j_longal, Neptune_vtk_longal_File);
    dump_longal("H2SCloud", h2s_cloud, r_mix, j_longal, Neptune_vtk_longal_File);
    dump_longal("H2SIce", h2s_ice, r_mix, j_longal, Neptune_vtk_longal_File);
    dump_longal("w_h2s", w_h2s, 1.0, j_longal, Neptune_vtk_longal_File);
    dump_longal("j_h2s", j_h2s, 1.0, j_longal, Neptune_vtk_longal_File);
    dump_longal("jT_h2s", jT_h2s, 1.0, j_longal, Neptune_vtk_longal_File);
    dump_longal("massflux_h2s", massflux_h2s, 1.0, j_longal, Neptune_vtk_longal_File);
    dump_longal("difflux_h2s", difflux_h2s, 1.0, j_longal, Neptune_vtk_longal_File);

    dump_longal("NH3", nh3, r_mix, j_longal, Neptune_vtk_longal_File);
    dump_longal("NH3Cloud", nh3_cloud, r_mix, j_longal, Neptune_vtk_longal_File);
    dump_longal("NH3Ice", nh3_ice, r_mix, j_longal, Neptune_vtk_longal_File);
    dump_longal("w_nh3", w_nh3, 1.0, j_longal, Neptune_vtk_longal_File);
    dump_longal("j_nh3", j_nh3, 1.0, j_longal, Neptune_vtk_longal_File);
    dump_longal("jT_nh3", jT_nh3, 1.0, j_longal, Neptune_vtk_longal_File);
    dump_longal("massflux_nh3", massflux_nh3, 1.0, j_longal, Neptune_vtk_longal_File);
    dump_longal("difflux_nh3", difflux_nh3, 1.0, j_longal, Neptune_vtk_longal_File);

    dump_longal("NH4SH", nh4sh, r_mix_plus, j_longal, Neptune_vtk_longal_File);
    dump_longal("w_nh4sh", w_nh4sh, r_mix_plus, j_longal, Neptune_vtk_longal_File);
    dump_longal("j_nh4sh", j_nh4sh, r_mix_plus, j_longal, Neptune_vtk_longal_File);
    dump_longal("jT_nh4sh", jT_nh4sh, r_mix_plus, j_longal, Neptune_vtk_longal_File);
    dump_longal("massflux_nh4sh", massflux_nh4sh, r_mix_plus, j_longal, Neptune_vtk_longal_File);
    dump_longal("difflux_nh3", difflux_nh4sh, r_mix_plus, j_longal, Neptune_vtk_longal_File);


    dump_longal("PressureDyn", p_dyn, 1e3, j_longal, Neptune_vtk_longal_File);
    dump_longal("PressureStat", p_stat, 1.0, j_longal, Neptune_vtk_longal_File);
    dump_longal("rho_mix", rho_mix, 1.0, j_longal, Neptune_vtk_longal_File);

    dump_longal("CoriolisForce", CoriolisForce, 1.0, j_longal, Neptune_vtk_longal_File);
    dump_longal("CentrifugalForce", CentrifugalForce, 1e3, j_longal, Neptune_vtk_longal_File);
    dump_longal("BuoyancyForce", BuoyancyForce, 1.0, j_longal, Neptune_vtk_longal_File);
    dump_longal("PresGradForce", PresGradForce, 1.0, j_longal, Neptune_vtk_longal_File);

    dump_longal("Q_Latent", Q_Latent, 1.0, j_longal, Neptune_vtk_longal_File);
    dump_longal("Q_Sensible", Q_Sensible, 1.0, j_longal, Neptune_vtk_longal_File);

    Neptune_vtk_longal_File <<  "VECTORS u-w-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Neptune_vtk_longal_File << u.x[i][j_longal][k] << " " 
                << y << " " << w.x[i][j_longal][k] << endl;
        }
    }
    Neptune_vtk_longal_File.close();
    cout << "   File:  " << "Nept_longal_" 
        << j_longal << "_" << n << ".vtk" 
        << "  has been written to Directory:  " << output_path << endl;
    return;
}
/*
 * 
*/
void cNeptuneModel::paraview_sphere_vts(int n){
    using namespace ParaViewNeptune;
    double x, y, z, sinthe, sinphi, costhe, cosphi;
    double r_mix_plus = r_mix * 1e6;
    string Neptune_sphere_vts_File_Name = output_path + "/Neptune_sphere_" 
        + std::to_string(n) + ".vts";
    ofstream Neptune_sphere_vts_File;
    Neptune_sphere_vts_File.precision(4);
    Neptune_sphere_vts_File.setf(ios::fixed);
    Neptune_sphere_vts_File.open(Neptune_sphere_vts_File_Name);
    if (!Neptune_sphere_vts_File.is_open()){
        cerr << "ERROR: could not open paraview_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Neptune_sphere_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
    Neptune_sphere_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
    Neptune_sphere_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Neptune_sphere_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Neptune_sphere_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Temperature PressureDyn PressureStat  H2S H2SCloud H2SIce NH3 NH3Cloud NH3Ice H2O H2OCloud H2OIce \">\n"  << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        sinphi = sin( phi.z[k]);
        cosphi = cos( phi.z[k]);
        for(int j = 0; j < jm; j++){
            sinthe = sin(the.z[j]);
            costhe = cos(the.z[j]);
            for(int i = 0; i < im; i++){
                aux_u.x[i][j][k] = sinthe * cosphi * u.x[i][j][k] + costhe * cosphi * v.x[i][j][k] - sinphi * w.x[i][j][k];
                aux_v.x[i][j][k] = sinthe * sinphi * u.x[i][j][k] + sinphi * costhe * v.x[i][j][k] + cosphi * w.x[i][j][k];
                aux_w.x[i][j][k] = costhe * u.x[i][j][k] - sinthe * v.x[i][j][k];
                Neptune_sphere_vts_File << aux_u.x[i][j][k] << " " << aux_v.x[i][j][k] << " " << aux_w.x[i][j][k]  << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << t.x[i][j][k] * t_ref - 273.15 << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PressureDyn\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << 1e3 * p_dyn.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
/*
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PressureStat\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << 1e-3 * p_stat.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
*/
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2O\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << h2o.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2S\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << h2s.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH3\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << nh3.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH4SH\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << r_mix_plus * nh4sh.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2OCloud\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << 1e3 * h2o_cloud.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
/*
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2SCloud\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << 1e3 * h2s_cloud.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
*/
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH3Cloud\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << 1e3 * nh3_cloud.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
/*
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << 1e3 * h2.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
*/
    Neptune_sphere_vts_File <<  "\n"  << endl;

    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2OIce\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << 1e3 * h2o_ice.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
/*
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2SIce\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << 1e3 * h2s_ice.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
*/
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH3Ice\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << 1e3 * nh3_ice.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
    Neptune_sphere_vts_File <<  "\n"  << endl;
/*
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"HE\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << 1e3 * he.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
    Neptune_sphere_vts_File <<  "\n"  << endl;
*/
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"u-Component\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << u_0 * aux_u.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"v-Component\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << u_0 * aux_v.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"w-Component\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << u_0 * aux_w.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
/*
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Seamount\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Neptune_sphere_vts_File << SeaMount.x[i][j][k] << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
*/
    Neptune_sphere_vts_File <<  "\n"  << endl;
    Neptune_sphere_vts_File <<  "    </DataArray>\n"  << endl;
    Neptune_sphere_vts_File <<  "   </PointData>\n" << endl;
    Neptune_sphere_vts_File <<  "   <Points>\n"  << endl;
    Neptune_sphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                x = rad.z[i] * sin( the.z[j])* cos(phi.z[k]);
                y = rad.z[i] * sin( the.z[j])* sin(phi.z[k]);
                z = rad.z[i] * cos( the.z[j]);
                Neptune_sphere_vts_File << x << " " << y << " " << z  << endl;
            }
            Neptune_sphere_vts_File <<  "\n"  << endl;
        }
        Neptune_sphere_vts_File <<  "\n"  << endl;
    }
    Neptune_sphere_vts_File <<  "    </DataArray>\n"  << endl;
    Neptune_sphere_vts_File <<  "   </Points>\n"  << endl;
    Neptune_sphere_vts_File <<  "  </Piece>\n"  << endl;
    Neptune_sphere_vts_File <<  " </StructuredGrid>\n"  << endl;
    Neptune_sphere_vts_File <<  "</VTKFile>\n"  << endl;
    Neptune_sphere_vts_File.close();
    cout << "   File:  " << "Nept_sphere_" 
        << n << ".vts" << "  has been written to Directory:  " 
        << output_path << endl;
}
/*
 * 
*/
void cNeptuneModel::NeptunePlotData(){
    string Name_PlotData_File = output_path + "/PlotData_Neptune.xyz";
    ofstream PlotData_File;
    PlotData_File.precision(4);
    PlotData_File.setf(ios::fixed);
    PlotData_File.open(Name_PlotData_File);
    if(!PlotData_File.is_open()){
        cerr << "ERROR: could not open PlotData file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    PlotData_File << "lons(deg)" << ", " << "lats(deg)" << ", " 
        << "topography" << ", " << "v-velocity(m/s)" << ", " 
        << "w-velocity(m/s)" << ", " << "velocity-mag(m/s)" << ", " 
        << "temperature(Celsius)" << ", " << "water_vapour(g/kg)" 
        << ", " << "precipitation(mm)" << ", " 
        <<  "precipitable water(mm)" << endl;
    double vel_mag;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            vel_mag = sqrt(pow(v.x[0][j][k] * u_0, 2) + pow(w.x[0][j][k] * u_0, 2));
//            PlotData_File << k << " " << j << " " << SeaMount.x[0][j][k] << " " 
            PlotData_File << k << " " << j << " "
                << v.x[0][j][k] * u_0 << " " << w.x[0][j][k] * u_0 << " " 
                << vel_mag << " " << t.x[0][j][k] * t_ref - t_ref << " " 
                << h2o.x[0][j][k] << " "<< nh3.x[0][j][k] <<  endl;
        }
    }
    PlotData_File.close();
    return;
}


