#include <iostream>
#include <fstream>
#include <vector>

#include "profile_3d.h"

using namespace std;

int main(int argc, char *argv[])
{
    if(argc != 3)
    {
        cout << "Usage: generate_profile_3d.e participant_list binary_collision_list" << endl;
        exit(0);
    }

    // filename for the participant and binary collision lists
    string participant_filename = argv[1];
    string binary_collision_filename = argv[2];

    cout << "read in files: " << participant_filename 
         << " and " << binary_collision_filename << endl;
    // read in participant list
    vector<participant_info> part_list;
    ifstream part_if(participant_filename.c_str());
    while (!part_if.eof())
    {
        participant_info temp_part;
        part_if >> temp_part.x >> temp_part.y >> temp_part.id;
        part_list.push_back(temp_part);
    }
    part_list.pop_back();

    // read in binary collision list
    vector<participant_info> binary_list;
    ifstream binary_if(binary_collision_filename.c_str());
    while (!binary_if.eof())
    {
        participant_info temp_part;
        binary_if >> temp_part.x >> temp_part.y;
        binary_list.push_back(temp_part);
    }
    binary_list.pop_back();

    cout << "Npart = " << part_list.size() << endl;
    cout << "Ncoll = " << binary_list.size() << endl;
    
    // output grid information
    int grid_nx = 261;
    int grid_ny = 261;
    int grid_neta = 101;
    double grid_dx = 0.1;
    double grid_dy = 0.1;
    double grid_deta = 0.1;

    double ecm = 19.6;
    // flag to turn on fluctuation
    int random_flag = 1;
    
    profile_3d test(part_list, binary_list, grid_nx, grid_ny, grid_neta, 
                    grid_dx, grid_dy, grid_deta, ecm, random_flag);

    test.generate_3d_profile();
    test.output_3d_rhob_profile("test.dat");
    

    return(0);
}
