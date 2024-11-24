#include <fstream>
#include <sstream>
#include <random>
#include <iostream>
#include <string>
#include <cassert>
#include <iomanip>
#include <climits>

#define SCALE 100

int main(int argc, char *argv[])
{
    if (argc != 4) {
        printf("./parse_data_osm [input file] [data file] [normalized data file]\n");
        exit(-1);
    }
    std::string filename = std::string(argv[1]);

    std::ifstream fin(filename);
    std::string line;

    std::fstream fout(argv[2], std::ios::out);
    uint32_t cnt = 0;

    double glb_low0 = std::numeric_limits<double>::max();
    double glb_low1 = std::numeric_limits<double>::max();
    double glb_high0 = std::numeric_limits<double>::min();
    double glb_high1 = std::numeric_limits<double>::min();
    
    while (std::getline(fin, line)) {

        size_t found_pos = line.find("POLYGON");

        double low0, low1, high0, high1;
        bool first_time = true;

        while (found_pos != std::string::npos) {
            
            line = line.substr(found_pos + 7); // POLYGON has 7 char

            std::istringstream iss2(line);

            std::string polygon_string;
            std::getline(iss2, polygon_string, '(');
            std::getline(iss2, polygon_string, '(');

            std::string location_string;
            std::getline(iss2, location_string, ')');

            std::istringstream iss3(location_string);

            std::string coord_string;

            while (std::getline(iss3, coord_string, ','))
            {
                std::istringstream iss4(coord_string);

                std::string coord;
                double val;

                iss4 >> coord;  
                while (coord[0] == '(') {
                    coord = coord.substr(1);
                }
                val = std::stod(coord.c_str());
                if (first_time) {
                    low0 = val;
                    high0 = val;
                } else {
                    low0 = std::min(low0, val);
                    high0 = std::max(high0, val);
                }
                
                iss4 >> coord;
                while (coord[coord.size() - 1] == ')') {
                    coord = coord.substr(0, coord.size() - 1);
                }
                val = std::stod(coord.c_str());
                if (first_time) {
                    low1 = val;
                    high1 = val;
                    first_time = false;
                } else {
                    low1 = std::min(low1, val);
                    high1 = std::max(high1, val);
                }
            }

            found_pos = line.find("POLYGON");
        }

        fout << "1 " << cnt << " " << std::setprecision(15) << low0 << " " << low1 << " " << high0 << " " << high1 << '\n';
        glb_low0 = std::min(glb_low0, low0);
        glb_low1 = std::min(glb_low1, low1);
        glb_high0 = std::max(glb_high0, high0);
        glb_high1 = std::max(glb_high1, high1);
        cnt ++;
    }
    fin.close();
    fout.close();

    printf("Generate %u objects...\n", cnt);
    printf("Start Normalizing Coordinates...\n");
    printf("glb_low0: %.5f, glb_low1: %.5f, glb_high0: %.5f, glb_high1: %.5f\n", glb_low0, glb_low1, glb_high0, glb_high1);


    double gap0 = glb_high0 - glb_low0;
    double gap1 = glb_high1 - glb_low1;


    uint32_t header, idx;
    double low0, low1, high0, high1;
    fin.open(argv[2]);
    fout.open(argv[3], std::ios::out);
    while (std::getline(fin, line)) {

        std::istringstream iss(line);
        iss >> header >> idx >> low0 >> low1 >> high0 >> high1;
        assert(header == 1);

        low0 = (low0 - glb_low0)/gap0 * SCALE;
        low1 = (low1 - glb_low1)/gap1 * SCALE;
        high0 = (high0 - glb_low0)/gap0 * SCALE;
        high1 = (high1 - glb_low1)/gap1 * SCALE;
        fout << "1 " << idx << " " << std::setprecision(15) << low0 << " " << low1 << " " << high0 << " " << high1 << '\n';
    }
    fin.close();
    fout.close();
}