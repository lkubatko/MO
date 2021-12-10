// extract_in_info.cpp
// g++ -std=c++11 -o extract_in_info extract_in_info.cpp


#include <iostream>
#include <fstream>
#include <string>

void read_res_files(const char* filename)
{
    std::ifstream file(filename);
    std::string outputfile0(filename), outputfile1(filename);
    outputfile0 += ".in.txt";
    outputfile1 += ".par_child.txt";
    std::ofstream output0(outputfile0.c_str());
    std::ofstream output1(outputfile1.c_str());
    
    if (!output0.is_open())
    {
        std::cout << "Error: Cannot open output file 0. " << std::endl;
        return;
    }
    
    if (!output1.is_open())
    {
        std::cout << "Error: Cannot open output file 1. " << std::endl;
        return;
    }
    
    std::ofstream myfile();
    if (file.is_open())
    {
        std::string line;
        getline(file, line);
        
        std::string delimiter = "in";
        std::string delimiter1 = ":";
        std::string delimiter2 = ",";
        std::string delimiter4 = ")";

        int pos, pos1, pos2, pos3;
        std::string token0, token1;
        while ((pos = line.find(delimiter)) != std::string::npos)
        {
	    line.erase(0, pos);
	    pos = line.find(delimiter);
            pos1 = line.find(delimiter1);
            pos2 = line.find(delimiter2);
	    pos3 = line.find(delimiter4);
            if (pos2 == std::string::npos)
                pos2 = pos3;
            else
                pos2 = std::min(pos2, pos3);

            token0 = line.substr(pos, pos1-pos);
            token1 = line.substr(pos1+1, pos2 - pos1 - 1);
            output0 << token0 << "\t" << token1 << std::endl;
            line.erase(0, pos2+1);
        }       

        std::string delimiter3 = "=";
        while (getline(file, line))
        {
            if (line[0] == 'P')
            {
                pos = line.find(delimiter3);
                pos1 = line.find(delimiter2);
                token0 = line.substr(pos+2, pos1-pos-2);
                line.erase(0, pos1+1);
                
                pos = line.find(delimiter3);
                pos1 = line.find(delimiter2);
                token1 = line.substr(pos+2, pos1-pos-2);
                line.erase(0, pos1+1);
                
		while (getline(file, line))
		{
		  if (line[0] == 'c')
		    break;
                  
		  output1 << token0 << "\t" << token1 << "\t" << line << std::endl;      
		}                  
            }
        }
        
        file.close();
        output0.close();
        output1.close();
    }
    else
    {
        std::cout << "Error: Cannot open input file given. " << std::endl;
        return;
    }
}

int main(int argc, char** argv) {
   
    if (argc <= 1)
    {
        std::cout << "No input file given. " << std::endl;
        return 1;
    }
   
    read_res_files(argv[1]);
 
    return 0;
}




