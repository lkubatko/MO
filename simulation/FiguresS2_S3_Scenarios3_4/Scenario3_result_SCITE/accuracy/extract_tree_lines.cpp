// extract_tree_lines.coo

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_set>

typedef std::vector<std::vector<int>> graph;
graph buildGraph(int number_v, std::vector<std::pair<int, int>> &edges)
{
    graph g(number_v);
    for (auto e : edges)
    {
        g[e.first-1].push_back(e.second-1);
    }
    return g;
}


void read_gv_files(const char* filename, int &number_v, std::vector<std::pair<int, int>> &edges)
{
    std::ifstream file(filename);
    if (file.is_open())
    {
        std::string line;
        getline(file, line);
        getline(file, line);
        std::string delimiter = " -> ";
        std::string delimiter1 = ";";
        
        std::unordered_set<int> all_vertices;
        
        while (getline(file, line))
        {
            if (line[0] == '}')
                break;
            
            int pos = line.find(delimiter);
            std::string token = line.substr(0, pos);
            int first = std::stoi(token);
            line.erase(0, pos+delimiter.length());
            
            pos = line.find(delimiter1);
            token = line.substr(0, pos);
            int second = std::stoi(token);
            
            std::pair<int, int> edge = std::make_pair(first, second);
            edges.push_back(edge);
            
            if (all_vertices.find(first) == all_vertices.end())
            {
                number_v++;
                all_vertices.insert(first);
            }
            if (all_vertices.find(second) == all_vertices.end())
            {
                number_v++;
                all_vertices.insert(second);
            }
        }
        
        file.close();
    }
    else
    {
        std::cout << "Error: Cannot open input file given. " << std::endl;
        return;
    }
}

std::vector<int> computeInDegrees(graph &g)
{
    std::vector<int> indegrees(g.size(), 0);
    for (auto adj : g)
    {
        for (auto v : adj)
        {
            indegrees[v]++;
        }
    }
    return indegrees;
}

void find_all_branches(graph &g, int id, std::vector<int> branch, std::vector<std::vector<int>> &all_branches)
{
    if (id >= g.size())
        return;
    
    if (g[id].empty())
    {
        all_branches.push_back(branch);
        return;
    }
    
    for (auto v : g[id])
    {
        branch.push_back(v);
        find_all_branches(g, v, branch, all_branches);
        branch.pop_back();
    }
}


void output_branch_to_file(const char* filename, int root_id, std::vector<std::vector<int>> &all_branches)
{
    std::ofstream myfile(filename);
    if(myfile.is_open())
    {
        for (auto branch : all_branches)
        {
            myfile<<root_id+1<<" ";
            for (auto v : branch)
            {
                myfile << v+1 << " ";
            }
            myfile<<"\n";
        }
        myfile.close();
    }
    else
    {
        std::cout << "Error: Cannot open output file given. " << std::endl;
    }
}

int main(int argc, char** argv) {
   
   if (argc <= 1)
   {
       std::cout << "No input file given. " << std::endl;
       return 1;
   }
   
   std::vector<std::pair<int, int>> edges;
   int num_v = 0;
   read_gv_files(argv[1], num_v, edges);
 
   graph g = buildGraph(num_v, edges);
   
   std::vector<int> indegrees = computeInDegrees(g);
    
   int root_id;
   for (; root_id < num_v; root_id++)
   {
      if (!indegrees[root_id])
      {
          break;
      }
   }
        
   if (root_id == num_v)
      std::cout << "Root cannot find. " << std::endl;
    
   std::vector<std::vector<int>> all_branches;
   std::vector<int> branch;
   
   find_all_branches(g, root_id, branch, all_branches);   

   output_branch_to_file(argv[2], root_id, all_branches);

   return 0;
}

