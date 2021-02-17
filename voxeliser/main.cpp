#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "Point.h"
#include "Rows.h"
#include "VoxelGrid.h"

float signed_volume(const Point &a, const Point &b, const Point &c, const Point &d) {
  // to do
  return 0;
}

bool intersects(const Point &orig, const Point &dest, const Point &v0, const Point &v1, const Point &v2) {
  // to do
  return 0;
}

int main(int argc, const char * argv[]) {
  const char *file_in = "bag_bk.obj";
  const char *file_out = "vox.obj";
  float voxel_size = 1.0;

  // Read file
  std::vector<Point> vertices;
  std::vector<std::vector<unsigned int>> faces;
  // to do
  std:: string  input =  file_in;
  input = "../" + input;
  std::cout << "Reading file: " << file_in << std::endl;
  std::ifstream infile(input.c_str(), std::ifstream::in);
  if (!infile)
  {
    std::cerr << "Input file not found.\n";
    return false;
  }
  std::string cursor;
  std::string line = "";
  std::getline(infile, line);
  while (line != "")
  {
    std::istringstream linestream(line);

    linestream >> cursor;
    double x,y,z;
    int h,j,k;
    if(cursor == "v")
    {
      for(int i = 0; i < 3; i++)
      {
        linestream >> cursor;
        if(i == 0)
        {
          x = std::stod(cursor);
        }
        else if(i==1)
        {
          y = std::stod(cursor);
        }
        else if( i ==2)
        {
          z = std:: stod(cursor);
        }
      }
      auto p = Point(x,y,z);
      vertices.push_back(p);
    } 
    else if (cursor == "f")
    {
      std::vector<unsigned int> fa;
      for(int i = 0; i < 3; i++)
      {
        linestream >> cursor;
        if(i == 0)
        {
          h = std::stoi(cursor);
          fa.push_back(h);
        }
        else if(i==1)
        {
          j = std::stoi(cursor);
          fa.push_back(i);
        }
        else if( i ==2)
        {
          k = std:: stoi(cursor);
          fa.push_back(i);
        }
      }
      faces.push_back(fa);
      
    }
    std::getline(infile, line);
  }
  std::cout << "Number of vertices: " << vertices.size() << std:: endl;
  std::cout << "Number of faces: " << faces.size() << std:: endl;


  // Create grid
  Rows rows;
  // to do
  VoxelGrid voxels(rows.x, rows.y, rows.z);
  
  // Voxelise
  for (auto const &triangle: faces) {
    // to do
  }
  
  // Fill model
  // to do
  
  // Write voxels
  std::ofstream myfile;
  std:: string  output = file_out;
  output = "../" + output;
  myfile.open(output);
  std::cout << "Writing file: " << file_out << std::endl;
  myfile << "output\n";
  myfile.close();
  // to do
  return 0;
}
