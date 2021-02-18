#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include <math.h> 

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

int read_object(std::string  input,   std::vector<Point> &vertices, std::vector<std::vector<unsigned int>> &faces)
{
  std::cout << "Reading file: " << input << std::endl;
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
  return 0;
}

void boundary_box(std::vector<Point> &boundary,std::vector<Point> &vertices )
{
  Point min_pt = vertices[0], max_pt = vertices[0];
  for (int i = 0 ; i < vertices.size(); i++)
  {
    if(vertices[i].x < min_pt.x)
    {
      min_pt.x = vertices[i].x;
    }
    if(vertices[i].y < min_pt.y)
    {
      min_pt.y = vertices[i].y;
    }
    if(vertices[i].z < min_pt.z)
    {
      min_pt.z = vertices[i].z;
    }
    if(vertices[i].x > max_pt.x)
    {
      max_pt.x = vertices[i].x;
    }
    if(vertices[i].y > max_pt.y)
    {
      max_pt.y = vertices[i].y;
    }
    if(vertices[i].z > max_pt.z)
    {
      max_pt.z = vertices[i].z;
    }
  }
  boundary.push_back(min_pt);
  boundary.push_back(max_pt);
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
  //READ FILE
  read_object(input, vertices, faces);
  // Create grid
  std::vector<Point> boundary;
  boundary_box(boundary, vertices);
  std::cout << "Boundary points min :" << boundary[0].x << "," << boundary[0].y << "," << boundary[0].z << std::endl;
  std::cout << "Boundary points max :" << boundary[1].x << "," << boundary[1].y << "," << boundary[1].z << std::endl;
  
  int row_x = std::ceil((boundary[1].x - boundary[0].x)/voxel_size);
  int row_y = std::ceil((boundary[1].y - boundary[0].y)/voxel_size);
  int row_z = std::ceil((boundary[1].z - boundary[0].z)/voxel_size);

  Rows rows(row_x, row_y, row_z);

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
