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
  Point ad = a-d;
  Point bd = b-d;
  Point cd = c-d;
  Point bdcd = bd.cross(cd);
  float dotabcd = ad.dot(bdcd); 
  float vd  = 1.0/6.0 * dotabcd;
  return vd;
}

bool intersects(const Point &orig, const Point &dest, const Point &v0, const Point &v1, const Point &v2) {
  // point v0
  float v01 = signed_volume(orig,dest,v0,v1);
  float v02 = signed_volume(orig,dest,v0,v2);
  bool  sv0 = false;
  if ( v01 * v02 <= 0)
  {
    sv0 = true;
  }
  // point v1
  float v10 = signed_volume(orig,dest,v1,v0);
  float v12 = signed_volume(orig,dest,v1,v2);
  bool  sv1 = false;
  if ( v10 * v12 <= 0)
  {
    sv1 = true;
  }
  // point v2 
  float v20 = signed_volume(orig,dest,v2,v0);
  float v21 = signed_volume(orig,dest,v2,v1);
  bool  sv2 = false;
  if ( v20 * v21 <= 0)
  {
    sv2 = true;
  }

  if( sv0 && sv1 && sv2)
  {
    return true;
  }
  else
  {
    return false;
  }
  ;
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
          fa.push_back(j);
        }
        else if( i ==2)
        {
          k = std:: stoi(cursor);
          fa.push_back(k);
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

bool boundary_interior(Point &bPoint, std::vector<Point> &bTriangle)
{
  bool resp = false;
  if(bTriangle[0].x -1 <= bPoint.x && bTriangle[1].x+1 >=  bPoint.x
  && bTriangle[0].y -1<= bPoint.y && bTriangle[1].y+1  >=  bPoint.y 
  && bTriangle[0].z -1 <= bPoint.z && bTriangle[1].z +1>=  bPoint.z)
  {
    resp = true; 
  }
  return resp;
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
  int l = 0;
  int f =0;
  for (auto const &triangle: faces) {

    // to do
    std::vector<Point> ver_triangles;
    std::vector<Point> bound_triangles;
    //compute the bounding box of the triangle
    for (auto  const &vert: triangle)
    {
      Point p = vertices[vert-1];
      ver_triangles.push_back(p);
    }
    boundary_box(bound_triangles, ver_triangles);
    // get the voxels in the bounding box of the triangle
    for(int k =0 ; k < rows.z; k ++)
    {
      for(  int  j  = 0 ; j < rows.y; j ++)
      {
        for ( int i = 0 ; i < rows.x; i ++)
        {
          Point minC = Point(i*voxel_size + (boundary[0].x), j*voxel_size + (boundary[0].y), k*voxel_size+ (boundary[0].z));
          Point maxC = Point(minC.x+voxel_size, minC.y+voxel_size, minC.z+voxel_size);
          if( boundary_interior(minC, bound_triangles))
          {
            Point f1 = Point(minC.x+0.5, minC.y+0.5, minC.z);
            Point f2 = Point(minC.x+0.5, minC.y+0.5, minC.z+1);
            Point f3 = Point(minC.x+0.5, minC.y, minC.z+0.5);
            Point f4 = Point(minC.x+0.5, minC.y+1, minC.z+0.5);
            Point f5 = Point(minC.x, minC.y+0.5, minC.z+0.5);
            Point f6 = Point(minC.x+1, minC.y+0.5, minC.z+0.5);
            float t11 = signed_volume(ver_triangles[0], ver_triangles[1], ver_triangles[2], f2);
            float t12 = signed_volume(ver_triangles[0], ver_triangles[1], ver_triangles[2], f1);
            float t1 = t11 * t12;
            float t21 = signed_volume(ver_triangles[0], ver_triangles[1], ver_triangles[2], f6);
            float t22 = signed_volume(ver_triangles[0], ver_triangles[1], ver_triangles[2], f5);
            float t2 = t21 * t22;
            float t31 = signed_volume(ver_triangles[0], ver_triangles[1], ver_triangles[2], f4);
            float t32 = signed_volume(ver_triangles[0], ver_triangles[1], ver_triangles[2], f3);
            float t3 = t31 * t32;
            bool t4 = intersects(f2,f1,ver_triangles[0], ver_triangles[1], ver_triangles[2]);
            bool t5 = intersects(f6,f5,ver_triangles[0], ver_triangles[1], ver_triangles[2]);
            bool t6 = intersects(f4,f3,ver_triangles[0], ver_triangles[1], ver_triangles[2]);
            bool rr = false;
            if( (t1 <=0 && t4 ) || (t2 <=0 && t5) || (t3<=0 && t6) )
            {
              voxels(i,j,k) = 1;
              for(int q = 0; q <k ; q ++)
              {
                voxels(i,j,q) = 1;
                if (q !=0 )
                {
                  f++;
                }
              }
              l ++; 
              rr = true;
            }

          }
        }
      }
    }

  }
  std::cout <<  "Number of boxes on the boundary: " <<  l << std::endl; 
  std::cout <<  "Number of boxes inside : " << f << std::endl;   
  // Fill model
  // to do
  
  // Write voxels
  std::ofstream myfile;
  std:: string  output = file_out;
  output = "../" + output;
  myfile.open(output);
  std::cout << "Writing file: " << file_out << std::endl;
  for(int k =0 ; k < rows.z; k ++)
  {
    for(  int  j  = 0 ; j < rows.y; j ++)
    {
      for ( int i = 0 ; i < rows.x; i ++)
      {
        Point pt  = Point(i+boundary[0].x, j+boundary[0].y, k+boundary[0].z);
          myfile << "v "<< pt.x << " " << pt.y << " " << pt.z << "\n";
        
      }
    }
  }
  int t = 0;
  for(int k =0 ; k < rows.z-1; k ++)
  {
    for(  int  j  = 0 ; j < rows.y-1; j ++)
    {
      for ( int i = 0 ; i < rows.x-1; i ++)
      {
        if(voxels(i,j,k) == 1)
        {
          int v1 = i + j*rows.x + k*rows.x*rows.y +1;
          int v2 = i + (j+1)*rows.x + k*rows.x*rows.y +1;
          int v3 = i + 1 + (j+1)*rows.x + k*rows.x*rows.y +1;
          int v4 = i + 1 + j*rows.x + k*rows.x*rows.y  +1;
          int v5 = i + j*rows.x + (k+1)*rows.x*rows.y +1;
          int v6 = i + 1 + j*rows.x + (k+1)*rows.x*rows.y +1;
          int v7 = i + 1 + (j+1)*rows.x + (k+1)*rows.x*rows.y +1;
          int v8 = i + (j+1)*rows.x + (k+1)*rows.x*rows.y +1;
          //std:: cout << v1 << "," << v2 << "," << v3 << "," << v4 << "," << v5 << "," << v6 << "," << v7 << "," << v8 << pt <<std::endl;
          myfile << "f " <<  v1 << " " << v4 << " " << v3 << " " << v2   << "\n";
          myfile << "f " <<  v5 << " " << v6 << " " << v7 << " " << v8   << "\n";
          myfile << "f " <<  v1 << " " << v4 << " " << v6 << " " << v5   << "\n";
          myfile << "f " <<  v2 << " " << v3 << " " << v7 << " " << v8   << "\n";
          myfile << "f " <<  v1 << " " << v2 << " " << v8 << " " << v5   << "\n";
          myfile << "f " <<  v4 << " " << v3 << " " << v7 << " " << v6   << "\n";
          t ++;
        }
      }
    }
  }
  myfile.close();
  std::cout << t << std::endl;
  std::cout << "FINISHED Writing file: " << file_out << std::endl;
  // to do
  return 0;
}
