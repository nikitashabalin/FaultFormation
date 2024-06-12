#include <iostream>
#include <algorithm>
#include <string>
#include <math.h>
#include <vector>
#include <opencv2/opencv.hpp>
#include "opencv2/highgui/highgui.hpp"
#include <random>
#include <chrono>

using namespace std;

struct vector2d {
    float x, y;
};

struct vector3d {
    float x, y, z; 
};

struct HeightAndGradient {
    float height; 
    float gradientX; 
    float gradientY; 
};

struct Particle {
    //Construct Particle at Position
    Particle(vector2d _pos) { pos = _pos; }

    vector2d pos;
    vector2d speed = {0.0, 0.0};

    float volume = 1.0;   //This will vary in time
    float sediment = 0.0; //Fraction of Volume that is Sediment!
};

float vectorLength(vector2d speed) {
    return sqrt(speed.x * speed.x + speed.y * speed.y); 
}

vector3d normalizeVector(vector3d vec) {
    float len = sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);

    if (len != 0) {
        vec.x /= len;
        vec.y /= len;
        vec.z /= len;
    }

    return vec;
}

vector3d getNormal(float curX, float curY, vector<std::vector<float>>& map) {

    float scale = 60.0;

    vector3d vec1 = { (map[curX][curY] - map[curX + 1][curY]) * scale, 1.0, 0.0 };    //Positive X
    vector3d vec2 = { (map[curX - 1][curY] - map[curX][curY]) * scale, 1.0, 0.0 };    //Negative X
    vec1 = normalizeVector(vec1);
    vec2 = normalizeVector(vec2);

    vector3d vec3 = { 0.0, 1.0, (map[curX][curY] - map[curX][curY + 1]) * scale };    //Positive Y
    vector3d vec4 = { 0.0, 1.0, (map[curX][curY - 1] - map[curX][curY]) * scale };    //Negative Y
    vec3 = normalizeVector(vec3);
    vec4 = normalizeVector(vec4);

    vector3d n;
    n.x = vec1.x * 0.15; n.y = vec1.y * 0.15; n.z = vec1.z * 0.15;
    n.x += vec2.x * 0.15; n.y += vec2.y * 0.15; n.z += vec2.z * 0.15;
    n.x += vec3.x * 0.15; n.y += vec3.y * 0.15; n.z += vec3.z * 0.15;
    n.x += vec4.x * 0.15; n.y += vec4.y * 0.15; n.z += vec4.z * 0.15;

    //-----------------------------------------FOR DIAGONALS-----------------------------------------------//
    vector3d vecD1 = { scale * (map[curX][curY] - map[curX + 1][curY + 1]) / sqrt(2), sqrt(2), scale * (map[curX][curY] - map[curX + 1][curY + 1]) / sqrt(2) };     //Positive Y
    vector3d vecD2 = { scale * (map[curX][curY] - map[curX + 1][curY - 1]) / sqrt(2), sqrt(2), scale * (map[curX][curY] - map[curX + 1][curY - 1]) / sqrt(2) };     //Positive Y
    vector3d vecD3 = { scale * (map[curX][curY] - map[curX - 1][curY + 1]) / sqrt(2), sqrt(2), scale * (map[curX][curY] - map[curX - 1][curY + 1]) / sqrt(2) };     //Positive Y
    vector3d vecD4 = { scale * (map[curX][curY] - map[curX - 1][curY - 1]) / sqrt(2), sqrt(2), scale * (map[curX][curY] - map[curX - 1][curY - 1]) / sqrt(2) };     //Positive Y
    
    vecD1 = normalizeVector(vecD1);
    vecD2 = normalizeVector(vecD2);
    vecD3 = normalizeVector(vecD3);
    vecD4 = normalizeVector(vecD4); 

    n.x += vecD1.x * 0.1; n.y += vecD1.y * 0.1; n.z += vecD1.z * 0.1;
    n.x += vecD2.x * 0.1; n.y += vecD2.y * 0.1; n.z += vecD2.z * 0.1;
    n.x += vecD3.x * 0.1; n.y += vecD3.y * 0.1; n.z += vecD3.z * 0.1;
    n.x += vecD4.x * 0.1; n.y += vecD4.y * 0.1; n.z += vecD4.z * 0.1;

    return n;
}

int main()
{
    int iterations = 1000;
    int heightMapSize = 1024;
    srand(time(0)); 
    float delta = 0.03f; 
    float maxHeight = 1.0f; 

    

    std::cout << "Enter map size: ";
    std::cin >> heightMapSize;
    std::cout << "Enter number of iterations: ";
    std::cin >> iterations;
    std::cout << "Delta: ";
    std::cin >> delta;
    std::cout << "Max height [-1.0 ; 1.0]: ";
    std::cin >> maxHeight;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, heightMapSize); 

    std::vector<std::vector<float>> heightMap(heightMapSize, std::vector<float>(heightMapSize)); 
    auto start = std::chrono::high_resolution_clock::now();  
    for (int i = 0; i < heightMapSize; ++i) { 
        for (int j = 0; j < heightMapSize; ++j) { 
            heightMap[i][j] = 0.0;
        }
    }
    

    for (int a = 0; a < iterations; a++)
    {

        vector2d point1;
        point1.x = rand() % heightMapSize;
        point1.y = rand() % heightMapSize;

        vector2d point2;
        point2.x = rand() % heightMapSize;
        point2.y = rand() % heightMapSize;

        vector2d vector1;
        vector1.x = point2.x - point1.x;
        vector1.y = point2.y - point1.y;

        for (int y = 0; y < heightMapSize; ++y) {
            for (int x = 0; x < heightMapSize; ++x) {
                vector2d vector2;
                vector2.x = x - point1.x;
                vector2.y = y - point1.y;
                int crossProduct = vector1.x * vector2.y - vector1.y * vector2.x;

                if (crossProduct >= 0) {
                    heightMap[y][x] += delta;
                }   
                else {
                    heightMap[y][x] -= delta;
                }
                
          
            }
        }
    }

    bool a = 0;
    int iterationsNum = 200000;
    std::cout << "Erosion? " << "1 / 0" << std::endl;
    std::cin >> a;


    if (a) {
        std::cout << "Number of iterations: ";
        std::cin >> iterationsNum;
        float minVol = 0.01;
        float density = 1.0;
        float friction = 0.05;
        float depositionRate = 0.1;
        float evapRate = 0.001;
        float dt = 1.4;

        for (int iteration = 0; iteration < iterationsNum; ++iteration) {

            vector2d newpos;
            newpos.x = rand() % (heightMapSize - 1);
            newpos.y = rand() % (heightMapSize - 1);

            Particle drop(newpos);


            while (drop.volume > minVol) {

                vector2d ipos = drop.pos;

                if (ipos.x + 1 < 0 || ipos.x - 1 < 0 || ipos.x + 1 >= heightMapSize - 1 || ipos.y - 1 >= heightMapSize - 1 || ipos.y + 1 < 0 || ipos.y - 1 < 0 || ipos.y + 1 >= heightMapSize - 1 || ipos.y - 1 >= heightMapSize - 1) {
                    break;
                }
                vector3d n = getNormal(ipos.x, ipos.y, heightMap);

                vector2d force = { n.x, n.z };
                drop.speed.x += force.x * dt / (drop.volume * density);
                drop.speed.y += force.y * dt / (drop.volume * density);

                drop.pos.x += dt * drop.speed.x;
                drop.pos.y += dt * drop.speed.y;

                drop.speed.x *= (1.0 - dt * friction);
                drop.speed.y *= (1.0 - dt * friction);

                if (drop.pos.x < 0 || drop.pos.x >= heightMapSize - 1 || drop.pos.y < 0 || drop.pos.y >= heightMapSize - 1) {
                    break;
                }

                float maxSediment = drop.volume * vectorLength(drop.speed) * (heightMap[ipos.x][ipos.y] - heightMap[(int)drop.pos.x][(int)drop.pos.y]);
                if (maxSediment < 0.0) maxSediment = 0.0;
                float sDiff = maxSediment - drop.sediment;

                drop.sediment += dt * depositionRate * sDiff;
                heightMap[ipos.x][ipos.y] -= dt * drop.volume * depositionRate * sDiff;
                drop.volume *= (1.0 - dt * evapRate);

            }
        }
    }

    //---------------------------------------------------------------------------------------------//
    unsigned int* pixels = new unsigned int[heightMapSize * heightMapSize];

    for (int i = 0; i < heightMapSize; ++i) {
        for (int j = 0; j < heightMapSize; ++j) {

            float val = heightMap[i][j];

            if (val > maxHeight)
                val = maxHeight;
            else if (val < -1.0f)
                val = -1.0f;

            int color = ((val + 1.0f) * 32767.5f);
            pixels[i * heightMapSize + j] = color;

        } 
    } 

    heightMap.clear();

    std::string fileName = "FaultFormation";
    std::cout << "Enter file name: ";
    std::cin >> fileName;
    cv::Mat image(heightMapSize, heightMapSize, CV_16UC1, pixels); 
    cv::Mat blurredImage; 
    cv::GaussianBlur(image, blurredImage, cv::Size(35, 35), 0);
    cv::imwrite(fileName+".png", blurredImage);
    cv::imshow("16-bit Image from Array", blurredImage);
    cv::waitKey(0);
    return 0;
}