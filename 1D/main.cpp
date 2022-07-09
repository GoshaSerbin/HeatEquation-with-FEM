#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "eigen-3.4.0/Eigen/Eigen"

using namespace Eigen;
using namespace std;

typedef double T;

T M = 1e3;

struct point {
    T x, y, z;

    point(T X = 0, T Y = 0, T Z = 0) {
        x = X;
        y = Y;
        z = Z;
    }
};

point operator*(const T k, point p) {
    p.x *= k;
    p.y *= k;
    p.z *= k;
    return p;
}

point operator+(point p1, point p2) {
    p1.x += p2.x;
    p1.y += p2.y;
    p1.z += p2.z;
    return p1;
}

point operator-(point p1, point p2) {
    p1.x -= p2.x;
    p1.y -= p2.y;
    p1.z -= p2.z;
    return p1;
}

ostream& operator<< (ostream& out, const point& p)
{
    out << p.x << " " << p.y << " " << p.z;
    return out;
}

istream& operator>> (istream& in, point& p)
{
    in >> p.x;
    in >> p.y;
    in >> p.z;
    return in;
}

struct elem {
    int* nodesNum;
    int num;
};

// 2 or 4
int ReadNum(const char* fileName)
{
    int num;
    ifstream ifile;
    ifile.open(fileName);
    if (!ifile.is_open())
    {
        cout << "Error! Could not open file!" << endl;
        return -1;
    }
    ifile >> num;
    ifile.close();
    return num;
}

int ReadRequestFile(const char* fileName, int num, point* coords, int* elemNum, int& type) {

    ifstream ifile;
    ifile.open(fileName);
    if (!ifile.is_open())
    {
        cout << "Error! Could not open file!" << endl;
        return -1;
    }
    ifile >> num;       //2 or 4
    for (int i = 0; i < num; ++i)
    {
        ifile >> coords[i].x;
        ifile >> coords[i].y;
        ifile >> coords[i].z;
    }

    if (num == 2)
        ifile >> elemNum[0];
    else
    {
        ifile >> elemNum[0];
        ifile >> elemNum[1];
    }

    ifile >> type;
    ifile.close();
    return 0;
}

void Mesh2p(const point coord0, const point coord1, int np, point* nodes, int i0) {
    for (int i = 0; i < np; ++i)
    {
        T xi = -1 + 2. * i / (np - 1);//local coords

        T N1 = (1 - xi) / 2.;
        T N2 = (1 + xi) / 2.;
        nodes[i + i0] = N1 * coord0 + N2 * coord1;

        //T Li = i * 1.0 / (np - 1);
        //nodes[i + i0] = coord0 + Li * (coord1 - coord0);
    }
}

void MakeMesh2p(const int elemNum, const point* coords, const int type, int& NE, int& NP, int& NC, point*& nodes, elem*& elems, int**& CPN, int*& CPNsizes)
{

    NE = elemNum;
    if (type == 1) {
        NP = elemNum + 1;
    }
    else {
        NP = 2 * elemNum + 1;
    }

    NC = 2;
    nodes = new point[NP];


    Mesh2p(coords[0], coords[1], NP, nodes, 0);



    elems = new elem[NE];
    if (type == 1) {
        for (int i = 0; i < NE; ++i)
        {
            elems[i].num = 2;
            elems[i].nodesNum = new int[elems[i].num];
            elems[i].nodesNum[0] = i + 1;
            elems[i].nodesNum[1] = i + 2;
        }
    }
    else {
        for (int i = 0; i < NE; ++i)
        {
            elems[i].num = 3;
            elems[i].nodesNum = new int[elems[i].num];
            elems[i].nodesNum[0] = 2 * i + 1;
            elems[i].nodesNum[1] = 2 * i + 2;
            elems[i].nodesNum[2] = 2 * i + 3;
        }
    }

    CPNsizes = new int[2];
    CPNsizes[0] = 1; CPNsizes[1] = 1;

    CPN = new int*[2];
    CPN[0] = new int[1]; CPN[0][0] = 1;
    CPN[1] = new int[1]; CPN[1][0] = NP;
}



int WriteMeshFile(const char* name, const int NE, const int NP, const int NC, point* nodes, elem* elems, int** CPN, int* CPNsizes) {
    ofstream ofile;
    ofile.open(name);
    if (!ofile.is_open())
    {
        cout << "Error! Could not open file!" << endl;
        return -1;
    }

    ofile << NE << " " << NP << " " << NC << endl;
    for (int i = 0; i < NE; ++i)
    {
        ofile << i + 1 << " " << elems[i].num << " ";
        for (int j = 0; j < elems[i].num; ++j) {
            ofile << elems[i].nodesNum[j] << " ";
        }
        ofile << endl;
    }

    for (int i = 0; i < NP; ++i)
    {
        ofile << i + 1 << " " << nodes[i] << endl;
    }



    for (int i = 0; i < NC; ++i) {
        ofile << CPNsizes[i] << " ";
    }
    for (int i = 0; i < NC; ++i) {
        for (int j = 0; j < CPNsizes[i]; ++j) {
            ofile << CPN[i][j] << endl;
        }
    }

    ofile.close();
    return 0;
}


void CreateMeshFile(const char* RequestFileName, const char* MeshFileName) {

    const int num = ReadNum(RequestFileName);
    int type;                           //type = 1
    point* coords = new point[num];
    int elemNum[1];

    ReadRequestFile(RequestFileName, num, coords, elemNum, type);

    int NE, NP, NC;
    point* nodes = nullptr;
    elem* elems = nullptr;
    int** CPN = nullptr; 
    int* CPNsizes = nullptr;
    MakeMesh2p(elemNum[0], coords, type, NE, NP, NC, nodes, elems, CPN, CPNsizes);

    WriteMeshFile(MeshFileName, NE, NP, NC, nodes, elems, CPN, CPNsizes);


    delete[] coords;
    delete[] nodes;
    for (int i = 0; i < NE; ++i) {
        delete[] elems[i].nodesNum;
    }
    delete[] elems;
    for (int i = 0; i < NC; ++i) {
        delete[] CPN[i];
    }
    delete[] CPN;
    delete[] CPNsizes;
}

int ReadMeshFile(const char* fileName, int& NE, int& NP, int& NC, point* &nodes, elem* &elems, int** &CPN, int* &CPNsizes) {
    ifstream ifile;
    ifile.open(fileName);
    if (!ifile.is_open())
    {
        cout << "Error! Could not open file!" << endl;
        return -1;
    }

    ifile >> NE;
    ifile >> NP;
    ifile >> NC;

    nodes = new point[NP];
    elems = new elem[NE];

    for (int i = 0; i < NE; ++i)
    {
        int elemInd;
        ifile >> elemInd;
        elemInd--;
        ifile >> elems[elemInd].num;
        elems[elemInd].nodesNum = new int[elems[elemInd].num];
        for (int j = 0; j < elems[elemInd].num; ++j) {
            ifile >> elems[elemInd].nodesNum[j];
        }
    }

    for (int i = 0; i < NP; ++i)
    {
        int nodeInd;
        ifile >> nodeInd;
        nodeInd--;
        ifile >> nodes[nodeInd];
    }

    CPNsizes = new int[NC];
    for (int i = 0; i < NC; ++i) {
        
        ifile >> CPNsizes[i];
    }

    CPN = new int* [NC];
    for (int i = 0; i < NC; ++i) {
        CPN[i] = new int[CPNsizes[i]];
        for (int j = 0; j < CPNsizes[i]; ++j) {
            ifile >> CPN[i][j];
        }
    }

    ifile.close();
    return 0;
}


struct BoundaryCondition {
    int nodeInd;
    int bcType;
    T alpha;
    T beta;
    T gamma;

};

int ReadBoundaryFile(const char* fileName, BoundaryCondition* &bc, int& bcSize) {
    ifstream ifile;
    ifile.open(fileName);
    if (!ifile.is_open())
    {
        cout << "Error! Could not open file!" << endl;
        return -1;
    }
    ifile >> bcSize;
    bc = new BoundaryCondition[bcSize];
    for (int i = 0; i < bcSize; ++i) {
        ifile >> bc[i].nodeInd; bc[i].nodeInd--;
        ifile >> bc[i].bcType;//1 2 3 
        ifile >> bc[i].alpha;
        ifile >> bc[i].beta;
        ifile >> bc[i].gamma;
        if (bc[i].bcType == 1) {
            bc[i].bcType = 3;
            bc[i].alpha *= M;
            bc[i].beta = 1.0;
            bc[i].gamma *= M;
        }
    }
    ifile.close();
    return 0;
}



T k(T x) {
    return 10;
}

T f(T x) {
    return 0.5+x;
}

void Tabulate(T (*f)(T), point* nodes, const int NP, T* fArray ) {
    for (int i = 0; i < NP; ++i) {
        fArray[i] = f(nodes[i].x);
    }
}


T toLocal(T x, T x_i, T le) {
    return 2. / le * (x - x_i) - 1;
}

T toGlobal(T xi, T x_i, T le) {
    return le / 2. * (xi + 1) + x_i;
}

void Add(elem el, Matrix<double, 2, 2> localMatrix, Vector<double, 2> localRHS, Matrix<double, Dynamic, Dynamic>& systemMatrix, Vector<double, Dynamic>& systemRHS) {

    int i = el.nodesNum[0] - 1;
    int j = el.nodesNum[1] - 1;

    systemMatrix(i, i) += localMatrix(0, 0);
    systemMatrix(i, j) += localMatrix(0, 1);
    systemMatrix(j, i) += localMatrix(1, 0);
    systemMatrix(j, j) += localMatrix(1, 1);

    systemRHS(i) += localRHS(0);
    systemRHS(j) += localRHS(1);

}

//alpha T + beta T' = gamma
void ApplyBoundaryValues(BoundaryCondition* bc, int bcSize, T* kArray, Matrix<double, Dynamic, Dynamic>& systemMatrix, Vector<double, Dynamic>& systemRHS) {

    for (int i = 0; i < bcSize; ++i) {
        int nodeInd = bc[i].nodeInd;
        switch (bc[i].bcType) { //case 1: u' = M*(u-u0)
        case 2:
            systemRHS[nodeInd] += bc[i].gamma / bc[i].beta * kArray[nodeInd];

            break;
        case 3:
            systemRHS[nodeInd] += bc[i].gamma / bc[i].beta * kArray[nodeInd];
            systemMatrix(nodeInd, nodeInd) += bc[i].alpha / bc[i].beta * kArray[nodeInd];

            break;
        }
    }

}

int WriteResultFile(const char* fileName, Vector<double, Dynamic> result) {
    ofstream ofile;
    ofile.open(fileName);
    if (!ofile.is_open())
    {
        cout << "Error! Could not open file!" << endl;
        return -1;
    }
    for (int i = 0; i < result.size(); ++i) {
        ofile <<result(i) << endl;
    }
    

    return 0;
}

void Calculate(elem el, point* nodes, T* kArray, T* fArray, Matrix<double, 2, 2> &localMatrix, Vector<double, 2> &localRHS) {

    int nodeInd0 = el.nodesNum[0] - 1;
    int nodeInd1 = el.nodesNum[1] - 1;
    T le = nodes[nodeInd1].x - nodes[nodeInd0].x;

    T ke = (kArray[nodeInd0] + kArray[nodeInd1]) / 2.0;
    localMatrix = Matrix<double, 2, 2>{ {1, -1}, {-1, 1} };
    localMatrix *= ke / le;
    

    T fe = (fArray[nodeInd0] + fArray[nodeInd1]) / 2.0;
    localRHS = Vector<double, 2>{ 1, 1 };
    localRHS *= le / 2.0 * fe;
}

int main()
{
    const char* reqfilename = "MeshRequest.txt";
    const char* meshfilename = "Mesh.txt";
    const char* bcfilename = "BoundaryConditions.txt";
    const char* ofilename = "Results.txt";

    //CreateMeshFile(reqfilename, meshfilename);

    int NE, NP, NC;
    point* nodes = nullptr;
    elem* elems = nullptr;
    int** CPN = nullptr;
    int* CPNsizes = nullptr;

    ReadMeshFile(meshfilename, NE, NP, NC, nodes, elems, CPN, CPNsizes);   

    T* kArray = new T[NP];
    T* fArray = new T[NP];
    Tabulate(k, nodes, NP, kArray);
    Tabulate(f, nodes, NP, fArray);

    Matrix<double, Dynamic, Dynamic> systemMatrix = MatrixXd::Zero(NP, NP);
    Vector<double, Dynamic> systemRHS = VectorXd::Zero(NP);
    for (int i = 0; i < NE; ++i) {

        Matrix<double, 2, 2> localMatrix;
        Vector<double, 2> localRHS;
        Calculate(elems[i], nodes, kArray,fArray,localMatrix,localRHS);
        Add(elems[i], localMatrix, localRHS, systemMatrix, systemRHS);

    }

    //alpha T + beta T' = gamma

    BoundaryCondition* bc = nullptr;
    int bcSize = 0;
    ReadBoundaryFile(bcfilename, bc, bcSize);

    ApplyBoundaryValues(bc, bcSize, kArray, systemMatrix, systemRHS);

    //solve
    Vector<double, Dynamic> result = systemMatrix.colPivHouseholderQr().solve(systemRHS);

    WriteResultFile(ofilename, result);

    cout << "Done!\n";
    delete[] bc;
    delete[] nodes;
    for (int i = 0; i < NE; ++i) {
        delete[] elems[i].nodesNum;
    }
    delete[] elems;
    for (int i = 0; i < NC; ++i) {
        delete[] CPN[i];
    }
    delete[] CPN;
    delete[] CPNsizes;
    return 0;
}